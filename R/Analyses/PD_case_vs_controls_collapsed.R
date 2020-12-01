###########################################################
#' Script for creating intermediary data 
#' PD case vs controls collapsed measurements
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
############################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
library(doMC)
source("R/utils/populationAnalysisUtils.R")
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

#######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))
registerDoMC(detectCores())


#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "PD_case_vs_controls_collapsed.R"
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   intermediate = get_intermediate_data_ref(),
                   healthcode = get_healthcode_ref())
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path('R/Analyses', SCRIPT_NAME))
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder
MODEL_OUTPUT <- paste0("PD_case_vs_controls_collapsed_measurements_",
                       gsub(" ", "_", get("metadata")$user_group), ".h5")
ANNOTATIONS <- list(
  analysisType = "case vs controls",
  analysisSubtype = "collapsed measurements",
  dataSubtype = "dataMatrix",
  userSubset = get("metadata")$user_group,
  study = get("metadata")$study,
  pipelineStep= "intermediary data")


#######################################################
## Helper
#######################################################
get.required.data <- function(){
  healthCode <- read.csv(
    synGet(SYN_ID_REF$healthcode$case_vs_controls)$path, sep = "\t")
  tapHC  <- healthCode %>% dplyr::filter(activity == "tapping") %>% .$healthCode
  voiHC  <- healthCode %>% dplyr::filter(activity == "voice") %>% .$healthCode
  walkHC <- healthCode %>% dplyr::filter(activity == "walking") %>% .$healthCode
  resHC  <- healthCode %>% dplyr::filter(activity == "resting") %>% .$healthCode
  datTap <- read.csv(synGet(SYN_ID_REF$processed$tapping)$path, sep = "\t")%>% 
    mutate(PD = as.factor(PD))
  datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, sep = "\t")%>% 
    mutate(PD = as.factor(PD))
  datRes <- read.csv(synGet(SYN_ID_REF$processed$resting)$path,sep = "\t")%>% 
    mutate(PD = as.factor(PD))
  datWal <- read.csv(synGet(SYN_ID_REF$processed$walking)$path, sep = "\t")%>% 
    mutate(PD = as.factor(PD))
  tapFeatures <- FEATURE_LIST$tapping
  walkFeatures <- FEATURE_LIST$walking
  restFeatures <- FEATURE_LIST$resting
  voiceFeatures <- FEATURE_LIST$voice
  data.list <- list(
    tap = list(data = datTap, features = tapFeatures, hc = tapHC),
    walk = list(data = datWal, features = walkFeatures, hc = walkHC),
    voice = list(data = datVoi, features = voiceFeatures, hc = voiHC),
    rest = list(data = datRes, features = restFeatures, hc = resHC))
  return(data.list)
}

#' Pipeline for PD case vs controls using collapsed healthcodes
#' @param data activity features (with demographic information of age, education, gender, diagnosis)
#' @param features sensor features to choose from dataframe
#' @param subsample list of sampled healthcode for comparison (matched healthcodes)
#' @param nRuns how many train-test split that will be assessed
PD_case_vs_controls_collapsed_user_analysis_pipeline <- function(data, features, subsample, nRuns){
  ## obtain numeric education and gender variables 
  ## (will be need for the causality tests later, also
  ## glmnet implementation of ridge-regression does not handle factors)
  data$education2 <- NumericEducation(x = data$education)
  data$gender2 <- BinaryGender(x = data$gender)
  
  ## filter out unmatched participants
  datM <- data[data$healthCode %in% subsample,]
  
  ## get collapsed features (matched data)
  auxM <- CollapseFeatures(x = datM, 
                           labelName = "PD", 
                           covNames = c("age", "gender2", "education2"), 
                           subjectIdName = "healthCode", 
                           featNames = features)
  cdatM <- auxM$out
  
  ## get collapsed features (unmatched data) after filtering out
  ## participants with less than 5 records
  data.filtered <- FilterOutParticipantsWithFewRecords(dat = data, thr = 5)
  aux <- CollapseFeatures(x = data.filtered, 
                          labelName = "PD", 
                          covNames = c("age", "gender2", "education2"), 
                          subjectIdName = "healthCode", 
                          featNames = features)
  cdat <- aux$out
  featNames2 <- auxM$cfeatNames ## only sensor features
  featNames3 <- c(featNames2, "age", "gender2", "education2") ## sensor + demographics
  featNames4 <- c("age", "gender2", "education2") ## only demographics features
  
  cat("run matched", "\n")
  matched.analysis <- RunAnalyses(nRuns = nRuns, 
                                  featNames2, 
                                  featNames3, 
                                  featNames4, 
                                  fixedSeed = 123,
                                  dat = cdatM,
                                  respName = "PD", 
                                  nSplits = 2,
                                  negClassName = "FALSE",
                                  posClassName = "TRUE")
  
  cat("run unmatched", "\n")
  unmatched.analysis <- RunAnalyses(nRuns = nRuns, 
                                    featNames2, 
                                    featNames3, 
                                    featNames4, 
                                    fixedSeed = 123,
                                    dat = cdat,
                                    respName = "PD", 
                                    nSplits = 2,
                                    negClassName = "FALSE",
                                    posClassName = "TRUE")
  return(list(unmatched = unmatched.analysis, matched = matched.analysis))
}

main <- function(){
  
  #######################################################
  ## Instantiate hdf5 file for storing correlation matrix
  #######################################################
  unlink(MODEL_OUTPUT)
  h5createFile(MODEL_OUTPUT)
  
  #######################################################
  ## Map Table and Results 
  #######################################################
  data.list <- get.required.data()
  results <- plyr::llply(.data = data.list, 
                         .parallel = TRUE,
                         .fun = function(activity){
                           data <- PD_case_vs_controls_collapsed_user_analysis_pipeline(
                             data = activity$data, 
                             features = activity$features, 
                             subsample = activity$hc, 
                             nRuns = 100)}) %>% 
    purrr::map(names(.), function(x, .){
      build_list_df_to_h5(.[[x]], MODEL_OUTPUT, x)}, .)
  
  #######################################################
  ## Store Results to Synapse
  #######################################################
  f <- synapser::File(MODEL_OUTPUT, OUTPUT_FOLDER_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "PD vs Non PD collapsed measurements (>5 nrecords)",
    executed = GIT_URL, 
    used = c(SYN_ID_REF$processed$tap, 
             SYN_ID_REF$processed$voice, 
             SYN_ID_REF$processed$walk, 
             SYN_ID_REF$processed$rest,
             SYN_ID_REF$healthcode$case_vs_controls)))
  unlink(MODEL_OUTPUT)
}

tryCatch({
  #' create logger for pipeline
  sink('pipeline.log', append = TRUE)
  cat(paste0(
    "[",Sys.time(), "]", " Running ", SCRIPT_NAME), "\n\n")
  sink()
  #' run script
  main()
  #' store logger
  sink('pipeline.log', append = TRUE)
  cat(paste0("[",Sys.time(), "]", " Done Running ", SCRIPT_NAME), "\n\n")
  sink()
}, error = function(e) {
  sink("error.log")
  cat(paste0("[",Sys.time(), "] ", SCRIPT_NAME, " - ", e), "\n\n")
  sink()
  stop("Stopped due to error - Please check error.log")
})
  


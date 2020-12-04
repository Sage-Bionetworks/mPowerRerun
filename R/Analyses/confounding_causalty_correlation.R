###########################################################
#' Script for creating intermediary data 
#' for user confounding causalty correlation in .h5 format
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
library(plyr)
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
SCRIPT_NAME <- "confounding_causalty_correlation.R"
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   intermediate = get_intermediate_data_ref(),
                   healthcode = get_healthcode_ref())
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path('R/Analyses', SCRIPT_NAME))
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder
MODEL_OUTPUT <- paste0("demo_confounders_correlation_",
                       gsub(" ", "_", get("metadata")$user_group), ".h5")
ANNOTATIONS <- list(
  analysisType = "demographics confounders",
  analysisSubtype = "correlation test",
  userSubset = get("metadata")$user_group,
  pipelineStep= "intermediary data")

#######################################################
## Helper
#######################################################
get.required.data <- function(){
  healthCode <- read.csv(synGet(SYN_ID_REF$healthcode$case_vs_controls)$path, sep = "\t")
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


#' pipeline wrapper to run confounding_correlation_pipeline
#' @param data activity features (with demographic information of age, education, gender, diagnosis)
#' @param features sensor features to choose from dataframe
#' @param subsample list of sampled healthcode for comparison (matched healthcodes)
#' @param nRuns how many train-test split that will be assessed
#' @return A list of metrics (log(p-value), correlation) in matrix form for each activity 
confounding_correlation_analysis_pipeline <- function(data, features, subsample, nRuns){
  set.seed(1234567)
  splitSeeds <- sample(seq(1e+4, 1e+5), nRuns, replace = TRUE)
  data$education2 <- NumericEducation(x = data$education)
  data$gender2 <- BinaryGender(x = data$gender)
  data.filtered <- FilterOutParticipantsWithFewRecords(dat = data, thr = 5)
  aux <- CollapseAndShapeData(dat = data.filtered,
                              labelName = "PD", 
                              subjectIdName = "healthCode",
                              covNames = c("age", "gender2", "education2"),
                              featNames = features,
                              negClassName = "FALSE", 
                              posClassName = "TRUE")
  auxM <- CollapseAndShapeData(dat = data[data$healthCode %in% subsample,],
                               labelName = "PD", 
                               subjectIdName = "healthCode",
                               covNames = c("age", "gender2", "education2"),
                               featNames = features,
                               negClassName = "FALSE", 
                               posClassName = "TRUE")
  dat <- aux$dat
  datM <- auxM$dat
  featNamesC <- aux$featNamesC
  age.rf <- RunCondIndepTestsCorRf(dat,
                                   datM,
                                   nRuns,
                                   splitSeeds,
                                   labelName = "PD",
                                   featNames = featNamesC,
                                   confName = "age",
                                   negClassName = "FALSE",
                                   posClassName = "TRUE")
  gender.rf <- RunCondIndepTestsCorRf(dat,
                                      datM,
                                      nRuns,
                                      splitSeeds,
                                      labelName = "PD",
                                      featNames = featNamesC,
                                      confName = "gender2",
                                      negClassName = "FALSE",
                                      posClassName = "TRUE")
  education.rf <- RunCondIndepTestsCorRf(dat,
                                         datM,
                                         nRuns,
                                         splitSeeds,
                                         labelName = "PD",
                                         featNames = featNamesC,
                                         confName = "education2",
                                         negClassName = "FALSE",
                                         posClassName = "TRUE")
  age.rr <- RunCondIndepTestsCorRr(dat,
                                   datM,
                                   nRuns,
                                   splitSeeds,
                                   labelName = "PD",
                                   featNames = featNamesC,
                                   confName = "age",
                                   negClassName = "FALSE",
                                   posClassName = "TRUE")
  gender.rr <- RunCondIndepTestsCorRr(dat,
                                      datM,
                                      nRuns,
                                      splitSeeds,
                                      labelName = "PD",
                                      featNames = featNamesC,
                                      confName = "gender2",
                                      negClassName = "FALSE",
                                      posClassName = "TRUE")
  education.rr <- RunCondIndepTestsCorRr(dat,
                                         datM,
                                         nRuns,
                                         splitSeeds,
                                         labelName = "PD",
                                         featNames = featNamesC,
                                         confName = "education2",
                                         negClassName = "FALSE",
                                         posClassName = "TRUE")
  return(
    list(age_rf = age.rf,
         gender_rf = gender.rf, 
         education_rf = education.rf,
         age_rr = age.rr,
         gender_rr = gender.rr, 
         education_rr = education.rr))
}


main <- function(){

  #######################################################
  ## Instantiate hdf5 file for storing correlation matrix
  #######################################################
  unlink(MODEL_OUTPUT)
  h5createFile(MODEL_OUTPUT)
  
  #######################################################
  ## get data and map data
  #######################################################
  data.list <- get.required.data() 
  results <- plyr::llply(.data = data.list, 
                         .parallel = TRUE,
                         .fun = function(activity){
                           data <- confounding_correlation_analysis_pipeline(
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
    "Run Confounding Causalty Correlation Analysis",
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



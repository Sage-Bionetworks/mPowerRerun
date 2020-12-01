###########################################################
#' Script for creating intermediary data for N of 1 analysis
#' which revolves in assessing each user treatment vs time of day
#' and relative importances of each sensor features used 
#' 
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
############################################################
library(synapser)
library(stringi)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
library(doMC)
library(parallel)
source("R/utils/personalizedAnalysisUtils.R")
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

#######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))
registerDoMC(4)

#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "Nof1_analysis.R"
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   intermediate = get_intermediate_data_ref(),
                   healthcode = get_healthcode_ref())
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = 
                         file.path('R/Analyses', SCRIPT_NAME))
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder
MODEL_OUTPUT <- paste0(
  "N_of_1_assessment_", 
  gsub(" ", "_", get("metadata")$user_group), ".h5")
ANNOTATIONS <- list(
  analysisType = "n of 1 analysis",
  dataSubtype = "dataMatrix",
  userSubset = get("metadata")$user_group,
  study = get("metadata")$study,
  pipelineStep= "intermediary data",
  analysisSubtype = c("treatment vs time-of-day", "relative importances"))
PRE_MEDICATION_THRESH <- 15
POST_MEDICATION_THRESH <- 15
P_VAL_THRESHOLD <- 0.05
TOD_THRESH <- 5

#######################################################
## Helper Functions
#######################################################
#' get all required data
get.required.data <- function(){
  healthCode <- read.csv(
    synGet(SYN_ID_REF$healthcode$n_of_one)$path, sep = "\t")
  tapHC  <- healthCode %>% dplyr::filter(activity == "tapping")
  voiHC  <- healthCode %>% dplyr::filter(activity == "voice")
  walkHC <- healthCode %>% dplyr::filter(activity == "walking")
  resHC  <- healthCode %>% dplyr::filter(activity == "resting")
  datTap <- read.csv(synGet(SYN_ID_REF$processed$tapping)$path, sep = "\t")%>% 
    mutate(PD = as.factor(PD), medTimepoint = as.factor(medTimepoint))
  datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, sep = "\t") %>% 
    mutate(PD = as.factor(PD), medTimepoint = as.factor(medTimepoint))
  datRes <- read.csv(synGet(SYN_ID_REF$processed$resting)$path,sep = "\t") %>% 
    mutate(PD = as.factor(PD), medTimepoint = as.factor(medTimepoint))
  datWal <- read.csv(synGet(SYN_ID_REF$processed$walking)$path, sep = "\t") %>% 
    mutate(PD = as.factor(PD), medTimepoint = as.factor(medTimepoint))
  tapFeatures <- FEATURE_LIST$tapping
  walkFeatures <- FEATURE_LIST$walking
  restFeatures <- FEATURE_LIST$resting
  voiceFeatures <- FEATURE_LIST$voice
  data.list <- list(
    tap = list(data = datTap, features = tapFeatures, tzData = tapHC),
    walk = list(data = datWal, features = walkFeatures, tzData = walkHC),
    voice = list(data = datVoi, features = voiceFeatures, tzData = voiHC),
    rest = list(data = datRes, features = restFeatures, tzData = resHC))
  return(data.list)
}

#' Pipeline for N of 1 Analysis
#' @param data activity features (with demographic information of age, education, gender, diagnosis)
#' @param features sensor features to choose from dataframe
#' @param tzData dataframe for healthcodes with time information and filtered
#'               healthcodes that have information of before-after medication 
#'               of 15 records in threshold
N_of_1_tod_vs_treatment_analysis_pipeline <- function(data, features, tzData){
  ## get the data for N-of-1
  data <- GetDataForNof1(data, PRE_MEDICATION_THRESH, POST_MEDICATION_THRESH) %>%
    IncludeUTCandLocalTimeVariables(., tzData) %>%
    dplyr::filter(tod >= TOD_THRESH) %>%
    GetDataForNof1(., PRE_MEDICATION_THRESH, POST_MEDICATION_THRESH) %>%
    LoessDetrendedFeatures(., features) %>%
    TransformFeatures(dat = ., featNames = features)
  return(list(arima = RunTreatmentVsTodEffectsArima(dat = data, featNames = features), 
              LmNw = RunTreatmentVsTodEffectsLmNw(dat = data, featNames = features),
              ri = RelativeImportances(data, features)))
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
                           data <- N_of_1_tod_vs_treatment_analysis_pipeline(
                             data = activity$data, 
                             features = activity$features, 
                             tzData = activity$tzData)})
  
  results %>% purrr::map(names(.), function(x, .){
    build_list_df_to_h5(.[[x]], MODEL_OUTPUT, x)}, .)
  
  
  #######################################################
  ## Store Results to Synapse
  #######################################################
  f <- synapser::File(MODEL_OUTPUT, OUTPUT_FOLDER_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "Run N of 1 assessment",
    executed = GIT_URL, 
    used = c(SYN_ID_REF$processed$tap, 
             SYN_ID_REF$processed$voice, 
             SYN_ID_REF$processed$walk, 
             SYN_ID_REF$processed$rest,
             SYN_ID_REF$healthcode$n_of_one)))
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

###########################################################
## script for cleaning raw features by merging activity with
## demographics and additional logic like inferring 
## worse hand performance from tapping assessment 
############################################################
library(synapser)
library(dplyr)
library(plyr)
library(purrr)
library(tidyverse)
library(tidyr)
library(lubridate)
library(stringr)
library(doMC)
library(mpowertools)
library(githubr)
library(config)
source("R/utils/projectUtils.R")
source("R/utils/featureEngineeringUtils.R")
source("R/utils/initializeVariables.R")

######################
## Configure
#####################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))

######################
## Global variables
#####################
FEATURE_LIST <- get_features()
SCRIPT_NAME <- "cleanActivityData.R"
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   processed = get_processed_features_ref(),
                   table = get_synapse_table_ref())
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FeatureEngineering", SCRIPT_NAME))
OUTPUT_FILENAME <- list()
OUTPUT_FILENAME$tapping <- paste0(
  "mpowertools_processed_tap_features_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
OUTPUT_FILENAME$resting <- paste0(
  "mpowertools_processed_rest_features_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
OUTPUT_FILENAME$walking <- paste0(
  "mpowertools_processed_walk_features_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
OUTPUT_FILENAME$voice <- paste0(
  "MATLAB_processed_voice_features_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")

DEMOGRAPHICS_LIST <-c(
  "age", "are.caretaker", "deep.brain.stimulation", "diagnosis.year",
  "education", "employment", "gender", "health.history", "healthcare.provider",
  "home.usage", "maritalStatus", "medical.usage", 
  "medical.usage.yesterday", "medication.start.year", "onset.year",
  "packs.per.day", "past.participation", "phone.usage",
  "race", "smartphone", "smoked", "surgery", "video.usage", "years.smoking")

#################
## Helper Functions
#################

#' function to retrieve/clean demographic information
#' of mPower users
get.demo <- function(){
  demo <- as.data.frame(synTableQuery(
    sprintf("SELECT * FROM %s", SYN_ID_REF$table$demo)))
  colnames(demo) <- gsub("_|-", ".", names(demo))
  if("inferred.diagnosis" %in% names(demo)){
    demo <- demo %>% 
      mutate(PD = demo$inferred.diagnosis) %>%
      dplyr::select(-c(professional.diagnosis, inferred.diagnosis)) %>%
      filter(dataGroups %in% c("parkinson", "control", NA))
  }else{
    demo <- demo %>% 
      dplyr::rename("PD" = "professional.diagnosis")
  }
  ## clean demographics data
  demo <- demo %>%
    dplyr::select(healthCode, PD, all_of(DEMOGRAPHICS_LIST)) %>%
    dplyr::filter((!is.infinite(age) & age <= 110) | is.na(age)) %>%
    plyr::ddply(.(healthCode), .fun = function(x){
      x$age = mean(x$age, na.rm = TRUE)
      x$PD = ifelse(
        length(unique(x$PD)) > 1, 
        NA, unique(x$PD))
      return(x[1,])})
  return(demo)
}

#' function to merge required demographic information
#' @param activity.data activity data
#' @param demo demographic data
#' @param features list of features being used 
merge.activity.with.demo <- function(activity.data, demo, features){
  merged.data <- activity.data %>%
    dplyr::select(recordId, 
                  healthCode, 
                  appVersion, 
                  phoneInfo, 
                  createdOn, 
                  all_of(features),
                  medTimepoint) %>% 
    dplyr::left_join(demo, by = c("healthCode"))
  return(merged.data)
}

#' function to clean tap data with
#' inferrence on worse performing hands
clean.tap.data <- function(){
  datTap <- read.delim(
    synGet(SYN_ID_REF$raw$tapping)$path, 
    sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::filter(error != "") %>% 
    dplyr::inner_join(worse.tap.hands.association(.), 
                      by = c("healthCode")) %>%
    dplyr::filter(is.na(tappingHands) | tappingHands == inferred.worse.hand)
  return(datTap)
}

#' function to clean rest data
clean.rest.data <- function(){
  datRes <- read.delim(
    synGet(SYN_ID_REF$raw$resting)$path, 
    sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::filter(error != "")
  return(datRes)
}

#' function to clean walk data
clean.walk.data <- function(){
  emptyPedometer <- as.data.frame(synTableQuery(paste(
    "SELECT recordId FROM", 
    SYN_ID_REF$table$walking, 
    "WHERE \"pedometer_walking_outbound.json.items\" is not null")))$recordId
  datWal <- read.delim(
    synGet(SYN_ID_REF$raw$walking)$path, 
    sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::filter(
      (error_walking_features == "") & (recordId %in% emptyPedometer)) 
  return(datWal)
}

#' function to clean voice data
clean.voice.data <- function(){
  keepRecs <- as.data.frame(synTableQuery(paste(
    "SELECT recordId FROM", 
    SYN_ID_REF$table$voice, 
    "WHERE \"audio_audio.m4a\" is not null")))$recordId
  datVoi <- read.delim(synGet(SYN_ID_REF$raw$voice)$path, 
                       sep = "\t", stringsAsFactors = FALSE) %>% 
    filter(recordId %in% keepRecs) %>%
    filter(Filtering == "frequencyFiltered") %>% 
    dplyr::select(-c("Filtering"))
  return(datVoi)
}

#' function store data to synapse using previous annotation as information
#' and update it with some new annotations
#' @param activity.data activity data
#' @param annotations.list named list of previous data annotation
store.to.synapse <- function(activity.data, 
                             annotations.list,
                             description){
  annotations.list <- annotations.list %>% purrr::map(., .f = function(x){x[[1]]})
  annotations.list$pipelineStep <- "processed"
  activity <- annotations.list$task
  output.file <- OUTPUT_FILENAME[[activity]]
  activity.data %>% write.table(., 
                                output.file, 
                                sep="\t", 
                                row.names=F, quote=F)
  file <- synapser::File(output.file, 
                         parent=SYN_ID_REF$processed$output_folder)
  file$annotations <- annotations.list
  synStore(
    file, activity = Activity(
      description,
      executed = GIT_URL, 
      used = c(SYN_ID_REF$raw[[activity]],
               SYN_ID_REF$table$demo)))
  unlink(output.file)
}
  
main <- function(){
  # get demographics
  demo <- get.demo()
  
  ## clean, merge, and store to synapse
  clean.tap.data() %>% 
    merge.activity.with.demo(
      ., demo, FEATURE_LIST$tapping) %>% 
    store.to.synapse(
      ., synGetAnnotations(SYN_ID_REF$raw$tapping), "merge demo, parse worse hand")
  
  ## clean rest data
  clean.rest.data() %>%
    merge.activity.with.demo(., demo, FEATURE_LIST$resting) %>%
    store.to.synapse(
      ., synGetAnnotations(SYN_ID_REF$raw$resting),"merge demo")
  
  clean.walk.data() %>% 
    merge.activity.with.demo(
      ., demo, FEATURE_LIST$walking) %>% 
    store.to.synapse(
      ., synGetAnnotations(SYN_ID_REF$raw$walking),"merge demo")
  
  clean.voice.data() %>% 
    merge.activity.with.demo(
      ., demo, FEATURE_LIST$voice) %>% 
    store.to.synapse(
      ., synGetAnnotations(SYN_ID_REF$raw$voice), "merge demo")
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

###################################################################
# script for running the cleaning and filtering of n of 1 healthcode
# with an nbefore and nafter medication of 15 records
####################################################################
library(data.table)
library(tidyr)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(synapser)
library(githubr)
library(purrr)
library(config)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

#######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))

#######################################################
## Instantiate Variables and Reference IDs
#######################################################
UTC_DATA_SYN_ID <- 'syn23512320'
SCRIPT_NAME <-  "Nof1FilterHealthcodes.R"
FEATURE_LIST <- get_features()
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   healthcode = get_healthcode_ref())
DEMOGRAPHICS_TBL_SYN_ID <- get("synapse_tables")$demo
WALK_TBL_SYN_ID <- SYN_ID_REF$processed$walking
TAP_TBL_SYN_ID <- SYN_ID_REF$processed$tapping
VOICE_TBL_SYN_ID <- SYN_ID_REF$processed$voice
REST_TBL_SYN_ID <- SYN_ID_REF$processed$resting
OUTPUT_FOLDER <- SYN_ID_REF$healthcode$output_folder
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FilterHealthCode", SCRIPT_NAME))
OUTPUT_SYN_ID <- SYN_ID_REF$healthcode$output_folder
OUTPUT_FILENAME <- paste0(
  "Nof1_filtered_cohort_",
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
ANNOTATIONS <- list(study = get("metadata")$study,
                    userSubset = get("metadata")$user_group,
                    dataSubtype = "processed",
                    analysisType = "n of 1 analysis",
                    pipelineStep = "healthcode subsampling",
                    digitalAssessmentDetails = c("tapping", "walking", "resting", "voice"))

#######################################################
## Helpers
#######################################################

#' get/filter demographics, remove duplicates diagnosis label
get.demo <- function(){
  demo <- as.data.frame(synTableQuery(
    sprintf('SELECT * FROM %s', DEMOGRAPHICS_TBL_SYN_ID)))
  colnames(demo) <- gsub("_|-", ".", names(demo))
  demo <- demo %>%
    dplyr::select(healthCode, 
                  age, gender, 
                  professional.diagnosis, 
                  education) %>%
    dplyr::rename(PD = professional.diagnosis) %>%
    dplyr::mutate(gender = gsub('[^[:alnum:]]','', gender),
                  PD = gsub('[^[:alnum:]]','', PD)) %>%
    dplyr::filter(PD == TRUE)
  return(demo)
}

#' get processed features
get.features <- function(){
  ftrs = map(list(tapping = TAP_TBL_SYN_ID,
                  voice = VOICE_TBL_SYN_ID,
                  walking = WALK_TBL_SYN_ID,
                  resting = REST_TBL_SYN_ID), 
             function(x){
               read.csv(synGet(x)$path, sep = "\t")})
  return(ftrs)
}

#' get UTC informations
get.utc.info <- function(){
  utc.data <-  read.csv(synGet(UTC_DATA_SYN_ID)$path, 
                        sep = "\t") %>% 
    dplyr::select(healthCode, UTC_offset)
  return(utc.data)
}

#' get filtered healthcodes
get.filtered.healthcodes <- function(demo, ftrs, utc.data){
  filtered.hc.med <- plyr::llply(ftrs, .fun = function(x, demo, utc.data){
    # Filter healthCodes with atleast 15 records
    filtered.data = x %>%
      dplyr::group_by(healthCode) %>%
      dplyr::summarise(nbefore = sum(medTimepoint == "Immediately before Parkinson medication", 
                                     na.rm=TRUE),
                       nafter = sum(medTimepoint == "Just after Parkinson medication (at your best)",
                                    na.rm=TRUE)) %>%
      dplyr::inner_join(demo, by = c("healthCode")) %>%
      dplyr::filter(nbefore >= 15 &
                    nafter >= 15) %>%
      dplyr::inner_join(utc.data, by = "healthCode") 
  }, demo, utc.data)
  return(filtered.hc.med)
}

main <- function(){
  
  #######################################################
  ### get required data
  #######################################################
  demo <- get.demo()
  ftrs <- get.features()
  utc.data <- get.utc.info()
  filtered.hc.med <- get.filtered.healthcodes(demo, ftrs, utc.data)

  #######################################################
  ## Save to Synapse
  #######################################################
  ## store to synapse
  purrr::map(names(filtered.hc.med), 
             function(activity){
               filtered.hc.med[[activity]] %>% 
                 dplyr::select(healthCode, nbefore, nafter, 
                               age, gender, education, PD, UTC_offset) %>% 
                 dplyr::mutate(activity = activity)}) %>%
    purrr::reduce(., rbind) %>% 
    write.table(., OUTPUT_FILENAME, sep="\t", row.names=F, quote=F)
  
  f <- synapser::File(OUTPUT_FILENAME, parent = OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(
    f, activity = Activity(
      'filter n of 1 healthcode',
      executed = GIT_URL,
      used = setNames(SYN_ID_REF$processed, NULL) %>% unlist()))
  unlink(OUTPUT_FILENAME)
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




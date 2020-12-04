######################################################
# script for running case controls matching healthcodes
# with criteria of users with at least 5 records
#######################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(MatchIt)
source("R/utils/populationAnalysisUtils.R")
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

############################################
# Configure synapse, github repo, and config
############################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))
set.seed(1234567)

############################################
# Get Variable references
############################################
SYN_ID_REF <- list(healthcode = get_healthcode_ref(),
                  processed = get_processed_features_ref())
FEATURE_LIST <- get_features()
SCRIPT_NAME <- "pdCaseControlsMatching.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FilterHealthCode", SCRIPT_NAME))
OUTPUT_SYN_ID <- SYN_ID_REF$healthcode$output_folder
OUTPUT_FILENAME <- paste0(
  "PD_case_vs_controls_matched_cohort_",
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
ANNOTATIONS <- list(analysisType = "case vs controls",
                    pipelineStep = "healthcode subsampling",
                    userSubset = get("metadata")$user_group)

############################################
## helper function
############################################
get.required.data <- function(){
  ## read dataset
  datTap <- read.csv(synGet(SYN_ID_REF$processed$tap)$path, sep = "\t")
  datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, sep = "\t")
  datRes <- read.csv(synGet(SYN_ID_REF$processed$rest)$path, sep = "\t")
  datWal <- read.csv(synGet(SYN_ID_REF$processed$walk)$path, sep = "\t")
  
  ## get features
  tapFeatures <- FEATURE_LIST$tapping
  walkFeatures <- FEATURE_LIST$walking
  restFeatures <- FEATURE_LIST$resting
  voiceFeatures <- FEATURE_LIST$voice
  
  ## store to named list
  data.list <- list(tapping = list(data = datTap, 
                               features = tapFeatures),
                    walking = list(data = datWal, 
                                features = walkFeatures),
                    voice = list(data = datVoi, 
                                 features = voiceFeatures),
                    resting = list(data = datRes, 
                                features = restFeatures))
  return(data.list)
}


main <- function(){
  ## match each activity
  matched.healthcodes.each.activities <- get.required.data() %>% 
    plyr::llply(., function(activity){
                  data <- PD_case_vs_controls_matching(
                    activity$data %>% tidyr::drop_na(age, gender), 
                    activity$features, 
                    thresh = 5)})
  
  ## store to synapse
  purrr::map(names(matched.healthcodes.each.activities), 
      function(activity){
        matched.healthcodes.each.activities[[activity]] %>% 
          dplyr::select(healthCode, gender, age) %>% 
          dplyr::mutate(activity = activity)}) %>%
    purrr::reduce(., rbind) %>% 
    write.table(., OUTPUT_FILENAME, sep="\t", row.names=F, quote=F)
    
  f <- synapser::File(OUTPUT_FILENAME, parent = OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(
    f, activity = Activity(
    'case vs controls >5 records',
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
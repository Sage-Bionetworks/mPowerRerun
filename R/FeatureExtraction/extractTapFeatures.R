###############################################################################
## Script for extracting raw tap features (no filters)
## with metadata of phoneinfo, appversion, recordId, createdOn, and medTimepoint
###############################################################################
library(synapser)
library(plyr)
library(dplyr)
library(jsonlite)
library(tidyr)
library(lubridate)
library(doMC)
library(mpowertools)
library(githubr)
library(config)
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
TAP_FEATURES <- get_features()$tapping
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   table = get_synapse_table_ref())
TAPPING_ACTIVITY_TABLE_SYN_ID <- SYN_ID_REF$table$tapping
OUTPUT_FOLDER_ID <- SYN_ID_REF$raw$output_folder
SCRIPT_NAME <- "extractTapFeatures.R"
OUTPUT_FILE <- paste0(
  "mpowertools_raw_tap_features_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FeatureExtraction", SCRIPT_NAME))
ANNOTATIONS <- list(task = "tapping", 
                    pipelineStep = "raw",
                    userSubset = get("metadata")$user_group)

### helper function
extract.tap.features <- function(cols, table, features){
  
  #######################################################
  # Query Table from synapse
  #######################################################
  actv_tapping_syntable <- synTableQuery(
    sprintf("select * from %s where %s IS NOT NULL", 
        table, paste0('"', cols, '"')))

  actv_tapping<- as.data.frame(actv_tapping_syntable)
  actv_tapping <- actv_tapping %>%
    dplyr::select(-c(ROW_ID, ROW_VERSION))
  actv_tapping$idx <- rownames(actv_tapping)

  attr(actv_tapping$createdOn, 'tzone') <- 'UTC'
  actv_tapping <- actv_tapping %>%
    mutate( createdOn = as.character(createdOn))

  ##############################################
  # Download JSON Files
  ##############################################
  tappingJsonFiles <- synDownloadTableColumns(actv_tapping_syntable, cols)
  tappingJsonFiles <- data.frame(tapping_json_fileId = names(tappingJsonFiles),
                                 tapping_json_file = as.character(tappingJsonFiles))
  actv_tapping <- base::merge(actv_tapping,tappingJsonFiles,
                              by.x= cols,
                              by.y= "tapping_json_fileId", all=T)

  actv_tapping <- actv_tapping %>%
    mutate(tapping_json_file=as.character(tapping_json_file))

  ######################################
  # Feature Extraction
  ######################################
  tappingFeatures <- ddply(.data=actv_tapping,
                           .variables=c("recordId", 
                                        "appVersion",
                                        "createdOn", 
                                        "healthCode",
                                        "phoneInfo", 
                                        "medTimepoint"),
                           .parallel=TRUE,
                           .fun = function(row) {
                             tryCatch({
                               mpowertools::getTappingFeatures(row$tapping_json_file)
                             }, error=function(err){
                               print(row$idx)
                               stop(paste0('error: ', row$tapping_json_file))
                             })
                           })
  tappingFeatures <- tappingFeatures %>%
    dplyr::select(recordId,
                  healthCode,
                  createdOn,
                  appVersion,
                  phoneInfo,
                  all_of(features),
                  medTimepoint,
                  error) %>%
    mutate_at(features, function(x){as.numeric(x)})
  return(tappingFeatures)
}


main <- function(){
  #######################################################
  ## Map Results
  #######################################################
  tap.cols <- list()
  tap.cols$left <- "tapping_left.json.TappingSamples"
  tap.cols$right <- "tapping_right.json.TappingSamples"
  tap.cols$default <- "tapping_results.json.TappingSamples"
  
  features <- plyr::llply(
    .data = tap.cols, 
    .parallel = FALSE,
    .fun = function(colname){
      tryCatch({
        features <- extract.tap.features(
          colname, TAPPING_ACTIVITY_TABLE_SYN_ID, TAP_FEATURES) 
      }, error = function(e){
        return(data.frame())
      })})
  
  tappingFeatures <- dplyr::bind_rows(
    features$default, 
    features$left %>% 
      mutate(tappingHands = "left"), 
    features$right %>% 
      mutate(tappingHands = "right"))

  #####################################
  # Final Data
  #####################################
  write.table(tappingFeatures, OUTPUT_FILE, sep="\t", row.names=F, quote=F)
  f <- synapser::File(OUTPUT_FILE, OUTPUT_FOLDER_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    executed = GIT_URL, 
    used = c(TAPPING_ACTIVITY_TABLE_SYN_ID)))
  unlink(OUTPUT_FILE)
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







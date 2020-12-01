###############################################################################
## Script for extracting raw rest features (no filters)
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
## Get Credentials and parallel register
#######################################################
synLogin()
config::get()
setGithubToken(
    readLines(get("git")$path))
registerDoMC(detectCores())

#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "extractRestFeatures.R"
FEATURE_LIST <- get_features()
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   table = get_synapse_table_ref())
REST_ACTIVITY_TABLE_SYN_ID <- SYN_ID_REF$table$walking
OUTPUT_FOLDER_ID <- SYN_ID_REF$raw$output_folder
OUTPUT_FILE <- paste0(
    "mpowertools_raw_rest_features_", 
    gsub(" ", "_", get("metadata")$user_group), ".tsv")
GIT_URL <- getPermlink(
    getRepo(get("git")$repo), 
    repositoryPath = file.path('R/FeatureExtraction',  SCRIPT_NAME))
ANNOTATIONS <- list(study = get("metadata")$study,
                    userSubset = get("metadata")$user_group,
                    consortium = "mHealth",
                    pipelineStep = "raw",
                    dataType = "sensor",
                    dataSubtype = "raw",
                    analysisType = "",
                    digitalAssessmentDetails = "resting",
                    digitalAssessmentCategory = "gait",
                    dataCollectionMethod = "active",
                    sensorType = c("accelerometer", "gyroscope"),
                    devicePlatform = "iOS",
                    deviceLocation = c("flat surface", "pocket"),
                    dataAccessInstructions = "syn23277418/wiki/607032")



main <- function(){
    
    #######################################################
    ## Query Table from synapse
    #######################################################
    actv_rest_syntable <- synTableQuery(    
        sprintf("select * from %s", 
                REST_ACTIVITY_TABLE_SYN_ID))
    actv_rest <- as.data.frame(actv_rest_syntable)
    actv_rest <- actv_rest %>% 
        dplyr::select(-c(ROW_ID, ROW_VERSION))
    actv_rest$idx <- rownames(actv_rest)
    attr(actv_rest$createdOn, 'tzone') <- 'UTC'
    actv_rest <- actv_rest %>% mutate( createdOn = as.character(createdOn))
    restFeatures <- FEATURE_LIST$walking
    
    ######################
    # Download JSON Files
    ######################
    # rest JSON files
    rest_json_files <- synDownloadTableColumns(actv_rest_syntable, 
                                               "deviceMotion_walking_rest.json.items")
    rest_json_files <- data.frame(rest_json_fileId =names(rest_json_files),
                                  rest_json_file = as.character(rest_json_files))
    actv_rest <- base::merge(actv_rest,rest_json_files, 
                             by.x="deviceMotion_walking_rest.json.items", 
                             by.y="rest_json_fileId")
    actv_rest <- actv_rest %>% mutate(rest_json_file=as.character(rest_json_file))
    
    actv_rest <- actv_rest %>% dplyr::select(-accel_walking_outbound.json.items,
                                             -pedometer_walking_return.json.items,
                                             -accel_walking_return.json.items,
                                             -deviceMotion_walking_outbound.json.items,
                                             -deviceMotion_walking_return.json.items,
                                             -pedometer_walking_outbound.json.items)
    
    ####################
    # Feature Extraction
    #####################
    #extract Rest features
    restFeatures <- plyr::ddply(.data=actv_rest, 
                          .variables=c("recordId", 
                                       "appVersion",
                                       "createdOn", 
                                       "healthCode",
                                       "phoneInfo", 
                                       "medTimepoint"), 
                          .parallel=TRUE,
                          .fun = function(row) { 
                              tryCatch({ mpowertools::getRestFeatures(row$rest_json_file)},
                                       error = function(err){
                                           print(paste0('Unable to process ', row$rest_json_file))
                                           stop(err) })}) %>% 
        dplyr::select(recordId, healthCode, createdOn, 
                      appVersion, phoneInfo, all_of(FEATURE_LIST$rest), 
                      medTimepoint, error) %>% 
        dplyr::mutate_at(FEATURE_LIST$rest, function(x){as.numeric(x)})
    
    ####################
    # Final Data
    ####################
    write.table(restFeatures, OUTPUT_FILE, sep="\t", row.names=F, quote=F, na="")
    f <- synapser::File(OUTPUT_FILE, OUTPUT_FOLDER_ID)
    f$annotations <- ANNOTATIONS
    synStore(f, activity = Activity(executed = GIT_URL, 
                                    used = c(REST_ACTIVITY_TABLE_SYN_ID)))
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


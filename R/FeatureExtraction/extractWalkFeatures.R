###############################################################################
## Script for extracting raw walk features (no filters)
## with metadata of phoneinfo, appversion, 
## recordId, createdOn, and medTimepoint
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
SCRIPT_NAME <- "extractWalkFeatures.R"
FEATURE_LIST <- get_features()
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   table = get_synapse_table_ref())
WALK_ACTIVITY_TABLE_SYN_ID <- SYN_ID_REF$table$walking
OUTPUT_FOLDER_ID <- SYN_ID_REF$raw$output_folder
OUTPUT_FILE <- paste0(
    "mpowertools_raw_walk_features_", 
    gsub(" ", "_", get("metadata")$user_group), ".tsv")
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FeatureExtraction", SCRIPT_NAME))
ANNOTATIONS <- list(task = "walking", 
                    pipelineStep = "raw",
                    userSubset = get("metadata")$user_group)


main <- function(){
    #######################################################
    ## Query Table from synapse
    #######################################################
    actv_walk_syntable <- synTableQuery(    
        sprintf("select * from %s", WALK_ACTIVITY_TABLE_SYN_ID))
    actv_walk <- as.data.frame(actv_walk_syntable)
    actv_walk <- actv_walk %>% dplyr::select(-c(ROW_ID, ROW_VERSION))
    actv_walk$idx <- rownames(actv_walk)
    attr(actv_walk$createdOn, 'tzone') <- 'UTC'
    actv_walk <- actv_walk %>% 
        mutate(createdOn = as.character(createdOn)) %>% 
        arrange(pedometer_walking_outbound.json.items, 
                deviceMotion_walking_outbound.json.items)
    
    ###################################
    # Download JSON Files
    ###################################
    #download outbound walking json file
    outbound_Walking_json_files <- synDownloadTableColumns(
        actv_walk_syntable, "deviceMotion_walking_outbound.json.items")
    outbound_Walking_json_files <- data.frame(
        outbound_Walking_json_fileId = names(outbound_Walking_json_files),
        outbound_Walking_json_file   = as.character(outbound_Walking_json_files))
    
    #download pedometer data
    outbound_pedometer_json_files <- synDownloadTableColumns(
        actv_walk_syntable, "pedometer_walking_outbound.json.items")
    outbound_pedometer_json_files <- data.frame(
        outbound_pedometer_json_fileId = names(outbound_pedometer_json_files),
        outbound_pedometer_json_file = as.character(outbound_pedometer_json_files))
    
    actv_walk  <- actv_walk %>% 
        dplyr::left_join(., outbound_Walking_json_files, 
                    by=c("deviceMotion_walking_outbound.json.items"="outbound_Walking_json_fileId")) %>%
        dplyr::left_join(., outbound_pedometer_json_files,
                    by = c("pedometer_walking_outbound.json.items"="outbound_pedometer_json_fileId")) %>% 
        mutate(outbound_pedometer_json_file=as.character(outbound_pedometer_json_file)) %>% 
        dplyr::select(-accel_walking_outbound.json.items,
                     -accel_walking_return.json.items,
                     -accel_walking_rest.json.items,
                     -deviceMotion_walking_rest.json.items,
                     -pedometer_walking_return.json.items,
                     -deviceMotion_walking_return.json.items)
    
    actv_walk$outbound_Walking_json_file <- as.character(actv_walk$outbound_Walking_json_file)
    actv_walk$outbound_pedometer_json_file <- as.character(actv_walk$outbound_pedometer_json_file)
    
    #################################
    # Feature Extraction
    #################################
    #extract pedometere features
    pedoFeatures <- ddply(.data=actv_walk, 
                          .variables=c("recordId", 
                                       "appVersion",
                                       "createdOn", 
                                       "healthCode",
                                       "phoneInfo", 
                                       "medTimepoint"), 
                          .parallel=TRUE,
                          .fun = function(row) { 
                              tryCatch({ mpowertools::getPedometerFeatures(row$outbound_pedometer_json_file)},
                                       error = function(err){
                                           print(paste0('Unable to process ', row$outbound_pedometer_json_file))
                                           stop(err) })  
                          })
    pedoFeatures['error_pedometer_features'] = pedoFeatures$error
    pedoFeatures$error <- NULL
    pedoFeatures$outbound_pedometer_json_file <- NULL
    
    
    #extract walking features
    walkFeatures <- ddply(.data=actv_walk, 
                          .variables=c("recordId", 
                                       "appVersion",
                                       "createdOn", 
                                       "healthCode",
                                       "phoneInfo", 
                                       "medTimepoint"), 
                          .parallel=T,
                          .fun = function(row) { 
                              tryCatch({ 
                                  res <- mpowertools::getWalkFeatures(row$outbound_Walking_json_file)
                              },
                              error = function(err){
                                  print(paste0('Unable to process ', row$outbound_Walking_json_file))
                                  stop(err) })  
                          })
    walkFeatures['error_walking_features'] = walkFeatures$error
    walkFeatures$error <- NULL
    walkFeatures$pedoJsonPath <- NULL
    walkFeatures$walkingJsonPath <- NULL
    
    combined.walk.features <- base::merge(pedoFeatures, 
                                walkFeatures, 
                                all.x = T, 
                                all.y = T,
                                sort  = F) %>% 
        dplyr::select(recordId, healthCode, 
                      createdOn, appVersion, phoneInfo,
                      all_of(FEATURE_LIST$walking), 
                      error_walking_features,
                      error_pedometer_features,
                      medTimepoint) %>%
        dplyr::mutate_at(FEATURE_LIST$walk, 
                         function(x){as.numeric(x)})
    
    ##########################
    # Final Data
    ##########################
    write.table(combined.walk.features, 
                OUTPUT_FILE, sep="\t", row.names=F, quote=F, na="")
    f <- synapser::File(OUTPUT_FILE, OUTPUT_FOLDER_ID)
    f$annotations <- list(
        study = get("metadata")$study,
        pipelineStep = "raw",
        dataType = "sensor",
        dataSubtype = "raw",
        digitalAssessmentDetails = "walking",
        digitalAssessmentCategory = "gait",
        dataCollectionMethod = "active",
        sensorType = "accelerometer",
        devicePlatform = "iOS",
        deviceLocation = "pocket")
    synStore(f, activity = Activity(executed = GIT_URL, 
                                    used = c(WALK_ACTIVITY_TABLE_SYN_ID)))
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







###############################################################################
## Script for extracting raw voice features (no filters)
## with metadata of phoneinfo, appversion, 
## recordId, createdOn, and medTimepoint
## Note: This is a script that joins healthcode that you have to a featurized
## table of voice features, as the matlab code is not reproducible.
###############################################################################
library(synapser)
library(dplyr)
library(lubridate)
library(doMC)
library(mpowertools)
library(githubr)
library(config)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

#######################################################
## Get Credentials
#######################################################
synLogin()
config::get()
setGithubToken(
    readLines(get("git")$path))

#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "getVoiceSubset.R"
FEATURE_LIST <- get_features()
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   table = get_synapse_table_ref())
VOICE_ACTIVITY_TABLE_SYN_ID <- SYN_ID_REF$table$voice
VOICE_FEATURES <- get("additional")$voice_features
OUTPUT_FOLDER_ID <- SYN_ID_REF$raw$output_folder
OUTPUT_FILE <- paste0(
    "MATLAB_voice_features_", 
    gsub(" ", "_", get("metadata")$user_group), ".tsv")
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FeatureExtraction", SCRIPT_NAME))
ANNOTATIONS <- list(
    userSubset = get("metadata")$user_group,
    pipelineStep = "raw",
    task = "voice")


main <- function(){
    #######################################################
    ## Merge Table healthcodes with pre-build feature datset
    #######################################################
    datVoi <- read.delim(synGet(
        VOICE_FEATURES)$path, 
        sep = "\t", stringsAsFactors = FALSE)
    
    voi.tbl.med.info <- as.data.frame(synTableQuery(
        paste("SELECT recordId, medTimepoint FROM", 
              VOICE_ACTIVITY_TABLE_SYN_ID)))
    
    datVoi <- datVoi %>% 
        dplyr::select(-c("medTimepoint")) %>%
        dplyr::inner_join(voi.tbl.med.info, 
                          by = c("recordId")) 
    
    ##########################
    # Final Data
    ##########################
    write.table(datVoi, OUTPUT_FILE, sep="\t", row.names=F, quote=F, na="")
    f <- synapser::File(OUTPUT_FILE, OUTPUT_FOLDER_ID)
    f$annotations <- ANNOTATIONS
    synStore(f, activity = Activity(executed = GIT_URL, 
                                    used = c(VOICE_ACTIVITY_TABLE_SYN_ID, VOICE_FEATURES)))
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





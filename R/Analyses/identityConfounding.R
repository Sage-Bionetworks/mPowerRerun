###########################################################
#' Script for running assessment of identity confounding
#' throught running permutation on record-wise splits of the
#' mPower data
#' 
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
############################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(githubr)
library(rhdf5)
library(doMC)
library(plyr)
source("R/utils/identityConfoundingUtils.R")
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

#######################################################
## configure credential and register parallel
#######################################################
synLogin()
config::get()
setGithubToken(
    readLines(get("git")$path))
registerDoMC(detectCores())


#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "identityConfounding.R"
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   intermediate = get_intermediate_data_ref(),
                   healthcode = get_healthcode_ref())
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(getRepo(get("git")$repo),
                       repositoryPath = file.path('R/Analyses', SCRIPT_NAME))
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder
MODEL_OUTPUT <- paste0(
    "identity_confounding_assessment_", 
    gsub(" ", "_", get("metadata")$user_group), ".h5")
ANNOTATIONS <- list(
    analysisType = "identity confounding",
    dataSubtype = "dataMatrix",
    userSubset = get("metadata")$user_group,
    study = get("metadata")$study,
    pipelineStep= "intermediary data")

#######################################################
## helper function
#######################################################

## pipeline to run identity confounding
run_identity_confounding_pipeline <- function(data, features, matched.hc, activity.type){
    nperm <- 1000
    nRuns <- 30
    verbose <- FALSE ## shows progress bar
    parallel <- TRUE ## st out the parallel option
    set.seed(123)
    
    ## read data
    dat <- data %>% 
        dplyr::filter(healthCode %in% matched.hc) %>%
        dplyr::mutate(healthCode = as.factor(healthCode)) %>% 
        dplyr::select(healthCode, all_of(features), PD) %>%
        na.omit(.)
    featNames <- features
    dat$healthCode <- factor(dat$healthCode)
    myseeds <- sample(10000:100000, nRuns, replace = TRUE)
    
    statsRWS <- matrix(NA, nRuns, 4)
    colnames(statsRWS) <- c("auc", "medianDRNull", "permPvalDR", "approxVar")
    
    drRWS <- matrix(NA, nperm, nRuns)
    colnames(drRWS) <- paste("run", seq(nRuns), sep = "")
    
    for (i in seq(nRuns)) {
        cat(activity.type, i, "\n")
        set.seed(myseeds[i])
        recordSplit <- GetIdxTrainTestSplitByRecord(dat, nSplits = 2)
        ####
        cat("compute AUC", "\n")
        aucRWS <- GetAUC(dat = dat, 
                         idxTrain = recordSplit$idxTrain, 
                         idxTest = recordSplit$idxTest, 
                         subjectIdName = "healthCode", 
                         labelName = "PD", 
                         featNames = featNames,
                         negClassName = "FALSE", 
                         posClassName = "TRUE")
        statsRWS[i, "auc"] <- aucRWS$aucObs
        statsRWS[i, "approxVar"] <- aucRWS$approxVar["v"]
        cat("run permutation tests", "\n")
        drRWS[, i] <- DRPermDistrAUC(dat = dat, 
                                     idxTrain = recordSplit$idxTrain, 
                                     idxTest = recordSplit$idxTest, 
                                     nperm = nperm, 
                                     subjectIdName = "healthCode", 
                                     labelName = "PD", 
                                     featNames = featNames,
                                     negClassName = "FALSE", 
                                     posClassName = "TRUE",
                                     verbose = verbose,
                                     parallel = parallel)
        statsRWS[i, "medianDRNull"] <- median(drRWS[, i], na.rm = TRUE)
        statsRWS[i, "permPvalDR"] <- sum(drRWS[, i] >= statsRWS[i, "auc"])/nperm
    }
    return(list(statsRWS = statsRWS %>% data.frame(.), 
                drRWS = drRWS %>% data.frame(.)))
}

get.required.data <- function(){
    healthCode <- read.csv(
        synGet(SYN_ID_REF$healthcode$identity_confounding)$path, sep = "\t")
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

###########################################################################
## Instantiate Required Matched healthcode, feature data, list of features
###########################################################################
main <- function(){
    ## process confounding analysis
    data.list <- get.required.data() 
    process.confounding <- data.list %>%
        plyr::llply(.data = ., 
                    .parallel = FALSE,
                    .fun = function(activity){
                        data <- run_identity_confounding_pipeline(
                            data = activity$data, 
                            features = activity$features, 
                            matched.hc = activity$hc,
                            activity.type = activity$type)}) %>% 
        purrr::map(names(.), function(x, .){
            build_list_df_to_h5(.[[x]], MODEL_OUTPUT, x)}, .)
    
    ## store to synapse
    f <- synapser::File(MODEL_OUTPUT, OUTPUT_FOLDER_ID)
    f$annotations <- ANNOTATIONS
    synStore(f, activity = Activity(
        "process identity confounding analysis for figures",
        executed = GIT_URL, 
        used = c(
            SYN_ID_REF$processed$tap, 
            SYN_ID_REF$processed$voice, 
            SYN_ID_REF$processed$walk, 
            SYN_ID_REF$processed$rest,
            SYN_ID_REF$healthcode$identity_confounding)))
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




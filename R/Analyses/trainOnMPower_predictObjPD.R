#######################################################
#' This script will rerun the prediction of ObjectivePD
#' using the model trained on mPower, it will give out
#' model results for each user in .tsv file and serialized model
#' used from trainin on the mPower user-collapsed data
#'  
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
#######################################################

#######################################################
## Library Imports
#######################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
library(config)
library(randomForest)
library(plyr)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

#######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(readLines(
    file.path(path.expand("~"), "git_token.txt")))
set.seed(123456789)

#######################################################
## Global Variables
#######################################################
SCRIPT_NAME <- "trainOnMPower_predictObjPD.R"
GIT_URL <- getPermlink(
    getRepo(get("git")$repo,
            ref="branch", 
            refName=get("git")$branch), 
    repositoryPath = file.path('R/Analyses', SCRIPT_NAME))
SYN_ID_REF <- list(
    processed = get_processed_features_ref(),
    healthcode = get_healthcode_ref(),
    intermediate = get_intermediate_data_ref(),
    objectivePD = get_obj_pd_ref())
FEATURE_LIST <- get_features()
OUTPUT_FILENAME <- paste0("PD_Probabilities_ObjectivePD_Users_",
                          gsub(" ", "_", get("metadata")$user_group), ".tsv")
ANNOTATIONS <- list(
    analysisType = "combined model",
    userSubset = get("metadata")$user_group,
    pipelineStep= "intermediary data",
    analysisSubtype = "confidence score on ObjectivePD users")


#######################################################
## Helpers
#######################################################
## Get collapsed features on the ObjectivePD data.
GetCollapsedFeaturesOPD <- function(dat, respName, featNames, covNames = NULL, negClassName, posClassName) {
    ids <- na.omit(as.character(unique(dat$healthCode))) 
    nids <- length(ids)
    nfeat <- length(featNames)
    out <- data.frame(matrix(NA, nids, 2*nfeat + 2))
    colnames(out) <- c("healthCode", respName, paste(featNames, "md", sep = "."), paste(featNames, "iqr", sep = "."))
    rownames(out) <- ids
    for (i in seq(nids)) {
        #cat(i, "\n")
        sdat <- dat[which(dat$healthCode == ids[i]),]
        out[i, "healthCode"] <- ids[i]
        out[i, respName] <- as.character(sdat[1, respName])
        out[i, 3:(nfeat+2)] <- apply(sdat[, featNames], 2, median, na.rm = TRUE)
        out[i, (nfeat+3):(2*nfeat+2)] <- apply(sdat[, featNames], 2, IQR, na.rm = TRUE)
    }
    out[, respName] <- factor(out[, respName])
    if (!is.null(covNames)) {
        covDat <- dat[!duplicated(dat$healthCode), c("healthCode", covNames)]
        out <- data.frame(out, covDat[match(out$healthCode, covDat$healthCode),])
    }
    return(out)
}

## Get the mPower data 
## (features + PD + medTimepoint + selected covariates)
GetmPowerData <- function(activity, covNames) {
    healthCode <- read.csv(
        synGet(SYN_ID_REF$healthcode$case_vs_controls)$path, sep = "\t")
    if (activity == "tap") {
        tapHC  <- healthCode %>% 
            dplyr::filter(activity == "tapping") %>% .$healthCode
        featNames <- FEATURE_LIST$tapping
        datTap <- read.delim(synGet(SYN_ID_REF$processed$tapping)$path, 
                             sep = "\t", stringsAsFactors = TRUE) %>% 
            dplyr::mutate(PD = as.factor(PD))
        datTap <- datTap[datTap$healthCode %in% tapHC,]
        dat <- datTap[, c("healthCode", featNames, "PD", "medTimepoint", covNames)]
        dat$healthCode <- factor(dat$healthCode)
        dat <- na.omit(dat)
    }
    if (activity == "voi") {
        voiHC  <- healthCode %>% 
            dplyr::filter(activity == "voice") %>% .$healthCode
        featNames <- FEATURE_LIST$voice
        datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, 
                           sep = "\t", stringsAsFactors = TRUE) %>% 
            dplyr::mutate(PD = as.factor(PD))
        datVoi <- datVoi[datVoi$healthCode %in% voiHC,]
        dat <- datVoi[, c("healthCode", featNames, "PD", "medTimepoint", covNames)]
        dat$healthCode <- factor(dat$healthCode)
        dat <- na.omit(dat)
    }
    if (activity == "wal") {
        walHC <- healthCode %>% 
            dplyr::filter(activity == "walking") %>% .$healthCode
        featNames <- FEATURE_LIST$walking
        datWal <- read.csv(synGet(SYN_ID_REF$processed$walking)$path, 
                           sep = "\t", stringsAsFactors = TRUE)%>% 
            dplyr::mutate(PD = as.factor(PD))
        datWal <- datWal[datWal$healthCode %in% walHC,]
        dat <- datWal[, c("healthCode", featNames, "PD", "medTimepoint", covNames)]
        dat$healthCode <- factor(dat$healthCode)
        dat <- na.omit(dat)
    }
    if (activity == "res") {
        resHC  <- healthCode %>% 
            dplyr::filter(activity == "resting") %>% .$healthCode
        featNames <- FEATURE_LIST$resting
        datRes <- read.csv(synGet(SYN_ID_REF$processed$resting)$path,
                           sep = "\t", stringsAsFactors = TRUE)%>% 
            dplyr::mutate(PD = as.factor(PD))
        datRes <- datRes[datRes$healthCode %in% resHC,]
        dat <- datRes[, c("healthCode", featNames, "PD", "medTimepoint", covNames)]
        dat$healthCode <- factor(dat$healthCode)
        dat <- na.omit(dat)
    }
    return(list(dat = dat, featNames = featNames))
}


## Collapse the mPower data (i.e., for each healthCode get the median and IQR across all records).
##
GetCollapsedFeaturesmPower <- function(dat, respName, featNames, covNames = NULL, negClassName, posClassName) {
    ids <- na.omit(as.character(unique(dat$healthCode))) 
    nids <- length(ids)
    nfeat <- length(featNames)
    out <- data.frame(matrix(NA, nids, 2*nfeat + 2))
    colnames(out) <- c("healthCode", respName, paste(featNames, "md", sep = "."), paste(featNames, "iqr", sep = "."))
    rownames(out) <- ids
    for (i in seq(nids)) {
        #cat(i, "\n")
        sdat <- dat[which(dat$healthCode == ids[i]),]
        out[i, "healthCode"] <- ids[i]
        out[i, respName] <- as.character(sdat[1, respName])
        out[i, 3:(nfeat+2)] <- apply(sdat[, featNames], 2, median, na.rm = TRUE)
        out[i, (nfeat+3):(2*nfeat+2)] <- apply(sdat[, featNames], 2, IQR, na.rm = TRUE)
    }
    out[, respName] <- factor(out[, respName])
    if (!is.null(covNames)) {
        covDat <- dat[!duplicated(dat$healthCode), c("healthCode", covNames)]
        out <- data.frame(out, covDat[match(out$healthCode, covDat$healthCode),])
    }
    return(out)
}

FilterOutParticipantsWithFewRecords <- function(dat, thr) {
    aux <- table(dat$healthCode)
    aux <- sort(aux, decreasing = TRUE)
    keep <- names(aux)[which(aux >= thr)]
    return(dat[dat$healthCode %in% keep,])
}


## For an especified activity, load the mPower data
## train a random forest on it and then predict on the
## ObjectivePD data.
GetRfPrediction <- function(activity, odatC, opdHC) {
    aux <- GetmPowerData(activity, covNames = NULL)
    dat <- aux$dat
    hc <- unique(dat$healthCode)
    ## filter out healthCodes that are also in ObjectivePD 
    hc <- setdiff(hc, opdHC) 
    dat <- dat[dat$healthCode %in% hc,]
    featNames <- aux$featNames
    featNames2 <- c(paste(aux$featNames, "md", sep = "."), 
                    paste(aux$featNames, "iqr", sep = "."))
    
    dat <- FilterOutParticipantsWithFewRecords(dat = dat, thr = 5)
    cdat <- GetCollapsedFeaturesmPower(dat = dat, 
                                       respName = "PD", 
                                       featNames = featNames, 
                                       covNames = NULL, 
                                       negClassName = FALSE, 
                                       posClassName = TRUE)
    myform <- as.formula(paste("PD ~ ", paste(featNames2, collapse = " + "), sep = ""))
    fit <- randomForest(myform, cdat, ntree = 1000)
    serialized.model.path <- sprintf(
        "serializedModel/%s_activity_randomforest_model_excludeObjPD.rds", activity)
    saveRDS(fit, serialized.model.path)
    predict(readRDS(serialized.model.path), 
            newdata = odatC[, featNames2], type = "prob")[, "TRUE"]
}


## Shape ObjectivePD clinical data.
MapRecordHealthCodeAndRocCode <- function(cVarNames = NULL, 
                                          covNames = NULL, 
                                          version = NULL) {
    cVar <- read.delim(synGet(SYN_ID_REF$objectivePD$clinical, 
                              version = version)$path, 
                       header = TRUE, sep = "\t")
    cVar <- cVar[, c("record", "Do.you.have.Parkinson.disease.", 
                     cVarNames, covNames)]
    urec <- unique(cVar$record)
    nrec <- length(urec)
    nvar <- length(cVarNames)
    ncov <- length(covNames)
    ucVar <- data.frame(matrix(NA, nrec, 2 + nvar + ncov))
    for (i in seq(nrec)) {
        sdat <- cVar[cVar$record == urec[i],]
        ucVar[i, 1] <- urec[i]
        ucVar[i, 2] <- as.character(sdat[1, 2])
        for (j in seq(nvar)) {
            ucVar[i, 2 + j] <- median(sdat[, cVarNames[j]], na.rm = TRUE)
        }
        if (ncov > 0) {
            for (j in seq(ncov)) {
                ucVar[i, 2 + nvar + j] <- sdat[1, covNames[j]]
            } 
        }
    }
    names(ucVar) <- c("record", "PD", cVarNames, covNames)
    sampleIds <- read.delim(synGet(SYN_ID_REF$objectivePD$mapping)$path, 
                          header = TRUE, sep = "\t")
    names(sampleIds) <- c("healthCode", "rocCode")
    sampleIds <- sampleIds %>%
        dplyr::filter(rocCode != "ROC00TEST")
    rocCode2 <- as.numeric(substr(sampleIds[, 2], 5, 6))
    sampleIds <- data.frame(sampleIds, rocCode2)
    dat <- data.frame(ucVar, sampleIds[match(ucVar$record, rocCode2),])
    dat$record <- as.character(dat$record)
    dat$rocCode2 <- as.character(dat$rocCode2)
    dat$healthCode <- as.character(dat$healthCode)
    dat$PD <- as.character(dat$PD)
    return(dat[!is.na(dat$healthCode),])
}


## Merge clinical and feature data.
IncludeClinicalData <- function(fdat, cdat) {
    ids <- unique(cdat$healthCode)
    fdat <- fdat[fdat$healthCode %in% ids,]
    Dat <- data.frame(matrix(NA, nrow(fdat), ncol(cdat)))
    names(Dat) <- names(cdat)
    nids <- length(ids)
    for (i in seq(nids)) {
        idx <- which(fdat$healthCode == ids[i])
        Dat[idx,] <- cdat[i,]
    }
    return(data.frame(Dat, fdat))
}


## Merge the predicted scores.
MergeOutputs <- function(oTap, oVoi, oWal, oRes) {
    aux <- unique(c(names(oTap), names(oVoi), names(oWal), names(oRes)))
    all <- data.frame(matrix(NA, length(aux), 5))
    all[, 1] <- aux
    colnames(all) <- c("healthCode", "tapping", "voice", "walk", "rest")
    all[match(names(oTap), aux), "tapping"] <- oTap
    all[match(names(oVoi), aux), "voice"] <- oVoi
    all[match(names(oWal), aux), "walk"] <- oWal
    all[match(names(oRes), aux), "rest"] <- oRes
    return(all)
}


#########################################
#########################################
#########################################
main <- function(){
    cnms <- c(
        "UPDRS.Part.3..Motor.Exam..clinician.rated..Score.OFF.STATE",
        "Estimated.Hoehn.and.Yahr.stage.OFF.STATE",
        "UPDRS.Total.Score.OFF.STATE")
    
    metC <- MapRecordHealthCodeAndRocCode(cVarNames = cnms, covNames = NULL)
    
    ## Get ObjectivePD features.
    tapFeat <- read.delim(synGet(SYN_ID_REF$objectivePD$tapping)$path, 
                          sep = "\t", stringsAsFactors = FALSE)
    tapHCopd <- names(table(tapFeat$healthCode))
    
    voiFeat <- read.delim(synGet(SYN_ID_REF$objectivePD$voice)$path, 
                          sep = "\t", stringsAsFactors = FALSE)
    voiHCopd <- names(table(voiFeat$healthCode))
    
    walFeat <- read.delim(synGet(SYN_ID_REF$objectivePD$walking)$path, 
                          sep = "\t", stringsAsFactors = FALSE)
    walHCopd <- names(table(walFeat$healthCode))
    
    resFeat <- read.delim(synGet(SYN_ID_REF$objectivePD$resting)$path, 
                          sep = "\t", stringsAsFactors = FALSE)
    resHCopd <- names(table(resFeat$healthCode))
    
    
    tap <- IncludeClinicalData(fdat = tapFeat, cdat = metC)
    tapC <- GetCollapsedFeaturesOPD(dat = tap, 
                                    respName = "PD", 
                                    featNames = FEATURE_LIST$tap, 
                                    covNames = NULL, 
                                    negClassName = "No", 
                                    posClassName = "Yes")
    
    
    voi <- IncludeClinicalData(fdat = voiFeat, cdat = metC)
    voiC <- GetCollapsedFeaturesOPD(dat = voi, 
                                    respName = "PD", 
                                    featNames = FEATURE_LIST$voice, 
                                    covNames = NULL, 
                                    negClassName = "No", 
                                    posClassName = "Yes")
    
    
    wal <- IncludeClinicalData(fdat = walFeat, cdat = metC)
    walC <- GetCollapsedFeaturesOPD(dat = wal, 
                                    respName = "PD", 
                                    featNames = FEATURE_LIST$walk, 
                                    covNames = NULL, 
                                    negClassName = "No", 
                                    posClassName = "Yes")
    
    
    res <- IncludeClinicalData(fdat = resFeat, cdat = metC)
    resC <- GetCollapsedFeaturesOPD(dat = res, 
                                    respName = "PD", 
                                    featNames = FEATURE_LIST$rest, 
                                    covNames = NULL, 
                                    negClassName = "No", 
                                    posClassName = "Yes")
    dir.create("serializedModel")
    tapPred <- GetRfPrediction(activity = "tap", odatC = tapC, opdHC = tapHCopd)
    walPred <- GetRfPrediction(activity = "wal", odatC = walC, opdHC = walHCopd)
    resPred <- GetRfPrediction(activity = "res", odatC = resC, opdHC = resHCopd)
    voiPred <- GetRfPrediction(activity = "voi", odatC = voiC, opdHC = voiHCopd)
    
    MergeOutputs(oTap = tapPred, 
                 oVoi = voiPred, 
                 oWal = walPred, 
                 oRes = resPred) %>% 
        write.table(., OUTPUT_FILENAME, 
                    sep = "\t", quote = F,row.names = FALSE)
    
    file <- synapser::File(OUTPUT_FILENAME, 
                           parent=SYN_ID_REF$intermediate$output_folder)
    file$annotations <- ANNOTATIONS
    synStore(
        file, activity = Activity(
            'retrieve objective PD scores',
            executed = GIT_URL, 
            used = c(SYN_ID_REF$processed %>% purrr::flatten_chr(),
                     SYN_ID_REF$healthcode$case_vs_controls,
                     SYN_ID_REF$objectivePD %>% purrr::flatten_chr())))
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





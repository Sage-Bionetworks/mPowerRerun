###########################################################
#' This script is used to check consistency of PD vs Non PD assessment
#' it will not output anything to Synapse, but use this
#' to assess on whether the analysis is viable
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
############################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
source("R/utils/populationAnalysisUtils.R")
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

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
SCRIPT_NAME <- "check_consistency_PD_vs_nonPD.R"
SYN_ID_REF <- list(processed = get_processed_features_ref(),
                   intermediate = get_intermediate_data_ref())
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(
    getRepo(get("git")$repo),
    repositoryPath = file.path('R/Analyses', SCRIPT_NAME))
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder

#########################################
## Helper function
#########################################
#' For each feature, computes the 
#' differences between the median of case and 
#' the median of controls
GetCaseControlDirection <- function(cdat, featureNames) {
    idxPD <- which(cdat$PD == TRUE)
    idxNonPD <- which(cdat$PD == FALSE)
    nfeat <- length(featureNames)
    out <- matrix(NA, nfeat, 1)
    rownames(out) <- featureNames
    colnames(out) <- "PD - nonPD"
    for (i in seq(nfeat)) {
        medianPD <- median(cdat[idxPD, featureNames[i]], na.rm = TRUE)
        medianNonPD <- median(cdat[idxNonPD, featureNames[i]], na.rm = TRUE)
        out[i, 1] <- medianPD - medianNonPD
    }
    return(out)
}


#' get the medication effects for all participants (for
#' a selected feature)
GetMedicationEffects <- function(betaList, featName) {
    featcol <- which(rownames(betaList[[1]]) == featName)
    ids <- names(betaList)
    nids <- length(ids)
    betas <- matrix(NA, nids, 1)
    rownames(betas) <- ids
    colnames(betas) <- c("effect")
    for (i in seq(nids)) {
        betas[i, 1] <- as.numeric(betaList[[i]][featcol, 2])
    }
    return(betas)
}

#' get the average medication effect 
#' across all participants for all features
MedicationEffectDirection <- function(betaList, featureNames) {
    nfeat <- length(featureNames)
    out <- matrix(NA, nfeat, 1)
    rownames(out) <- featureNames
    colnames(out) <- c("effect")
    for (i in seq(nfeat)) {
        effs <- GetMedicationEffects(betaList, featureNames[i])
        out[i, ] <- as.vector(apply(effs, 2, mean))
    }
    return(out)
}


## get the signs of the effects
## (and flips the sign of the
## medication effects)
GetEffectSigns <- function(x) {
    xx <- x
    xx[, 1] <- sign(x[, 1])
    xx[, 2] <- -sign(x[, 2])
    return(xx)
}


## remove the collapsed features based on IQRs
## (it only makes sense to compare the median of the
## collapsed values in this analysis)
RemoveIQRFeatures <- function(x) {
    nms <- rownames(x)
    aux <- unlist(lapply(strsplit(nms, split = "[.]"), function (x) {x[2]}))
    keep <- which(aux == "med")
    x <- x[keep,, drop = FALSE]
    nms <- unlist(lapply(strsplit(rownames(x), split = "[.]"), function (x) {x[1]}))
    rownames(x) <- nms
    return(x)
}


## compute the consistency between the case control differences
## and the medication effects
GetEffectSignConsistency <- function(cdat,
                                     betaList,
                                     imps) {
    imp <- names(apply(imps, 1, median))
    d1 <- GetCaseControlDirection(cdat, featureNames = imp)
    d1 <- RemoveIQRFeatures(x = d1)
    nms <- unique(unlist(lapply(strsplit(imp, split = "[.]"), function (x) {x[1]})))
    d2 <- MedicationEffectDirection(betaList, featureNames = nms)
    d3 <- cbind(d1, d2[match(rownames(d1), rownames(d2)), 1])
    d4 <- GetEffectSigns(x = d3)
    tb <- table(d4[, 1], d4[, 2])
    props <- c(tb[1,1]+tb[2,2], tb[1,2]+tb[2,1])/sum(tb)
    names(props) <- c("consistent", "inconsistent")
    return(list(tb = tb, props = props, effDir = d4))
}

get.required.data <- function(){
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
        tapping = list(data = datTap, 
                       features = tapFeatures, 
                       activity = "tap"),
        walking = list(data = datWal, 
                       features = walkFeatures,
                       activity = "walk"),
        voice = list(data = datVoi, 
                     features = voiceFeatures, 
                     activity = "voice"),
        resting = list(data = datRes, 
                       features = restFeatures,
                       activity = "rest"))
    return(data.list)
}


###########################################################
## GET COLLAPSED FEATURES (SAME DATA FROM PD VS NON-PD
## CLASSIFICATION USED TO GET THE DIRECTION BETWEEN 
## PD AND NON-PD COHORT)
###########################################################

get_feature_consistency_pipeline <- function(data, features, activity, n_of_1_weights, case_vs_controls_weights){
    ## obtain numeric education and gender variables 
    ## (will be need for the causality tests later, also
    ## glmnet implementation of ridge-regression does not handle factors)
    data$education2 <- NumericEducation(x = data$education)
    data$gender2 <- BinaryGender(x = data$gender)
    ## get collapsed features (unmatched data) after filtering out
    ## participants with less than 5 records
    dat <- FilterOutParticipantsWithFewRecords(dat = data, thr = 5)
    aux <- CollapseFeatures(x = dat, 
                            labelName = "PD", 
                            covNames = c("age", "gender2", "education2"), 
                            subjectIdName = "healthCode", 
                            featNames = features)
    collapsed.data <- aux$out
    betaList <- h5read(
        n_of_1_weights, 
        glue::glue("/{activity}/LmNw", 
                   activity = activity)) %>% 
        index_user_metrics(.) %>% .$lmBetas
    imps <- h5read(
        case_vs_controls_weights, 
        glue::glue("/{activity}/unmatched/rf2/imp", 
                   activity = activity)) %>% 
        index_features(.)
    return(GetEffectSignConsistency(
        collapsed.data, betaList, imps))
}

main <- function(){
    case_vs_controls_weights <- H5Fopen(
        synGet(SYN_ID_REF$intermediate$collapsed_pd_vs_nonpd)$path)
    n_of_1_weights <- H5Fopen(
        synGet(SYN_ID_REF$intermediate$n_of_1_analysis)$path)
    data.list <- get.required.data()
    results <- plyr::llply(.data = data.list, 
                           .fun = function(activity){
                               data <- get_feature_consistency_pipeline(
                                   data = activity$data, 
                                   features = activity$features, 
                                   activity = activity$activity, 
                                   n_of_1_weights, 
                                   case_vs_controls_weights)})
    
    effDir <- rbind(results$tapping$effDir, 
                    results$walking$effDir, 
                    results$voice$effDir, 
                    results$resting$effDir)
    tb <- table(effDir[, 1], effDir[, 2])
    print(tb)
    print(round(c(tb[1,1]+tb[2,2], tb[1,2]+tb[2,1])/sum(tb), 2))
    print(c(tb[1,1]+tb[2,2], tb[1,2]+tb[2,1]))
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







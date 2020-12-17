
#############################################
#' These next 4 functions are similar
#' to the functions we used before
#' to compute the UI-pvalues, but
#' instead of returning the p-values
#' they return the position of the 
#' most significant feature, that we use 
#' for the computation of the relative
#' importances.

#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
#############################################

#########################################################
# This script is used to generate supplementary figure 8-9,
# and figure for main text 2 of the mPower manuscript
##########################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")
source("R/utils/personalizedAnalysisUtils.R")

######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(
    readLines(get("git")$path))

######################################################
## Variables and Synapse References
#######################################################
SYN_ID_REF <- list(figures = get_figure_ref(),
                   table = get_synapse_table_ref(),
                   intermediate = get_intermediate_data_ref())
SCRIPT_NAME <-  "Nof1_figures.R"
GIT_URL <- getPermlink(
    getRepo(get("git")$repo,
            ref="branch", 
            refName=get("git")$branch), 
    repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
N_OF_1_ANALYSIS_SYN_ID <- SYN_ID_REF$intermediate$n_of_1_analysis

ANNOTATIONS <- list(
    analysisType = "n of 1 analysis",
    analysisSubtype = "treatment vs tod - relative importance",
    userSubset = tolower(get("metadata")$user_group), 
    pipelineStep = "figures")

MAIN_FIGURE_2 <- paste0("mPower_",
    gsub(" ", "_", get("metadata")$user_group), 
    "_main_text_figure_2",".png")

SUPPL_FIGURE_8 <- paste0("mPower_",
    gsub(" ", "_", get("metadata")$user_group), 
    "_supplementary_figure_8",".png")

SUPPL_FIGURE_9 <- paste0("mPower_",
    gsub(" ", "_", get("metadata")$user_group), 
    "_supplementary_figure_9",".png")

#' create logger for pipeline
sink('pipeline.log', append = TRUE)
cat(paste0(
    "[",Sys.time(), "]", " Running ", SCRIPT_NAME), "\n\n")
sink()

## correction across individuals
myMtMethod <- "BH" 
## correction across the 5 conditional independence tests per feature
mtMethod1 = "BH" 
## correction across the features
mtMethod2 = "BH"
## significance threshold
thr <- 0.05

UItestsTopPutativeTreat <- function(xL, mtMethod1 = "BH", mtMethod2 = "BH", thr = 0.05) {
    ids <- names(xL)
    nids <- length(ids)
    nfeat <- nrow(xL[[1]])
    topPvals <- matrix(NA, nids, 1)
    rownames(topPvals) <- ids
    colnames(topPvals) <- "pvalUItest"
    for (i in seq(nids)) {
        topPvals[i, 1] <- ParticipantTopPutativeTreatPval(xL[[i]], mtMethod1, mtMethod2, thr)
    }
    
    topPvals
}


UItestsTopPutativeTod <- function(xL, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
    ids <- names(xL)
    nids <- length(ids)
    nfeat <- nrow(xL[[1]])
    topPvals <- matrix(NA, nids, 1)
    rownames(topPvals) <- ids
    colnames(topPvals) <- "pvalUItest"
    for (i in seq(nids)) {
        topPvals[i, 1] <- ParticipantTopPutativeTodPval(xL[[i]], mtMethod1, mtMethod2, thr)
    }
    
    topPvals
}


ParticipantTopPutativeTreatPval <- function(x, mtMethod1 = "BH", mtMethod2 = "BH", thr = 0.05) {
    ev <- ParticipantEvidence(x, mtMethod1, thr)$ev
    ## We only adjust for tod if we detect a consistent 
    ## tod (or tod and treat) effect. We do not adjust
    ## for tod if: 
    ## evidence == "treat effect"
    ## evidence == "inconsistent treat effect"
    ## evidence == "inconsistent treat and tod effect"
    ## evidence == "inconsistent tod effect"
    ## In these cases we replace evidence by NA,
    ## so that we can get the indexes of the cases where
    ## we will not adjust for tod (idxMarg) by asking which cases 
    ## are NA, and get the indexes of the cases where we will 
    ## adjust for tod (idxcond) by asking which cases are
    ## not NA.
    ev[ev == "treat effect"] <- NA 
    ev[ev == "inconsistent treat effect"] <- NA 
    ev[ev == "inconsistent treat and tod effect"] <- NA 
    ev[ev == "inconsistent tod effect"] <- NA
    idxMarg <- which(is.na(ev))
    idxCond <- which(!is.na(ev))
    nfeat <- length(ev)
    pvals <- rep(NA, nfeat)
    names(pvals) <- names(ev)
    pvals[idxMarg] <- x[idxMarg, "a(Y,X)=0"]
    pvals[idxCond] <- x[idxCond, "a(Y,X|T)=0"]
    ## multiple testing correction across all features
    pvals <- p.adjust(pvals, method = mtMethod2)
    
    which.min(pvals) ## position of the union-intersection p-value
}


ParticipantTopPutativeTodPval <- function(x, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
    ev <- ParticipantEvidence(x, mtMethod1, thr)$ev
    ## We only adjust for treat if we detect a consistent 
    ## treat (or tod and treat) effect.
    ev[ev == "tod effect"] <- NA 
    ev[ev == "inconsistent tod effect"] <- NA 
    ev[ev == "inconsistent treat and tod effect"] <- NA 
    ev[ev == "inconsistent treat effect"] <- NA 
    nfeat <- length(ev)
    pvals <- rep(NA, nfeat)
    names(pvals) <- names(ev)
    idxMarg <- which(is.na(ev))
    idxCond <- which(!is.na(ev))
    pvals[idxMarg] <- x[idxMarg, "a(Y,T)=0"]
    pvals[idxCond] <- x[idxCond, "a(Y,T|X)=0"]
    pvals <- p.adjust(pvals, method = mtMethod2)
    
    which.min(pvals)
}


## grab the relative importance of the most significant feature
## from the output of the relative importance code
GetRelativeImportanceOfMostSignificantFeature <- function(ri, topFeatTre, topFeatTod) {
    ids <- names(ri)
    nids <- length(ids)
    out <- matrix(NA, nids, 2)
    rownames(out) <- ids
    colnames(out) <- c("tre", "tod")
    for (i in seq(nids)) {
        out[i, 1] <- ri[[i]][topFeatTre[i], "medTimepoint"]
        out[i, 2] <- ri[[i]][topFeatTod[i], "tod"]
    }
    
    out
}


## separates the relative importance values between individuals
## responders and non-responders
GetRespondersVsNonRespondersRelativeImportance <- function(ui, ri, fdr.thr = 0.05) {
    thr <- -log(fdr.thr, 10)
    responderRI <- c()
    nonResponderRI <- c()
    for (i in seq(4)) {
        idxR <- which(ui[, i] >= thr)
        idxNR <- which(ui[, i] < thr)
        responderRI <- c(responderRI, ri[idxR, i])
        nonResponderRI <- c(nonResponderRI, ri[idxNR, i])
    }
    
    list(responderRI = responderRI,
         nonResponderRI = nonResponderRI)
}


## make colors transparent (grabed it from the web)
makeTransparent = function(..., alpha=0.5) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
}

######################################
## shape relative importance scores
######################################
n_of_1_analysis <- H5Fopen(synGet(N_OF_1_ANALYSIS_SYN_ID)$path)

outTapLmNw <- index_user_metrics(h5read(n_of_1_analysis, "/tap/LmNw"))
outTapArima <- index_user_metrics(h5read(n_of_1_analysis, "/tap/arima"))

outVoiLmNw <- index_user_metrics(h5read(n_of_1_analysis, "/voice/LmNw"))
outVoiArima <- index_user_metrics(h5read(n_of_1_analysis, "/voice/arima"))

outResLmNw <- index_user_metrics(h5read(n_of_1_analysis, "/rest/LmNw"))
outResArima <- index_user_metrics(h5read(n_of_1_analysis, "/rest/arima"))

outWalLmNw <- index_user_metrics(h5read(n_of_1_analysis, "/walk/LmNw"))
outWalArima <- index_user_metrics(h5read(n_of_1_analysis, "/walk/arima"))

riTap <- h5read(n_of_1_analysis, "/tap/ri")
riVoi <- h5read(n_of_1_analysis, "/voice/ri")
riWal <- h5read(n_of_1_analysis, "/walk/ri")
riRes <- h5read(n_of_1_analysis, "/rest/ri")

######################################
## shape relative importance scores
######################################

## grab relative importances for the most
## significant tapping feature
outLmNw <- outTapLmNw
outArima <- outTapArima
uiTreLm <- UItestsTopPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsTopPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsTopPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsTopPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsTopPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsTopPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
topFeatTreTapNW <- uiTre$uiTreNw
topFeatTodTapNW <- uiTod$uiTodNw
topRiTapNW <- GetRelativeImportanceOfMostSignificantFeature(
    riTap, topFeatTreTapNW, topFeatTodTapNW)
topFeatTreTapLM <- uiTre$uiTreLm
topFeatTodTapLM <- uiTod$uiTodLm
topRiTapLM <- GetRelativeImportanceOfMostSignificantFeature(riTap, topFeatTreTapLM, topFeatTodTapLM)
topFeatTreTapAR <- uiTre$uiTreAr
topFeatTodTapAR <- uiTod$uiTodAr
topRiTapAR <- GetRelativeImportanceOfMostSignificantFeature(riTap, topFeatTreTapAR, topFeatTodTapAR)


## grab relative importances for the most
## significant voice feature
outLmNw <- outVoiLmNw
outArima <- outVoiArima
uiTreLm <- UItestsTopPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsTopPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsTopPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsTopPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsTopPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsTopPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
topFeatTreVoiNW <- uiTre$uiTreNw
topFeatTodVoiNW <- uiTod$uiTodNw
topRiVoiNW <- GetRelativeImportanceOfMostSignificantFeature(riVoi, topFeatTreVoiNW, topFeatTodVoiNW)
topFeatTreVoiLM <- uiTre$uiTreLm
topFeatTodVoiLM <- uiTod$uiTodLm
topRiVoiLM <- GetRelativeImportanceOfMostSignificantFeature(riVoi, topFeatTreVoiLM, topFeatTodVoiLM)
topFeatTreVoiAR <- uiTre$uiTreAr
topFeatTodVoiAR <- uiTod$uiTodAr
topRiVoiAR <- GetRelativeImportanceOfMostSignificantFeature(riVoi, topFeatTreVoiAR, topFeatTodVoiAR)


## grab relative importances for the most
## significant walk feature
outLmNw <- outWalLmNw
outArima <- outWalArima
uiTreLm <- UItestsTopPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsTopPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsTopPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsTopPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsTopPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsTopPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
topFeatTreWalNW <- uiTre$uiTreNw
topFeatTodWalNW <- uiTod$uiTodNw
topRiWalNW <- GetRelativeImportanceOfMostSignificantFeature(riWal, topFeatTreWalNW, topFeatTodWalNW)
topFeatTreWalLM <- uiTre$uiTreLm
topFeatTodWalLM <- uiTod$uiTodLm
topRiWalLM <- GetRelativeImportanceOfMostSignificantFeature(riWal, topFeatTreWalLM, topFeatTodWalLM)
topFeatTreWalAR <- uiTre$uiTreAr
topFeatTodWalAR <- uiTod$uiTodAr
topRiWalAR <- GetRelativeImportanceOfMostSignificantFeature(riWal, topFeatTreWalAR, topFeatTodWalAR)


## grab relative importances for the most
## significant rest feature
outLmNw <- outResLmNw
outArima <- outResArima
uiTreLm <- UItestsTopPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsTopPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsTopPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsTopPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsTopPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsTopPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
topFeatTreResNW <- uiTre$uiTreNw
topFeatTodResNW <- uiTod$uiTodNw
topRiResNW <- GetRelativeImportanceOfMostSignificantFeature(riRes, topFeatTreResNW, topFeatTodResNW)
topFeatTreResLM <- uiTre$uiTreLm
topFeatTodResLM <- uiTod$uiTodLm
topRiResLM <- GetRelativeImportanceOfMostSignificantFeature(riRes, topFeatTreResLM, topFeatTodResLM)
topFeatTreResAR <- uiTre$uiTreAr
topFeatTodResAR <- uiTod$uiTodAr
topRiResAR <- GetRelativeImportanceOfMostSignificantFeature(riRes, topFeatTreResAR, topFeatTodResAR)


## Newey-West results
riTreNw <- MergeOutputs(oTap = topRiTapNW[, "tre", drop = F], 
                        oVoi = topRiVoiNW[, "tre", drop = F], 
                        oWal = topRiWalNW[, "tre", drop = F], 
                        oRes = topRiResNW[, "tre", drop = F], 
                        colName = "tre")
riTodNw <- MergeOutputs(oTap = topRiTapNW[, "tod", drop = F], 
                        oVoi = topRiVoiNW[, "tod", drop = F], 
                        oWal = topRiWalNW[, "tod", drop = F], 
                        oRes = topRiResNW[, "tod", drop = F], 
                        colName = "tod")


## naive standard linear model results
riTreLm <- MergeOutputs(oTap = topRiTapLM[, "tre", drop = F], 
                        oVoi = topRiVoiLM[, "tre", drop = F], 
                        oWal = topRiWalLM[, "tre", drop = F], 
                        oRes = topRiResLM[, "tre", drop = F], 
                        colName = "tre")
riTodLm <- MergeOutputs(oTap = topRiTapLM[, "tod", drop = F], 
                        oVoi = topRiVoiLM[, "tod", drop = F], 
                        oWal = topRiWalLM[, "tod", drop = F], 
                        oRes = topRiResLM[, "tod", drop = F], 
                        colName = "tod")


## ARIMA results
riTreAr <- MergeOutputs(oTap = topRiTapAR[, "tre", drop = F], 
                        oVoi = topRiVoiAR[, "tre", drop = F], 
                        oWal = topRiWalAR[, "tre", drop = F], 
                        oRes = topRiResAR[, "tre", drop = F], 
                        colName = "tre")
riTodAr <- MergeOutputs(oTap = topRiTapAR[, "tod", drop = F], 
                        oVoi = topRiVoiAR[, "tod", drop = F], 
                        oWal = topRiWalAR[, "tod", drop = F], 
                        oRes = topRiResAR[, "tod", drop = F], 
                        colName = "tod")





######################################
## shape union-intersection p-values
######################################


## get corrected union-intersection p-values from the tapping outputs 
outLmNw <- outTapLmNw
outArima <- outTapArima
uiTreLm <- UItestsPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
uiTreTap <- apply(uiTre, 2, p.adjust, myMtMethod)
uiTodTap <- apply(uiTod, 2, p.adjust, myMtMethod)


## get corrected union-intersection p-values from the voice outputs 
outLmNw <- outVoiLmNw
outArima <- outVoiArima
uiTreLm <- UItestsPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
uiTreVoi <- apply(uiTre, 2, p.adjust, myMtMethod)
uiTodVoi <- apply(uiTod, 2, p.adjust, myMtMethod)


## get corrected union-intersection p-values from the walk outputs 
outLmNw <- outWalLmNw
outArima <- outWalArima
uiTreLm <- UItestsPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
uiTreWal <- apply(uiTre, 2, p.adjust, myMtMethod)
uiTodWal <- apply(uiTod, 2, p.adjust, myMtMethod)


## get corrected union-intersection p-values from the balance outputs 
outLmNw <- outResLmNw
outArima <- outResArima
uiTreLm <- UItestsPutativeTreat(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTreNw <- UItestsPutativeTreat(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTreAr <- UItestsPutativeTreat(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTre <- data.frame(uiTreLm, uiTreNw, uiTreAr)
names(uiTre) <- c("uiTreLm", "uiTreNw", "uiTreAr")
uiTodLm <- UItestsPutativeTod(x = outLmNw$lmPvals, mtMethod1, mtMethod2, thr = thr)
uiTodNw <- UItestsPutativeTod(x = outLmNw$nwPvals, mtMethod1, mtMethod2, thr = thr)
uiTodAr <- UItestsPutativeTod(x = outArima$arimaPvals, mtMethod1, mtMethod2, thr = thr)
uiTod <- data.frame(uiTodLm, uiTodNw, uiTodAr)
names(uiTod) <- c("uiTodLm", "uiTodNw", "uiTodAr")
uiTreRes <- apply(uiTre, 2, p.adjust, myMtMethod)
uiTodRes <- apply(uiTod, 2, p.adjust, myMtMethod)


## get proportions of detected effects
fTap <- GetProportions(uiTreTap, uiTodTap, thr = 0.05, r = 2, task = "tapping")
fVoi <- GetProportions(uiTreVoi, uiTodVoi, thr = 0.05, r = 2, task = "voice")
fWal <- GetProportions(uiTreWal, uiTodWal, thr = 0.05, r = 2, task = "walk")
fRes <- GetProportions(uiTreRes, uiTodRes, thr = 0.05, r = 2, task = "rest")
ff <- cbind(fTap, fVoi, fWal, fRes)
ff[1:2,]


#################################################
## get demographics data and merge it with the
## union intersection p-values
#################################################

## get demographics data (contains the disease onset year, which will be
## correlated with the union-intersection p-values in Figure 2 and
## Supplenetary Figure 8).
demo <- synTableQuery("SELECT * FROM syn7222419")$asDataFrame()
names(demo) <- make.names(names(demo))

## the union-intersection p-values are already adjusted for multiple testing,
## so no further adjustment necessary
mtmethod <- "none" 

## merge union-intersection p-values and the demographics data
datUiTreTap <- CollateUipvalsAndDemo(uiTreTap, demo, mtmethod = mtmethod)
datUiTodTap <- CollateUipvalsAndDemo(uiTodTap, demo, mtmethod = mtmethod)
dim(datUiTreTap)

datUiTreVoi <- CollateUipvalsAndDemo(uiTreVoi, demo, mtmethod = mtmethod)
datUiTodVoi <- CollateUipvalsAndDemo(uiTodVoi, demo, mtmethod = mtmethod)
dim(datUiTreVoi)

datUiTreWal <- CollateUipvalsAndDemo(uiTreWal, demo, mtmethod = mtmethod)
datUiTodWal <- CollateUipvalsAndDemo(uiTodWal, demo, mtmethod = mtmethod)
dim(datUiTreWal)

datUiTreRes <- CollateUipvalsAndDemo(uiTreRes, demo, mtmethod = mtmethod)
datUiTodRes <- CollateUipvalsAndDemo(uiTodRes, demo, mtmethod = mtmethod)
dim(datUiTreRes)

datUiTre <- rbind(datUiTreTap, datUiTreVoi, datUiTreWal, datUiTreRes)
datUiTod <- rbind(datUiTodTap, datUiTodVoi, datUiTodWal, datUiTodRes)


## merge outputs of the standard linear regression analyses
uiTreLm <- MergeOutputs(oTap = uiTreTap, oVoi = uiTreVoi, oWal = uiTreWal, oRes = uiTreRes, colName = "uiTreLm")
uiTodLm <- MergeOutputs(oTap = uiTodTap, oVoi = uiTodVoi, oWal = uiTodWal, oRes = uiTodRes, colName = "uiTodLm")
uiTreLm <- -log(uiTreLm, 10)
uiTodLm <- -log(uiTodLm, 10)


## merge outputs of the regression with ARIMA error analyses
uiTreArima <- MergeOutputs(oTap = uiTreTap, oVoi = uiTreVoi, oWal = uiTreWal, oRes = uiTreRes, colName = "uiTreAr")
uiTodArima <- MergeOutputs(oTap = uiTodTap, oVoi = uiTodVoi, oWal = uiTodWal, oRes = uiTodRes, colName = "uiTodAr")
uiTreArima <- -log(uiTreArima, 10)
uiTodArima <- -log(uiTodArima, 10)


## merge outputs of the Newey-West regression analyses
uiTreNw <- MergeOutputs(oTap = uiTreTap, oVoi = uiTreVoi, oWal = uiTreWal, oRes = uiTreRes, colName = "uiTreNw")
uiTodNw <- MergeOutputs(oTap = uiTodTap, oVoi = uiTodVoi, oWal = uiTodWal, oRes = uiTodRes, colName = "uiTodNw")
uiTreNw <- -log(uiTreNw, 10)
uiTodNw <- -log(uiTodNw, 10)

## order the results according to the linear regression results
o <- order(apply(uiTreLm, 1, mean, na.rm = TRUE), decreasing = TRUE)


## ploting parameters
mycex <- 1
rcolor <- "darkgreen"
nrcolor <- "#440154FF"
tapCol <- "#7E6148"
voiCol <- "#E64B36"
walCol <- "#8491B4"
resCol <- "#00A087"
alpha <- 0.05
cm <- 1.5
cp <- 0.9
cexleg2 <- 0.75




######################################
## Updated Figure 2 (main text)
######################################

trcolor <- makeTransparent(rcolor, alpha = 0.5)
tnrcolor <- makeTransparent(nrcolor, alpha = 0.5)
nc <- 15
myxlim <- c(0, 0.65)
myylim <- c(0, 14)

lmat <- matrix(c(1, 1, 3, 4,
                 2, 2, 5, 6), 2, 4, byrow = TRUE)

figpath <- ""

title <- MAIN_FIGURE_2
png(title, width = 2500, height = 1300, res = 200)
layout(lmat)
par(mar = c(3, 3.4, 1.75, 1) + 0.1, mgp = c(2, 0.75, 0))
plot(uiTreNw[o, "tapping"], col = "blue", pch = 19, xaxt = "n", xlab = "", 
     ylab = "-log10(FDR)", type = "n",
     main = "Medication", ylim = c(0, 35), cex.main = cm)
points(uiTreNw[o, "tapping"], col = tapCol, pch = 19, cex = mycex)
points(uiTreNw[o, "voice"], col = voiCol, pch = 19, cex = mycex)
points(uiTreNw[o, "walk"], col = walCol, pch = 19, cex = mycex)
points(uiTreNw[o, "rest"], col = resCol, pch = 19, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("Assays:", "tapping", "voice", "walk", "balance"), 
       text.col = c("black", tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
axis(side = 1, labels = FALSE, at = seq(nrow(uiTreNw)))
mtext(side = 1, "participants", las = 1, at = nrow(uiTreNw)/2, cex = cp, line = 1.5)
mtext(side = 2, "(a)", las = 1, at = 35, cex = 1, line = -2)
################
plot(uiTodNw[o, "tapping"], col = "blue", pch = 17, xaxt = "n", xlab = "", 
     ylab = "-log10(FDR)", type = "n",
     main = "Time-of-the-day", ylim = c(0, 35))
points(uiTodNw[o, "tapping"], col = tapCol, pch = 17, cex = mycex)
points(uiTodNw[o, "voice"], col = voiCol, pch = 17, cex = mycex)
points(uiTodNw[o, "walk"], col = walCol, pch = 17, cex = mycex)
points(uiTodNw[o, "rest"], col = resCol, pch = 17, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("Assays:", "tapping", "voice", "walk", "balance"), 
       text.col = c("black", tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
axis(side = 1, labels = FALSE, at = seq(nrow(uiTodNw)))
mtext(side = 1, "participants", las = 1, at = nrow(uiTodNw)/2, cex = cp, line = 1.5)
mtext(side = 2, "(d)", las = 1, at = 35, cex = 1, line = -2)
###
###
aux <- RespondersAndNonRespondersFits(dat = datUiTre, 
                                      respName = "logPvalUiNw", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiNw", 
                                      logPvalThr = -log(alpha, 10))
plot(aux$form, data = aux$dat, type = "n", main = "Medication", 
     ylab = "-log10(FDR)",
     xlab = "disease onset year", ylim = c(0, 35))
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 19, cex = cp)
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 19, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
mylegend <- c("", "Fluctuators", "Non-fluctuators", 
              "Estimate", round(aux$bR, 4), round(aux$bNR, 4), 
              "p-value", round(aux$pvalR, 4), round(aux$pvalNR, 4))
legend("top", ncol = 3, 
       legend = mylegend, 
       text.col = c("black", rcolor, nrcolor, 
                    "black", rcolor, nrcolor, 
                    "black", rcolor, nrcolor),
       cex = cexleg2)
mtext("(b)", side = 2, at = 25, cex = 1, las = 1, line = -2)
#################
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTreNw, ri = riTreNw, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim, 
     xlab = "relative importance", probability = TRUE, main = "Medication")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("Fluctuators", "Non-fluctuators"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.2)
mtext("(c)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
#################
aux <- RespondersAndNonRespondersFits(dat = datUiTod, 
                                      respName = "logPvalUiNw", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiNw", 
                                      logPvalThr = -log(alpha, 10))
plot(aux$form, data = aux$dat, type = "n", main = "Time-of-the-day", 
     ylab = "-log10(FDR)",
     xlab = "disease onset year", ylim = c(0, 35))
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 17, cex = cp)
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 17, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
mylegend <- c("", "Fluctuators", "Non-fluctuators", 
              "Estimate", round(aux$bR, 4), round(aux$bNR, 4), 
              "p-value", round(aux$pvalR, 4), round(aux$pvalNR, 4))
legend("top", ncol = 3, 
       legend = mylegend, 
       text.col = c("black", rcolor, nrcolor, 
                    "black", rcolor, nrcolor, 
                    "black", rcolor, nrcolor),
       cex = cexleg2)
mtext("(e)", side = 2, at = 25, cex = 1, las = 1, line = -2)
#################
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTodNw, ri = riTodNw, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim,
     xlab = "relative importance", probability = TRUE, main = "Time-of-the-day")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("Fluctuators", "Non-fluctuators"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.2)
mtext("(f)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
dev.off()

#############################
## Store to Synapse
#############################
f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
f$annotations <- ANNOTATIONS
synStore(f, activity = Activity(executed = GIT_URL, used = c(N_OF_1_ANALYSIS_SYN_ID)))
unlink(title)


######################################
## Updated Supplementary Figure 8
######################################
lmat <- matrix(c(1, 1, 5, 9,
                 2, 2, 6, 10,
                 3, 3, 7, 11,
                 4, 4, 8, 12), 4, 4, byrow = TRUE)

figpath <- ""

myylim <- c(0, 25)

title <- SUPPL_FIGURE_8
png(title, width = 2500, height = 2000, res = 200)
layout(lmat)
par(mar = c(0, 4, 1.75, 0.75) + 0.1, mgp = c(2, 0.75, 0))
plot(uiTreLm[o, "tapping"], col = "blue", pch = 19, xaxt = "n", xlab = "", 
     ylab = "-log10(treatment FDR)", type = "n",
     main = "Linear regression", ylim = c(0, 55), cex.main = cm)
points(uiTreLm[o, "tapping"], col = tapCol, pch = 19, cex = mycex)
points(uiTreLm[o, "voice"], col = voiCol, pch = 19, cex = mycex)
points(uiTreLm[o, "walk"], col = walCol, pch = 19, cex = mycex)
points(uiTreLm[o, "rest"], col = resCol, pch = 19, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("tapping", "voice", "walk", "balance"), 
       text.col = c(tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
text(54, 50, "Putative treatment effect", cex = 1.5)
mtext(side = 2, "(a)", las = 1, at = 52, cex = 1, line = -2)
################
par(mar = c(3, 4, 0, 0.75) + 0.1)
plot(uiTodLm[o, "tapping"], col = "blue", pch = 17, xaxt = "n", xlab = "", 
     ylab = "-log10(time-of-the-day FDR)", type = "n",
     main = "", ylim = c(0, 55))
points(uiTodLm[o, "tapping"], col = tapCol, pch = 17, cex = mycex)
points(uiTodLm[o, "voice"], col = voiCol, pch = 17, cex = mycex)
points(uiTodLm[o, "walk"], col = walCol, pch = 17, cex = mycex)
points(uiTodLm[o, "rest"], col = resCol, pch = 17, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("tapping", "voice", "walk", "balance"), 
       text.col = c(tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
mtext(side = 1, "participants", las = 1, at = nrow(uiTodLm)/2, cex = cp, line = 0.8)
text(54, 50, "Putative time-of-the-day effect", cex = 1.5)
mtext(side = 2, "(d)", las = 1, at = 52, cex = 1, line = -2)
###
###
par(mar = c(0, 4, 1.75, 0.75) + 0.1, mgp = c(2, 0.75, 0))
plot(uiTreArima[o, "tapping"], col = "blue", pch = 19, xaxt = "n", xlab = "", 
     ylab = "-log10(treatment FDR)", type = "n",
     main = "ARIMA", ylim = c(0, 55), cex.main = cm)
points(uiTreArima[o, "tapping"], col = tapCol, pch = 19, cex = mycex)
points(uiTreArima[o, "voice"], col = voiCol, pch = 19, cex = mycex)
points(uiTreArima[o, "walk"], col = walCol, pch = 19, cex = mycex)
points(uiTreArima[o, "rest"], col = resCol, pch = 19, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("tapping", "voice", "walk", "balance"), 
       text.col = c(tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
text(54, 50, "Putative treatment effect", cex = 1.5)
mtext(side = 2, "(g)", las = 1, at = 52, cex = 1, line = -2)
################
par(mar = c(3, 4, 0, 0.75) + 0.1)
plot(uiTodArima[o, "tapping"], col = "blue", pch = 17, xaxt = "n", xlab = "", 
     ylab = "-log10(time-of-the-day FDR)", type = "n",
     main = "", ylim = c(0, 55))
points(uiTodArima[o, "tapping"], col = tapCol, pch = 17, cex = mycex)
points(uiTodArima[o, "voice"], col = voiCol, pch = 17, cex = mycex)
points(uiTodArima[o, "walk"], col = walCol, pch = 17, cex = mycex)
points(uiTodArima[o, "rest"], col = resCol, pch = 17, cex = mycex)
abline(h = -log(alpha, 10), col = "black", lwd = 1)
legend("topright", legend = c("tapping", "voice", "walk", "balance"), 
       text.col = c(tapCol, voiCol, walCol, resCol), bty = "n", cex = 1.2)
mtext(side = 1, "participants", las = 1, at = nrow(uiTodArima)/2, cex = cp, line = 0.8)
text(54, 50, "Putative time-of-the-day effect", cex = 1.5)
mtext(side = 2, "(j)", las = 1, at = 52, cex = 1, line = -2)
###
###
par(mar = c(0, 4, 1.75, 0.75) + 0.1, mgp = c(2, 0.75, 0))
aux <- RespondersAndNonRespondersFits(dat = datUiTre, 
                                      respName = "logPvalUiLm", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiLm", 
                                      logPvalThr = -log(alpha, 10))
auxTreLM <- aux
plot(aux$form, data = aux$dat, type = "n", main = "Linear-regression", 
     ylab = "-log10(treatment FDR)",
     xlab = "disease onset year", ylim = c(0, 55), xaxt = "n")
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 19, cex = cp)
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 19, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
legend("topright", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n")
mtext("(b)", side = 2, at = 52, cex = 1, las = 1, line = -2)
####
par(mar = c(3, 4, 0, 0.75) + 0.1)
aux <- RespondersAndNonRespondersFits(dat = datUiTod, 
                                      respName = "logPvalUiLm", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiLm", 
                                      logPvalThr = -log(alpha, 10))
auxTodLM <- aux
plot(aux$form, data = aux$dat, type = "n", main = "", 
     ylab = "-log10(time-of-the-day FDR)",
     xlab = "disease onset year", ylim = c(0, 55))
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 17, cex = cp)
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 17, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
legend("topright", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n")
mtext("(e)", side = 2, at = 52, cex = 1, las = 1, line = -2)
###
###
par(mar = c(0, 4, 1.75, 0.75) + 0.1, mgp = c(2, 0.75, 0))
aux <- RespondersAndNonRespondersFits(dat = datUiTre, 
                                      respName = "logPvalUiAr", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiAr", 
                                      logPvalThr = -log(alpha, 10))
auxTreAR <- aux
plot(aux$form, data = aux$dat, type = "n", main = "ARIMA", 
     ylab = "-log10(treatment FDR)",
     xlab = "disease onset year", ylim = c(0, 55), xaxt = "n")
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 19, cex = cp)
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 19, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
legend("topright", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n")
mtext("(h)", side = 2, at = 52, cex = 1, las = 1, line = -2)
####
par(mar = c(3, 4, 0, 0.75) + 0.1)
aux <- RespondersAndNonRespondersFits(dat = datUiTod, 
                                      respName = "logPvalUiAr", 
                                      covName = "onset.year", 
                                      logPvalName = "logPvalUiAr", 
                                      logPvalThr = -log(alpha, 10))
auxTodAR <- aux
plot(aux$form, data = aux$dat, type = "n", main = "", 
     ylab = "-log10(time-of-the-day FDR)",
     xlab = "disease onset year", ylim = c(0, 55))
points(aux$dat[aux$dat$responder==TRUE, aux$covName], 
       aux$dat[aux$dat$responder==TRUE, aux$respName], 
       col = rcolor, pch = 17, cex = cp)
points(aux$dat[aux$dat$responder==FALSE, aux$covName], 
       aux$dat[aux$dat$responder==FALSE, aux$respName], 
       col = nrcolor, pch = 17, cex = cp)
abline(a = aux$aR, b = aux$bR, lwd = 2, col = rcolor)
abline(a = aux$aNR, b = aux$bNR, lwd = 2, col = nrcolor)
legend("topright", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n")
mtext("(k)", side = 2, at = 52, cex = 1, las = 1, line = -2)
###
par(mar = c(4, 4, 1.5, 0.75) + 0.1)
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTreLm, ri = riTreLm, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim, 
     xlab = "relative importance", probability = TRUE, main = "Lin.-regr. - medication")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.0)
mtext("(c)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
###
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTodLm, ri = riTodLm, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim, 
     xlab = "relative importance", probability = TRUE, main = "Lin.-regr. - ToD")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.0)
mtext("(f)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
###
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTreArima, ri = riTreAr, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim, 
     xlab = "relative importance", probability = TRUE, main = "ARIMA - medication")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.0)
mtext("(i)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
###
out <- GetRespondersVsNonRespondersRelativeImportance(ui = uiTodArima, ri = riTodAr, fdr.thr = 0.05)
hist(out$nonResponderRI, col = tnrcolor, xlim = myxlim, nclass = nc, ylim = myylim, 
     xlab = "relative importance", probability = TRUE, main = "ARIMA - ToD")
hist(out$responderRI, col = trcolor, nclass = nc, add = TRUE, probability = TRUE)
legend("right", legend = c("responders", "non-responders"), 
       text.col = c(rcolor, nrcolor), bty = "n", cex = 1.0)
mtext("(l)", side = 3, at = 0.6, cex = 1, las = 1, line = -2)
dev.off()

#############################
## Store to Synapse
#############################
f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
f$annotations <- ANNOTATIONS
synStore(f, activity = Activity(executed = GIT_URL, used = c(N_OF_1_ANALYSIS_SYN_ID)))
unlink(title)




######################################
## Updated Supplementary Figure 9
######################################

aux <- MergeOutputs(oTap = uiTreTap, oVoi = uiTreVoi, oWal = uiTreWal, oRes = uiTreRes, colName = "uiTreNw")
aux[aux <= 0.05] <- -1
aux[aux > 0.05] <- 1
aux[is.na(aux)] <- 0
colnames(aux) <- c("tapping", "voice", "walk", "balance")

lmat <- rbind( c(4,3,4), c(2,1,4) )
lhei <- c(1.0, 4)
lwid <- c(0.5, 4, 0.5)


title <- SUPPL_FIGURE_9
png(title, width = 1800, height = 3000, res = 300)
heatmap.2(aux, col = colorpanel(n = 124, low = "darkred", mid = "lightgrey", high = "darkblue"), 
          trace = "none", key = FALSE, margins = c(8, 3), Colv = TRUE, Rowv = TRUE,
          dendrogram = "none", main = "Treatment responders", ylab = "participants",
          xlab = "assessments", labRow = FALSE, colsep = seq(0, 4), rowsep = seq(0, 131),
          sepcolor = "black", sepwidth = 0.001,
          lmat = lmat, lhei = lhei, lwid = lwid)
legend("topright", legend = c("Responders", "Non-responders", "Not available"),
       fill = c("darkred", "darkblue", "lightgrey"), bty = "n")
dev.off()

#############################
## Store to Synapse
#############################
f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
f$annotations <- ANNOTATIONS
synStore(f, activity = Activity(executed = GIT_URL, used = c(N_OF_1_ANALYSIS_SYN_ID)))
unlink(title)
sink('pipeline.log', append = TRUE)
cat(paste0("[",Sys.time(), "]", " Done Running N of 1 Analysis"), "\n\n")
sink()







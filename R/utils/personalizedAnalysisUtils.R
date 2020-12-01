##############################################
## Utility Functions used for N of 1 Analysis
############################################

# Dependencies
library(randomForest)
library(pROC)
library(forecast)
library(gplots)
library(sandwich)
library(lmtest)
library(relaimpo)


####################################
## Helper Function
####################################

#' Retrieve Data for N of 1
#'
#' Function to retrieve data for N of 1 Analysis
#'
#' @param dat featurized dataset (dataframe, tibble)
#' @param beforeThr threshold of records before medication
#' @param beforeThr threshold of records after medication
#' @examples
#' data = GetDataForNof1(tap_data, 15, 15)
GetDataForNof1 <- function(dat, beforeThr, afterThr) {
  CountBeforeAfterMedicationPerParticipant <- function(dat) {
    participantIds <- unique(dat$healthCode)
    nParticipants <- length(participantIds)
    out <- data.frame(matrix(NA, nParticipants, 4))
    colnames(out) <- c("healthCode", "n", "nBefore", "nAfter")
    out[, 1] <- participantIds
    for (i in seq(nParticipants)) {
      #cat(i, "\n")
      pdat <- dat[which(dat$healthCode == participantIds[i]),]
      aux <- as.character(pdat$medTimepoint)
      out[i, "n"] <- nrow(pdat)
      out[i, "nBefore"] <- sum(aux == "Immediately before Parkinson medication", na.rm = TRUE)
      out[i, "nAfter"] <- sum(aux == "Just after Parkinson medication (at your best)", na.rm = TRUE)
    }
    out[order(out[, "nBefore"] + out[, "nAfter"], decreasing = TRUE),]
  }
  GetParticipants <- function(x, beforeThr, afterThr) {
    idx <- which(x$nBefore >= beforeThr & x$nAfter >= afterThr)
    x$healthCode[idx]
  }
  ## Get data from (self-reported) Parkinson's patients
  dat <- dat[dat$PD == TRUE,]
  ## Keep only the "before medication" and "after medication" data
  ## (drops the "at another time" data, as well as, the data missing 
  ## the medTimepoint response)
  idxBefore <- which(dat$medTimepoint == "Immediately before Parkinson medication")
  idxAfter <- which(dat$medTimepoint == "Just after Parkinson medication (at your best)")
  dat <- dat[sort(c(idxBefore, idxAfter)),]
  ## Make sure that there are only two levels for 
  ## the medTimepoint factor
  dat$medTimepoint <- factor(dat$medTimepoint)
  ## Count the number of activity tasks performed before and
  ## after medication
  countBA <- CountBeforeAfterMedicationPerParticipant(dat)
  ## Select participants that performed at least "beforeThr" 
  ## activity tasks before medication, and at least "afterThr" 
  ## tasks after medication
  tokeep <- GetParticipants(x = countBA, beforeThr = beforeThr, afterThr = afterThr)
  dat <- dat[dat$healthCode %in% tokeep,]  
  ## Replace infinite values by NA
  dat[sapply(dat, is.infinite)] <- NA
  return(dat)
}

#' Retrieve Time of Day Information
#'
#' Process the createdOn variable in order to extract the
#' time-of-the-day that the activity was performed. It takes
#' the local time on createdOn, transform it to UTC, and then
#' back to local time by removing an offset from the UTC time
#' based on the state of residency (time zone data) of the
#' participant.
#'
#' @param dat featurized dataset (dataframe, tibble)
#' @param tzdat dataframe containing state & timezone information
IncludeUTCandLocalTimeVariables <- function(dat, tzDat) {
  GetDayAndTimeOfDay <- function(x) {
    n <- length(x)
    tod <- rep(NA, n)
    day <- rep(NA, n)
    idx <- which(!is.na(x))
    x <- x[!is.na(x)]
    x <- as.character(x)
    n <- length(x)
    aux <- unlist(strsplit(x, " "))
    if (length(aux) == (2*n)) {
      auxDay <- aux[seq(1, 2*n-1, by = 2)]
      auxTime <- aux[seq(2, 2*n, by = 2)]
      auxTime <- as.numeric(unlist(strsplit(auxTime, ":")))
      ih <- seq(1, 3*n, by = 3)
      im <- seq(2, 3*n, by = 3)
      is <- seq(3, 3*n, by = 3)
      secs <- 3600 * auxTime[ih] + 60 * auxTime[im] + auxTime[is]
      tod[idx] <- secs/3600
      day[idx] <- unlist(auxDay)
    }
    if (length(aux) == (3*n)) {
      auxDay <- aux[seq(1, 3*n-2, by = 3)]
      auxTime <- aux[seq(2, 3*n-1, by = 3)]
      auxTime <- as.numeric(unlist(strsplit(auxTime, ":")))
      ih <- seq(1, 3*n, by = 3)
      im <- seq(2, 3*n, by = 3)
      is <- seq(3, 3*n, by = 3)
      secs <- 3600 * auxTime[ih] + 60 * auxTime[im] + auxTime[is]
      tod[idx] <- secs/3600
      day[idx] <- unlist(auxDay)
    }
    list(tod = tod, day = day)
  }
  ## rename the original createdOn time stamp
  ## as the createdOnUTC
  dat$createdOnUTC <- as.POSIXct(dat$createdOn, tz = "UTC")
  ## create day and time of the day variables in UTC
  auxTime <- GetDayAndTimeOfDay(dat$createdOnUTC)
  dat$todUTC <- auxTime$tod
  dat$dayUTC <- as.POSIXct(auxTime$day, tz = "UTC")
  ## get the day and tod values in local time
  ids <- as.character(unique(dat$healthCode))
  aux <- tzDat[match(ids, tzDat$healthCode), c("healthCode", "UTC_offset")]
  utcOffset <- rep(NA, nrow(dat))
  for (i in seq(length(ids))) {
    utcOffset[which(dat$healthCode %in% aux[i, "healthCode"])] <- aux[i, "UTC_offset"]
  }
  tod <- dat$todUTC + utcOffset
  ## negative values correspond to tasks done later in the
  ## day, that show up as done earlier in the next day in
  ## UTC time
  nextday <- which(tod < 0)
  ## fix the local tod
  tod[nextday] <- 24 + tod[nextday]
  ## fix the local day (go one day back)
  day <- dat$dayUTC
  day[nextday] <- day[nextday] - (24*60*60)
  ## include the local time variables
  ## (it still prints as UTC, though)
  dat$tod <- tod
  dat$day <- day
  dat$createdOn <- as.POSIXct(day + tod * 3600)
  return(dat)
}

#' Rank Quantile Transformation
#'
#' Performs a rank-quantile transformation of the features.
#'
#' @param dat featurized dataset (dataframe, tibble)
#' @param featNames list of features
TransformFeatures <- function(dat, featNames) {
  NormalTrans <- function(x) {
    n <- sum(!is.na(x))
    r <- rank(x, na.last = "keep")
    qnorm((r - 0.5)/n)
  }
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids)
  nfeat <- length(featNames)
  for (i in seq(nids)) {
    idx <- which(dat$healthCode == ids[i])
    for (j in seq(nfeat)) {
      dat[idx, featNames[j]] <- NormalTrans(dat[idx, featNames[j]])
    }
  }
  return(dat)
}


#' Loess Regression Feature Detrending
#'
#' Detrend the feature data using loess regression.
#'
#' @param dat featurized dataset (dataframe, tibble)
#' @param featNames list of features
LoessDetrendedFeatures <- function(dat, featNames) {
  participantIds <- unique(dat$healthCode)
  nParticipants <- length(participantIds)
  featColumns <- match(featNames, colnames(dat))
  for (i in seq(nParticipants)) {
    #cat(i, "\n")
    idx1 <- which(dat$healthCode == participantIds[i])
    pdat <- dat[idx1,]
    o <- order(pdat$createdOn)
    pdat <- pdat[o,]
    xaxis <- seq(nrow(pdat))
    for (j in featColumns) {
      aux <- loess(pdat[, j] ~ xaxis)
      pdat[as.numeric(names(aux$residuals)), j] <- aux$residuals + median(pdat[, j], na.rm = TRUE)
    }
    dat[idx1[o],] <- pdat
  }
  return(dat)
}


################################################
## temporal regression functions
################################################

#' Run all user Treatment vs Time of Day N-of-1 Analysis (Newey-West)
#'
#' Run the treatment versus tod personalized analyses
#' based on linear regression and Newey-West regression
#' for all healthCodes in a dataset, and stores the results.
#'
#' @param dat featurized activity dataset (dataframe, tibble)
#' @param featNames list of features
RunTreatmentVsTodEffectsLmNw <- function(dat, featNames) {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids)
  lmPvals <- vector(mode = "list", length = nids)
  names(lmPvals) <- ids
  nwPvals <- lmPvals
  lmBetas <- lmPvals
  lmResAcf <- lmPvals
  for (i in seq(nids)) {
    cat("participant: ", i, "\n")
    aux <- ParticipantTreatVsTodEffectsLmNw(dat, ids[i], featNames)
    lmBetas[[i]] <- aux$lmBetas
    lmPvals[[i]] <- aux$lmPvals
    nwPvals[[i]] <- aux$nwPvals
    lmResAcf[[i]] <- aux$lmResAcf
  }
  return(list(lmBetas = lmBetas,
              lmPvals = lmPvals,
              nwPvals = nwPvals,
              lmResAcf = lmResAcf))
}


#' Run all user Treatment vs Time of Day N-of-1 Analysis (ARIMA)
#'
#' Run the treatment versus tod personalized analyses
#' based on regression with ARIMA errors for all healthCodes
#' in a dataset, and stores the results.
#'
#' @param dat featurized activity dataset (dataframe, tibble)
#' @param featNames list of features
RunTreatmentVsTodEffectsArima <- function(dat, featNames) {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids)
  arimaPvals <- vector(mode = "list", length = nids)
  names(arimaPvals) <- ids
  arimaBetas <- arimaPvals
  arimaResAcf <- arimaPvals
  arimaModel <- arimaPvals
  for (i in seq(nids)) {
    cat("participant: ", i, "\n")
    aux <- ParticipantTreatVsTodEffectsArima(dat, ids[i], featNames)
    arimaBetas[[i]] <- aux$arimaBetas
    arimaPvals[[i]] <- aux$arimaPvals
    arimaResAcf[[i]] <- aux$arimaResAcf
    arimaModel[[i]] <- aux$arimaModel
  }

  return(list(arimaBetas = arimaBetas,
       arimaPvals = arimaPvals,
       arimaResAcf = arimaResAcf,
       arimaModel = arimaModel))
}


#' Run each user Treatment vs Time of Day N-of-1 Analysis (LmNw)
#'
#' Function that actually runs the treatment versus tod
#' personalized analyses based on linear regression and
#' Newey-West regression for a given specific healthCode.
#' @param dat featurized activity dataset (dataframe, tibble)
#' @param id healthcode id
#' @param featNames list of features
ParticipantTreatVsTodEffectsLmNw <- function(dat, id, featNames) {
  
  dat <- dat[dat$healthCode == id,]
  nfeat <- length(featNames)
  lmBetas <- data.frame(matrix(NA, nfeat, 5)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(lmBetas) <- featNames
  colnames(lmBetas) <- c("b.TX", "b.YX", "b.YT", "b.YX|T", "b.YT|X", "feature")
  
  lmPvals <- data.frame(matrix(NA, nfeat, 5)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(lmPvals) <- featNames
  colnames(lmPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0", "feature")
  
  nwPvals <- data.frame(matrix(NA, nfeat, 5)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(nwPvals) <- featNames
  colnames(nwPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0", "feature")
  
  lmResAcf <- data.frame(matrix(NA, nfeat, 4)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(lmResAcf) <- featNames
  colnames(lmResAcf) <- c("lag1_Y~X", "lag1_Y~X+T",  "pval_Y~X", "pval_Y~X+T", "feature")

  x <- dat$medTimepoint
  tod <- dat$tod
  fit1 <- lm(tod ~ x)
  su1 <- summary(fit1)$coefficients
  lmBetas[, "b.TX"] <- su1[2, 1]
  lmPvals[, "a(T,X)=0"] <- su1[2, 4]
  nwPvals[, "a(T,X)=0"] <- coeftest(fit1, vcov = NeweyWest)[2, 4]
  for (i in seq(nfeat)) {
    #cat(i, "\\n")
    y <- dat[, featNames[i]]
    ## lm fits
    fit2 <- lm(y ~ x)
    fit3 <- lm(y ~ tod)
    fit4 <- lm(y ~ x + tod)
    su2 <- summary(fit2)$coefficients
    su3 <- summary(fit3)$coefficients
    su4 <- summary(fit4)$coefficients
    lmBetas[i, "b.YX"] <- su2[2, 1]
    lmBetas[i, "b.YT"] <- su3[2, 1]
    lmBetas[i, "b.YX|T"] <- su4[2, 1]
    lmBetas[i, "b.YT|X"] <- su4[3, 1]
    lmPvals[i, "a(Y,X)=0"] <- su2[2, 4]
    lmPvals[i, "a(Y,T)=0"] <- su3[2, 4]
    lmPvals[i, "a(Y,X|T)=0"] <- su4[2, 4]
    lmPvals[i, "a(Y,T|X)=0"] <- su4[3, 4]
    ## lm with Newey-West covariance estimation
    nwPvals[i, "a(Y,X)=0"] <- coeftest(fit2, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,T)=0"] <- coeftest(fit3, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,X|T)=0"] <- coeftest(fit4, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,T|X)=0"] <- coeftest(fit4, vcov = NeweyWest)[3, 4]

    ## get residuals from lm fits
    res2 <- fit2$resid
    lmResAcf[i, "lag1_Y~X"] <- acf(res2, plot = FALSE)$acf[2]
    lmResAcf[i, "pval_Y~X"] <- Box.test(res2, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
    res4 <- fit4$resid
    lmResAcf[i, "lag1_Y~X+T"] <- acf(res4, plot = FALSE)$acf[2]
    lmResAcf[i, "pval_Y~X+T"] <- Box.test(res4, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
  }

  return(list(lmBetas = lmBetas,
       lmPvals = lmPvals,
       nwPvals = nwPvals,
       lmResAcf = lmResAcf))
}


#' Run each user Treatment vs Time of Day N-of-1 Analysis (LmNw)
#'
#' Function that actually runs the treatment versus tod
#' personalized analyses based on regression with ARIMA errors
#' Newey-West regression for a given specific healthCode.
#'
#' @param dat featurized activity dataset (dataframe, tibble)
#' @param id healthcode id
#' @param featNames list of features
ParticipantTreatVsTodEffectsArima <- function(dat, id, featNames) {
  arimaAnalysisTest1 <- function(y, xreg, lag = 5, ic = "bic", max.p = 10, max.q = 10) {
    aux <- auto.arima(x = y, xreg = xreg, ic = ic, approximation = TRUE, max.p = max.p, max.q = max.q)
    auxOrder <- arimaorder(aux)
    p <- auxOrder[1]
    d <- auxOrder[2]
    q <- auxOrder[3]
    eff <- aux$coef["xreg"]
    zstat <- abs(eff)/sqrt(aux$var.coef["xreg", "xreg"])
    pval <- 2 * pnorm(zstat, lower.tail = FALSE)
    list(eff = eff, pval = pval, p = p, d = d, q = q, autoArimaFit = aux)
  }
  arimaAnalysisTest2 <- function(y, xreg, lag = 5, ic = "bic", max.p = 10, max.q = 10) {
    aux <- auto.arima(x = y, xreg = xreg, ic = ic, approximation = TRUE, max.p = max.p, max.q = max.q)
    auxOrder <- arimaorder(aux)
    p <- auxOrder[1]
    d <- auxOrder[2]
    q <- auxOrder[3]
    ####
    eff1 <- aux$coef["xreg1"]
    zstat1 <- abs(eff1)/sqrt(aux$var.coef["xreg1", "xreg1"])
    pval1 <- 2 * pnorm(zstat1, lower.tail = FALSE)
    ####
    eff2 <- aux$coef["xreg2"]
    zstat2 <- abs(eff2)/sqrt(aux$var.coef["xreg2", "xreg2"])
    pval2 <- 2 * pnorm(zstat2, lower.tail = FALSE)

    list(eff1 = eff1, pval1 = pval1, eff2 = eff2, pval2 = pval2,
         p = p, d = d, q = q, autoArimaFit = aux)
  }
  dat <- dat[dat$healthCode == id,]
  nfeat <- length(featNames)

  arimaBetas <- data.frame(matrix(NA, nfeat, 5)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(arimaBetas) <- featNames
  colnames(arimaBetas) <- c("b.TX", "b.YX", "b.YT", "b.YX|T", "b.YT|X", "feature")
  
  arimaPvals <- data.frame(matrix(NA, nfeat, 5)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(arimaPvals) <- featNames
  colnames(arimaPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0", "feature")
  
  arimaResAcf <- data.frame(matrix(NA, nfeat, 4)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(arimaResAcf) <- featNames
  colnames(arimaResAcf) <- c("lag1_Y~X", "lag1_Y~X+T",  "pval_Y~X", "pval_Y~X+T", "feature")

  arimaModel <- data.frame(matrix(NA, nfeat, 6)) %>% 
    dplyr::mutate(feature = featNames)
  rownames(arimaModel) <- featNames
  colnames(arimaModel) <- c("p_Y~X", "d_Y~X", "q_Y~X", "p_Y~X+T", "d_Y~X+T", "q_Y~X+T", "feature")

  x <- dat$medTimepoint
  tod <- dat$tod

  afit1 <- try(arimaAnalysisTest1(y = tod, xreg = as.numeric(x) - 1,
                                  lag = lag, ic = "bic", max.p = 10, max.q = 10), silent = TRUE)
  if (!inherits(afit1, "try-error")) {
    arimaBetas[, "b.TX"] <- afit1$eff
    arimaPvals[, "a(T,X)=0"] <- afit1$pval
  }

  for (i in seq(nfeat)) {
    #cat(i, "\n")
    y <- dat[, featNames[i]]

    ## arima fits
    afit2 <- try(arimaAnalysisTest1(y = y, xreg = as.numeric(x) - 1), silent = TRUE)
    if (!inherits(afit2, "try-error")) {
      arimaBetas[i, "b.YX"] <- afit2$eff
      arimaPvals[i, "a(Y,X)=0"] <- afit2$pval
      arimaModel[i, "p_Y~X"] <- afit2$p
      arimaModel[i, "d_Y~X"] <- afit2$d
      arimaModel[i, "q_Y~X"] <- afit2$q
      res2 <- afit2$autoArimaFit$resid
      res2 <- res2[!is.na(res2)]
      arimaResAcf[i, "lag1_Y~X"] <- acf(res2, plot = FALSE)$acf[2]
      arimaResAcf[i, "pval_Y~X"] <- Box.test(res2, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
    }
    afit3 <- try(arimaAnalysisTest1(y = y, xreg = tod), silent = TRUE)
    if (!inherits(afit3, "try-error")) {
      arimaBetas[i, "b.YT"] <- afit3$eff
      arimaPvals[i, "a(Y,T)=0"] <- afit3$pval
    }
    xreg <- cbind(as.numeric(x) - 1, tod)
    colnames(xreg) <- NULL
    afit4 <- try(arimaAnalysisTest2(y = y, xreg = xreg), silent = TRUE)
    if (!inherits(afit4, "try-error")) {
      arimaBetas[i, "b.YX|T"] <- afit4$eff1
      arimaBetas[i, "b.YT|X"] <- afit4$eff2
      arimaPvals[i, "a(Y,X|T)=0"] <- afit4$pval1
      arimaPvals[i, "a(Y,T|X)=0"] <- afit4$pval2
      arimaModel[i, "p_Y~X+T"] <- afit4$p
      arimaModel[i, "d_Y~X+T"] <- afit4$d
      arimaModel[i, "q_Y~X+T"] <- afit4$q
      res4 <- afit4$autoArimaFit$resid
      res4 <- res4[!is.na(res4)]
      arimaResAcf[i, "lag1_Y~X+T"] <- acf(res4, plot = FALSE)$acf[2]
      arimaResAcf[i, "pval_Y~X+T"] <- Box.test(res4, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
    }
  }

  return(list(arimaBetas = arimaBetas,
       arimaPvals = arimaPvals,
       arimaResAcf = arimaResAcf,
       arimaModel = arimaModel))
}



################################################
## functions to shape and process the outputs
## of the treatments versus tod analyses
################################################

#' Run Union-Intersection Putative Treatment Effects
#'
#' Computes the union-intersection p-values for putative treatment effects for
#' each participant in the data frame xL.
#' mtMethod1 and mtMethod2 correspond to the methods used for the first and
#' second levels of multiple testing correction.
#'
#' @param xL user temporal regression p-values for the putatite treatment effects
#' @param mtMethod1 first level method for significance correction
#' @param mtMethod2 second level method for significance correction
#' @param thr level of significance (alpha)
UItestsPutativeTreat <- function(xL,
                                 mtMethod1 = "BH",
                                 mtMethod2 = "BH",
                                 thr = 0.05) {
  ids <- names(xL)
  nids <- length(ids)
  nfeat <- nrow(xL[[1]])
  pvals <- matrix(NA, nids, 1)
  rownames(pvals) <- ids
  colnames(pvals) <- "pvalUItest"
  for (i in seq(nids)) {
    pvals[i, 1] <- ParticipantPutativeTreatPval(xL[[i]], mtMethod1, mtMethod2, thr)
  }
  return(pvals)
}



#' Run Union-Intersection Putative Time of Day Effects
#'
#' Computes the union-intersection p-values for putative tod effects for
#' each participant in the data frame xL.
#' mtMethod1 and mtMethod2 correspond to the methods used for the first and
#' second levels of multiple testing correction.
#'
#' @param xL user temporal regression p-values for the putatite tod effects
#' @param mtMethod1 first level method for significance correction
#' @param mtMethod2 second level method for significance correction
#' @param thr level of significance (alpha)
UItestsPutativeTod <- function(xL,
                               mtMethod1 = "none",
                               mtMethod2 = "BH",
                               thr = 0.05) {
  ids <- names(xL)
  nids <- length(ids)
  nfeat <- nrow(xL[[1]])
  pvals <- matrix(NA, nids, 1)
  rownames(pvals) <- ids
  colnames(pvals) <- "pvalUItest"
  for (i in seq(nids)) {
    pvals[i, 1] <- ParticipantPutativeTodPval(xL[[i]], mtMethod1, mtMethod2, thr)
  }
  return(pvals)
}


#' Run adjustment for p-values for Treatment Effects
#'
#' Given a matrix x of temporal regression p-values for the putatite
#' treatment effects of a given participant (where the rows index the
#' features and columns index the conditional independence tests), this
#' function determines for each feature if we need to adjust for tod,
#' then selects the p-values from either the marginal/unadjusted or
#' conditional/tod-adjusted p-values, which will be used for the
#' computation of the union-intersection p-value.
##
#' mtMethod1 corresponds to the method used for the first level
#' of multiple testing correction (where for each feature we adjust
#' across the 5 conditional independence tests).
#'
#' mtMethod2 corresponds to the method used for the second level
#' of multiple testing correction (where we adjust across all
#' features).
#'
#' @param x user temporal regression p-values for the putatite treatment effects
#' @param mtMethod1 first level method for significance correction
#' @param mtMethod2 second level method for significance correction
#' @param thr level of significance (alpha)

ParticipantPutativeTreatPval <- function(x, mtMethod1 = "BH",
                                         mtMethod2 = "BH",
                                         thr = 0.05) {
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
  return(min(pvals)) ## union-intersection p-value
}


#' Run adjustment for p-values for Time of Day Effects
#'
#' Given a matrix x of temporal regression p-values for the putatite
#' tod effects of a given participant (where the rows index the
#' features and columns index the conditional independence tests), this
#' function determines for each feature if we need to adjust for treatment,
#' then selects the p-values from either the marginal/unadjusted or
#' conditional/treatment-adjusted p-values, which will be used for the
#' computation of the union-intersection p-value.
#'
#' mtMethod1 corresponds to the method used for the first level
#' of multiple testing correction (where for each feature we adjust
#' across the 5 conditional independence tests).
#'
#' mtMethod2 corresponds to the method used for the second level
#' of multiple testing correction (where we adjust across all
#' features).
#'
#' @param x user temporal regression p-values for the putatite treatment effects
#' @param mtMethod1 first level method for significance correction
#' @param mtMethod2 second level method for significance correction
#' @param thr level of significance (alpha)
ParticipantPutativeTodPval <- function(x, mtMethod1 = "none",
                                       mtMethod2 = "BH",
                                       thr = 0.05) {
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
  return(min(pvals))
}


#' Participant Evidence
#'
#' Determines if a feature shows evidence of treat, tod, or
#' both treat and tod putative effects looking at the conditional
#' independence patterns quantified by the p-values from the temporal
#' regression models.
#'
#' @param x user temporal regression p-values for the putatite treatment effects
#' @param mtMethod method for significance correction
#' @param thr level of significance (alpha)
ParticipantEvidence <- function(x, mtMethod = "BH", thr = 0.05) {
  PvalAdjustment <- function(x, mtMethod, byrow = TRUE) {
    if (byrow) {
      ax <- t(apply(x, 1, p.adjust, mtMethod))
    }
    else {
      ax <- apply(x, 2, p.adjust, mtMethod)
    }
    return(ax)
  }
  SummarizeEvidence <- function(x, thr = 0.05) {
    EvidenceRule <- function(pvals, thr) {
      ## 1: a(T,X)=0, 2: a(Y,X)=0, 3: a(Y,T)=0, 4: a(Y,X|T)=0, 5: a(Y,T|X)=0
      evidence <- NA
      if (sum(is.na(pvals)) == 0) {
        if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] <= thr & pvals[5] <= thr) {
          evidence <- "treat and tod effects"
        }
        if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] <= thr & pvals[5] > thr) {
          evidence <- "treat effect"
        }
        if(pvals[2] <= thr & pvals[3] > thr & pvals[4] <= thr & pvals[5] > thr) {
          evidence <- "treat effect"
        }
        if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] <= thr) {
          evidence <- "tod effect"
        }
        if(pvals[2] > thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] <= thr) {
          evidence <- "tod effect"
        }
        if(pvals[2] <= thr & pvals[3] > thr & pvals[4] > thr & pvals[5] > thr) {
          evidence <- "inconsistent treat effect"
        }
        if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] > thr) {
          evidence <- "inconsistent treat and tod effect"
        }
        if(pvals[2] > thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] > thr) {
          evidence <- "inconsistent tod effect"
        }
      }
      return(evidence)
    }
    evidence <- apply(x, 1, EvidenceRule, thr)
    return(data.frame(x, evidence, check.names = FALSE))
  }
  ## performs multiple testing adjustment across the 5
  ## conditional independence tests for each feature
  ## separately
  ax <- PvalAdjustment(x, mtMethod, byrow = TRUE)
  ## compute the evidence on the adjusted p-values
  ev <- as.character(SummarizeEvidence(ax, thr)[, "evidence"])
  return(list(ev = ev))
}


#' Merge demographic variables with the union-intersection p-values
#'
#' @param ui union intersection matrix of users
#' @param demo demographics data
#' @param mtmethod p-values correction method
CollateUipvalsAndDemo <- function(ui, demo, mtmethod = "none") {
  colnames(ui) <- c("uiLm", "uiNw", "uiAr")
  ids1 <- rownames(ui)
  ids2 <- demo$healthCode
  ids <- intersect(ids1, ids2)
  dat <- cbind(ui[match(ids, rownames(ui)),], demo[match(ids, demo$healthCode),])
  dat$medication.start.year[which(dat$medication.start.year == 0)] <- NA
  dat$logPvalUiLm <- -log(p.adjust(dat$uiLm, method = mtmethod), 10)
  dat$logPvalUiNw <- -log(p.adjust(dat$uiNw, method = mtmethod), 10)
  dat$logPvalUiAr <- -log(p.adjust(dat$uiAr, method = mtmethod), 10)
  return(dat)
}


#' Fit regression lines to the p-values of responders and non-responders
#' @param dat featurize activities (dataframe tibble)
#' @param respName responder names (healthCode)
#' @param covName covariate names
#' @param logPvalName name log P-value threshold
#' @param logPvalThr threshold for log P-value
RespondersAndNonRespondersFits <- function(dat,
                                           respName,
                                           covName,
                                           logPvalName,
                                           logPvalThr) {
  ## Replace infinite values by NA
  dat[sapply(dat, is.infinite)] <- NA
  form <- as.formula(paste(respName, covName, sep = " ~ "))
  ## all
  fit1 <- lm(form, data = dat)
  a1 <- fit1$coef[1]
  b1 <- fit1$coef[2]
  pval1 <- summary(fit1)$coefficients[2, 4]
  ## responders
  idx2 <- which(dat[, logPvalName] > logPvalThr)
  dat2 <- dat[idx2,]
  fit2 <- lm(form, data = dat2)
  a2 <- fit2$coef[1]
  b2 <- fit2$coef[2]
  pval2 <- summary(fit2)$coefficients[2, 4]
  ## non-responders
  idx3 <- which(dat[, logPvalName] <= logPvalThr)
  dat3 <- dat[idx3,]
  fit3 <- lm(form, data = dat3)
  a3 <- fit3$coef[1]
  b3 <- fit3$coef[2]
  pval3 <- summary(fit3)$coefficients[2, 4]
  ##
  responder <- rep(NA, nrow(dat))
  responder[idx2] <- TRUE
  responder[idx3] <- FALSE
  dat$responder <- responder

  return(list(dat = dat, form = form,
       respName = respName,
       covName = covName,
       a = a1, b = b1, pval = pval1,
       aR = a2, bR = b2, pvalR = pval2,
       aNR = a3, bNR = b3, pvalNR = pval3))
}


#' Merge Activity Outputs
#' @param oTap output matrix for tap
#' @param oVoi output matrix for voice
#' @param oWal output matrix for walk
#' @param oRes output matrix for rest
#' @param colName designated column for merging ("uiTreLm", "uiTreAr", etc.)
MergeOutputs <- function(oTap, oVoi, oWal, oRes, colName) {
  aux <- unique(c(rownames(oTap), rownames(oVoi), rownames(oWal), rownames(oRes)))
  all <- matrix(NA, length(aux), 4)
  rownames(all) <- aux
  colnames(all) <- c("tapping", "voice", "walk", "rest")
  all[match(rownames(oTap), aux), "tapping"] <- oTap[, colName]
  all[match(rownames(oVoi), aux), "voice"] <- oVoi[, colName]
  all[match(rownames(oWal), aux), "walk"] <- oWal[, colName]
  all[match(rownames(oRes), aux), "rest"] <- oRes[, colName]
  return(all)
}


#' Get Participants Proportion in Treatment & Tod or individual proportion
#'
#' Get proportions of participants showing treatment, tod, as well as,
#' treatment alone, tod alone, and treatment and tod together.
#'
#' @param uiTre union-intersection of treatment effect matrix
#' @param uiTod union-intersection of time of day effect matrix
GetProportions <- function(uiTre, uiTod, thr = 0.05, r = 2, task = NULL) {
  out <- matrix(NA, 5, 1)
  rownames(out) <- c("treat", "tod", "treat_alone", "tod_alone", "treat_and_tod")
  colnames(out) <- task
  n <- nrow(uiTre)
  out[1, 1] <- round(sum(uiTre[, 2] <= thr)/n, r)
  out[2, 1] <- round(sum(uiTod[, 2] <= thr)/n, r)
  out[3, 1] <- round(sum(uiTre[, 2] <= thr & uiTod[, 2] > thr)/n, r)
  out[4, 1] <- round(sum(uiTod[, 2] <= thr & uiTre[, 2] > thr)/n, r)
  out[5, 1] <- round(sum(uiTre[, 2] <= thr & uiTod[, 2] <= thr)/n, r)
  return(out)
}


GetRelativeImportance <- function(sdat, featNames) {
  nfeat <- length(featNames)
  out <- matrix(NA, nfeat, 3) 
  rownames(out) <- featNames
  colnames(out) <- c("R2", "medTimepoint", "tod")
  for (i in seq(nfeat)) {
    myformula <- as.formula(paste(featNames[i], "~ medTimepoint + tod"))
    fit <- lm(myformula, data = sdat)
    aux <- calc.relimp(fit, type = "lmg", rela = FALSE)
    out[i, 1] <- aux$R2
    out[i, 2:3] <- aux$lmg
  }
  return(out %>% 
           tibble::as_tibble(.) %>%
           dplyr::mutate(feature = featNames))
}


RelativeImportances <- function(dat, featNames) {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids) 
  relImp <- vector(mode = "list", length = nids)
  names(relImp) <- ids
  for (i in seq(nids)) {
    cat("participant: ", i, "\n")
    sdat <- dat[dat$healthCode == ids[i],]
    relImp[[i]] <- GetRelativeImportance(sdat, featNames)
  }
  return(relImp)
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
  return(out)
}


##########################################################
## Functions used for population analysis which revolves 
## around confounding assessment, case vs controls assessment,
## and feature variability comparison assessment
############################################################

# remove imports later
library(pROC)
library(randomForest)
library(ppcor)
library(glmnet)
library(energy)
library(lubridate)


#' Pipeline for PD case vs controls matching
#' 
#' This function will match healthcodes who contributed more than 5 records 
#' by separating male/female demographics and then match their distribution by pairing 
#' users with the other same age user in each of the gender cohort to reduce demographic
#' confounding towards Parkinson diagnosis
#' 
#' @param data activity features (with demographic information of age, education, gender, diagnosis)
#' @param features sensor features to choose from dataframe
#' @param plot (boolean) plot or no plot
#' 
#' @examples 
#' PD_case_vs_controls_matching(
#'     tap_features, c("tapInterval,...etc."), TRUE 
#'
#' @return 
#' A demographic dataframe of matched healthcodes
PD_case_vs_controls_matching  <- function(dat, features, threshold, plot = FALSE){
  dat$gender <- as.character(dat$gender)
  aux <- GetCollapsedData(x = dat, 
                          labelName = "PD", 
                          covNames = c("gender", "age"), 
                          subjectIdName = "healthCode", 
                          featNames = features)
  cdat <- aux$out
  
  ## include new column with the number of records counts
  nrecs <- table(dat$healthCode)
  cdat$nrecs <- nrecs[cdat$healthCode]
  
  ## keep of participants that contributed at least 5 records
  cdatf <- cdat[cdat$nrecs >= threshold,]
  
  ## perform the matching
  maux <- SeparateMaleFemaleAgeMatching(cdatf)
  
  if(plot == TRUE){
    boxplot(age ~ pd, data = maux)
    boxplot(age ~ pd + gender, data = maux)
  }
  return(maux)
}



#######################################################
## Helper Functions
#######################################################
#' Set Education into Ranked Variables
#'
#' Function that takes in a vector of education and change it into rankings
#' @param x education vector
#' @return ranked education vector
#' @examples
#' tap.data$education <- NumericEducation(tap.data$education)
NumericEducation <- function(x) {
  xx <- rep(NA, length(x))
  xx[x == "Some high school"] <- 1
  xx[x == "School Diploma/GED"] <- 2
  xx[x == "Some college"] <- 3
  xx[x == "2-year college degree"] <- 3
  xx[x == "4-year college degree"] <- 4
  xx[x == "Some graduate school"] <- 5
  xx[x == "Master's Degree"] <- 6
  xx[x == "Doctoral Degree"] <- 7
  return(xx)
}


#' Change Gender to Binary
#'
#' Function that takes in a vector of gender and change into binary variable
#' @param x gender vector
#' @return binary gender vector
#' @examples
#' tap.data$gender <- BinaryGender(tap.data$gender)
BinaryGender <- function(x) {
  xx <- rep(NA, length(x))
  xx[x == "Female"] <- 1
  xx[x == "Male"] <- 0
  return(xx)
}


#' Aggregate Features
#'
#' Function that takes in dataframe, labels, covariate names, healthcodes
#' and collapse it into aggregated features that will be used for population analysis
#'
#' @param x featurized dataset
#' @param labelName name of desired label (in our case PD)
#' @param covNames name of covariates  (demographics information)
#' @param subjectIdName healthcode of user
#' @param featNames list of features
#'
#' @return cleaned aggregated dataset
#' @examples
#' CollapseFeatures(x = matched.dataset, \
#' labelName = "PD", \
#' covNames = c("age", "gender2", "education2"), \
#' subjectIdName = "healthCode", \
#' featNames = features)
CollapseFeatures <- function(x, labelName, covNames, subjectIdName, featNames) {
  ids <- as.character(unique(x[, subjectIdName]))
  nids <- length(ids)
  nfeat <- length(featNames)
  nvar <- length(c(subjectIdName, labelName, covNames))
  out <- data.frame(matrix(NA, nids, 2 * nfeat + nvar))
  cfeatNames <- c(paste(featNames, "med", sep = "."), paste(featNames, "iqr", sep = "."))
  colnames(out) <- c(subjectIdName, labelName, covNames, cfeatNames)
  rownames(out) <- ids
  for (i in seq(nids)) {
    sdat <- x[which(x[, subjectIdName] == ids[i]),]
    out[i, subjectIdName] <- ids[i]
    out[i, labelName] <- as.character(sdat[1, labelName])
    out[i, covNames] <- as.character(sdat[1, covNames])
    out[i, (nvar + 1):(nfeat + nvar)] <- apply(sdat[, featNames], 2, median, na.rm = TRUE)
    out[i, (nfeat + nvar + 1):(2 * nfeat + nvar)] <- apply(sdat[, featNames], 2, IQR, na.rm = TRUE)
  }
  for (j in seq(length(covNames))) {
    out[which(out[, covNames[j]] == "NA"), covNames[j]] <- NA
    out[, covNames[j]] <- as.numeric(out[, covNames[j]])
  }
  return(list(out = out, cfeatNames = cfeatNames))
}


#' Filter Record by Threshold
#'
#' Function that takes in aggregated dataframe and filter user with record less than threshold
#'
#' @param dat aggregated features
#' @param thr threshold of records
#' @return Record Filtered Dataset
FilterOutParticipantsWithFewRecords <- function(dat, thr) {
  aux <- table(dat$healthCode)
  aux <- sort(aux, decreasing = TRUE)
  keep <- names(aux)[which(aux >= thr)]
  dat <- dat[dat$healthCode %in% keep,]
  dat$healthCode <- factor(dat$healthCode)
  return(dat)
}


#' Run Random Forest Classifications
#'
#' Run random forest classification
#' train and evaluate the random forest classifier
#' based on permutation and non-permutation
#' on multiple train/test data splits AUC,
#' accuracy and balanced accuracy
#'
#' @param nRuns number of runs for data splitting
#' @param dat aggregated features
#' @param respName metadata features
#' @param featNames activity features
#' @param negClassName class name for negative diagnosis
#' @param posClassName class name for positive diagnosis
#' @return A list of metrics of random forest classification performance metrics
RunRandomForestsC <- function(nRuns,
                              dat,
                              respName,
                              featNames,
                              nSplits = 2,
                              negClassName,
                              posClassName,
                              seeds = NULL,
                              ntrees = 500) {
  ComputeBinaryAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    sum(yhat == ytest)/length(yhat)
  }
  ComputeBalancedAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    tp <- sum(ytest == 1 & yhat == 1)
    pos <- sum(ytest == 1)
    tn <- sum(ytest == 0 & yhat == 0)
    neg <- sum(ytest == 0)
    bacc <- (tp/pos + tn/neg)/2
    bacc
  }
  TrainTestRandomSplitC <- function(dat, nSplits = 2) {
    n <- nrow(dat)
    test <- sample(seq(n), round(n/nSplits), replace = FALSE)
    train <- setdiff(seq(n), test)

    list(trainDat = dat[train,], testDat = dat[test,])
  }
  dat[, respName] <- factor(as.character(dat[, respName]),
                            levels = c(negClassName, posClassName))
  dat <- dat[, c("healthCode", respName, featNames)]
  dat <- na.omit(dat)
  metrics <- matrix(NA, nRuns, 3)
  colnames(metrics) <- c("Auc", "Acc", "BalAcc")
  metricsP <- matrix(NA, nRuns, 3)
  colnames(metricsP) <- c("Auc", "Acc", "BalAcc")
  nFeat <- length(featNames)
  imp <- matrix(NA, nFeat, nRuns)
  rownames(imp) <- colnames(dat)[-c(1, 2)]
  impP <- matrix(NA, nFeat, nRuns)
  rownames(impP) <- rownames(imp)
  pvals <- matrix(NA, nRuns, 1)
  colnames(pvals) <- c("pval")
  myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
  for (i in seq(nRuns)) {
    cat("run ", i, "\n")
    if (!is.null(seeds)) {
      set.seed(seeds[i])
    }
    splitDat <- TrainTestRandomSplitC(dat, nSplits)
    trainDat <- splitDat$trainDat[, -1]
    testDat <- splitDat$testDat[, -1]
    ## random forest fit
    fit <- randomForest::randomForest(myFormula, data = trainDat, ntree = ntrees)
    if (!inherits(fit, "try-error")) {
      imp[, i] <- fit$importance[, 1]
      predProbs <- predict(fit, testDat[, -1, drop = FALSE], type = "prob")
      ytest <- testDat[, 1]
      rocb <- pROC::roc(ytest, predProbs[, posClassName], direction = "<",
                  levels = c(negClassName, posClassName))
      metrics[i, "Auc"] <- pROC::auc(rocb)[1]
      metrics[i, "Acc"] <- ComputeBinaryAccuracy(ytest, predProbs[, posClassName], thr = 0.5)
      metrics[i, "BalAcc"] <- ComputeBalancedAccuracy(ytest, predProbs[, posClassName], thr = 0.5)
      pvals[i, "pval"] <- TestAUC(metrics[i, "Auc"], ytest, negClassName, posClassName)
    }
    ## random forest fit to permutated label data
    trainDatP <- trainDat
    trainDatP[, 1] <- trainDatP[sample(nrow(trainDatP), replace = FALSE), 1]
    fitP <- try(randomForest::randomForest(myFormula, data = trainDatP), silent = TRUE)
    if (!inherits(fitP, "try-error")) {
      impP[, i] <- fitP$importance[, 1]
      predProbsP <- predict(fitP, testDat[, -1, drop = FALSE], type = "prob")
      ytest <- testDat[sample(nrow(testDat), replace = FALSE), 1]
      rocbP <- pROC::roc(ytest, predProbsP[, posClassName], direction = "<",
                   levels = c(negClassName, posClassName))
      metricsP[i, "Auc"] <- pROC::auc(rocbP)[1]
      metricsP[i, "Acc"] <- ComputeBinaryAccuracy(ytest, predProbsP[, posClassName], thr = 0.5)
      metricsP[i, "BalAcc"] <- ComputeBalancedAccuracy(ytest, predProbsP[, posClassName], thr = 0.5)
    }
  }
  return(list(metrics = data.frame(metrics, check.names = FALSE),
              metricsP = data.frame(metricsP, check.names = FALSE),
              imp = data.frame(imp, check.names = FALSE) %>% 
                dplyr::mutate(feature = rownames(.)), 
              impP = data.frame(impP, check.names = FALSE) %>%
                dplyr::mutate(feature = rownames(.)),
              pvals = data.frame(pvals, check.names = FALSE)))
}


#' Run Ridge Regression Classifications
#'
#' Run random forest classification
#' train and evaluate the random forest classifier
#' based on permutation and non-permutation
#' on multiple train/test data splits AUC,
#' accuracy and balanced accuracy
#'
#' @param nRuns number of runs for data splitting
#' @param dat aggregated features
#' @param respName target label variable
#' @param featNames activity features
#' @param negClassName class name for negative diagnosis (FALSE)
#' @param posClassName class name for positive diagnosis (TRUE)
#' @return A list of metrics of random forest classification performance metrics
RunRidgeRegrC <- function(nRuns,
                          dat,
                          respName,
                          featNames,
                          nSplits = 2,
                          negClassName,
                          posClassName,
                          seeds = NULL) {
  ComputeBinaryAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    sum(yhat == ytest)/length(yhat)
  }
  ComputeBalancedAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    tp <- sum(ytest == 1 & yhat == 1)
    pos <- sum(ytest == 1)
    tn <- sum(ytest == 0 & yhat == 0)
    neg <- sum(ytest == 0)
    bacc <- (tp/pos + tn/neg)/2
    bacc
  }
  TrainTestRandomSplitC <- function(dat, nSplits = 2) {
    n <- nrow(dat)
    test <- sample(seq(n), round(n/nSplits), replace = FALSE)
    train <- setdiff(seq(n), test)

    list(trainDat = dat[train,], testDat = dat[test,])
  }
  dat[, respName] <- factor(as.character(dat[, respName]),
                            levels = c(negClassName, posClassName))
  dat <- dat[, c("healthCode", respName, featNames)]
  dat <- na.omit(dat)
  metrics <- matrix(NA, nRuns, 3)
  colnames(metrics) <- c("Auc", "Acc", "BalAcc")
  metricsP <- matrix(NA, nRuns, 3)
  colnames(metricsP) <- c("Auc", "Acc", "BalAcc")
  nFeat <- length(featNames)
  imp <- matrix(NA, nFeat, nRuns)
  rownames(imp) <- colnames(dat)[-c(1, 2)]
  impP <- matrix(NA, nFeat, nRuns)
  rownames(impP) <- rownames(imp)
  pvals <- matrix(NA, nRuns, 1)
  colnames(pvals) <- c("pval")
  myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
  for (i in seq(nRuns)) {
    cat("run ", i, "\n")
    if (!is.null(seeds)) {
      set.seed(seeds[i])
    }
    splitDat <- TrainTestRandomSplitC(dat, nSplits)
    trainDat <- splitDat$trainDat[, -1]
    testDat <- splitDat$testDat[, -1]
    xtrain <- scale(trainDat[, featNames])
    xtest <- scale(testDat[, featNames])
    ytrain <- trainDat[, respName]
    ytest <- testDat[, respName]

    ## ridge regr fit
    fit <- try(glmnet::cv.glmnet(x = xtrain, y = ytrain, family = "binomial",
                         alpha = 0), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      predProbs <- predict(fit, newx = xtest, type = "response", s = "lambda.min")
      rocObj <- pROC::roc(ytest, predProbs[, 1], direction = "<",
                    levels = c(negClassName, posClassName))
      metrics[i, "Auc"] <- pROC::auc(rocObj)[1]
      metrics[i, "Acc"] <- ComputeBinaryAccuracy(ytest, predProbs[, 1], thr = 0.5)
      metrics[i, "BalAcc"] <- ComputeBalancedAccuracy(ytest, predProbs[, 1], thr = 0.5)
      pvals[i, "pval"] <- TestAUC(metrics[i, "Auc"], ytest, negClassName, posClassName)
    }
    ## ridge regr fit to permutated label data
    ytrainP <- ytrain[sample(length(ytrain), replace = FALSE)]
    ytestP <- ytest[sample(length(ytest), replace = FALSE)]
    fitP <- try(glmnet::cv.glmnet(x = xtrain, y = ytrainP, family = "binomial",
                          alpha = 0), silent = TRUE)
    if (!inherits(fitP, "try-error")) {
      predProbsP <- predict(fitP, newx = xtest, type = "response", s = "lambda.min")
      rocObjP <- pROC::roc(ytestP, predProbsP[, 1], direction = "<",
                     levels = c(negClassName, posClassName))
      metricsP[i, "Auc"] <- pROC::auc(rocObjP)[1]
      metricsP[i, "Acc"] <- ComputeBinaryAccuracy(ytestP, predProbsP[, 1], thr = 0.5)
      metricsP[i, "BalAcc"] <- ComputeBalancedAccuracy(ytestP, predProbsP[, 1], thr = 0.5)
    }
  }
  return(list(metrics = data.frame(metrics, check.names = FALSE),
              metricsP = data.frame(metricsP, check.names = FALSE),
              pvals = data.frame(pvals, check.names = FALSE)))
}

#' Run Classification Assessment
#'
#' Entrypoint for random forest classsification and ridge regression classification function.
#' Run analyses on 3 different set of features (sensor features alone, demographics alone, and both)
#' for both random forest and ridge regression, assessing AUC, balanced accuracy and accuracy.
#'
#' @param nRuns number of runs for data splitting
#' @param featNames2 feature set 1 (sensor only feature)
#' @param featNames3 feature set 2 (sensor plus demographics)
#' @param featNames4 feature set 3 (only demographcis)
#' @param dat featurized dataset
#' @param respName column of target variable
#' @param nSplits number of split for training and testing
#' @param negClassName name of negative class name in binary target variable
#' @param posClassName name of positive class name in binary target variable
#'
#' @return
#' A list of result metrics for each analysis assessment
#'
#' @examples
#' result <- RunAnalyses(nRuns = nRuns, c('numberTaps'),
#' c('numberTaps', 'age', 'gender2', 'education2'),
#' c('age', 'gender2', 'education2'),
#' fixedSeed = 123,
#' dat = matched.tap.data,
#' respName = 'PD',
#' nSplits = 2,
#' negClassName = FALSE,
#' posClassName = TRUE)
RunAnalyses <- function(nRuns,
                        featNames2,
                        featNames3,
                        featNames4,
                        fixedSeed = 123,
                        dat,
                        respName = "PD",
                        nSplits = 2,
                        negClassName = "FALSE",
                        posClassName = "TRUE") {

  cat("rf, 2", "\n")
  set.seed(fixedSeed)
  rf2 <- RunRandomForestsC(nRuns = nRuns,
                           dat = dat,
                           respName = respName,
                           featNames = featNames2,
                           nSplits = nSplits,
                           negClassName = negClassName,
                           posClassName = posClassName)

  cat("rf, 3", "\n")
  set.seed(fixedSeed)
  rf3 <- RunRandomForestsC(nRuns = nRuns,
                           dat = dat,
                           respName = respName,
                           featNames = featNames3,
                           nSplits = nSplits,
                           negClassName = negClassName,
                           posClassName = posClassName)

  cat("rf, 4", "\n")
  set.seed(fixedSeed)
  rf4 <- RunRandomForestsC(nRuns = nRuns,
                           dat = dat,
                           respName = respName,
                           featNames = featNames4,
                           nSplits = nSplits,
                           negClassName = negClassName,
                           posClassName = posClassName)

  cat("rr, 2", "\n")
  set.seed(fixedSeed)
  rr2 <- RunRidgeRegrC(nRuns = nRuns,
                       dat = dat,
                       respName = respName,
                       featNames = featNames2,
                       nSplits = nSplits,
                       negClassName = negClassName,
                       posClassName = posClassName)

  cat("rr, 3", "\n")
  set.seed(fixedSeed)
  rr3 <- RunRidgeRegrC(nRuns = nRuns,
                       dat = dat,
                       respName = respName,
                       featNames = featNames3,
                       nSplits = nSplits,
                       negClassName = negClassName,
                       posClassName = posClassName)

  cat("rr, 4", "\n")
  set.seed(fixedSeed)
  rr4 <- RunRidgeRegrC(nRuns = nRuns,
                       dat = dat,
                       respName = respName,
                       featNames = featNames4,
                       nSplits = nSplits,
                       negClassName = negClassName,
                       posClassName = posClassName)

  return(list(rf2 = rf2,
              rf3 = rf3,
              rf4 = rf4,
              rr2 = rr2,
              rr3 = rr3,
              rr4 = rr4))
}


#' AUC Check Assessment
#'
#' Analytical test for checking whether AUC is better than random results
#'
#' @param aucObs observed AUC values from permutation/normal assessment (from matrix)
#' @param ytest test y variables for assessment
#' @param negClassName name of negative class name in binary target variable
#' @param posClassName name of positive class name in binary target variable
TestAUC <- function(aucObs, ytest, negClassName, posClassName) {
  GetNormApproxVarAUC <- function(ytest, negClassName, posClassName) {
    ytest <- factor(ytest)
    ylevels <- levels(ytest)
    n1 <- sum(ytest == negClassName)
    n2 <- sum(ytest == posClassName)
    n <- n1 + n2
    v <- (n + 1)/(12 * n1 * n2)
    c(v = v, n = n, nNeg = n1, nPos = n2)
  }
  v <- GetNormApproxVarAUC(ytest, negClassName, posClassName)["v"]
  return(pnorm(aucObs, 0.5, sqrt(v), lower.tail = FALSE))
}


#' Collapse and Shape Data
#'
#' Collapse and clean the data (enforce that label is a factor
#' and return data frame containing only the healthCode,
#' label, covariates and feature variables).
#'
#' @param dat featurized dataset
#' @param labelName name of desired label (in our case PD)
#' @param covNames name of covariates  (demographics information)
#' @param subjectIdName healthcode of user
#' @param featNames list of features
#' @examples
#' CollapseAndShapeData(dat = tap.data.filtered,
#' labelName = "PD",
#' subjectIdName = "healthCode",
#' covNames = c("age", "gender2", "education2"),
#' featNames = tap.features,
#' negClassName = FALSE,
#' posClassName = TRUE)
CollapseAndShapeData <- function(dat,
                                 labelName = "PD",
                                 subjectIdName = "healthCode",
                                 covNames,
                                 featNames,
                                 negClassName = "FALSE",
                                 posClassName = "TRUE") {
  cat("collapsing the features", "\n")
  aux <- CollapseFeatures(x = dat,
                          labelName,
                          covNames,
                          subjectIdName,
                          featNames)
  cat("cleaning the data", "\n")
  dat <- aux$out
  featNamesC <- aux$cfeatNames
  dat[, labelName] <- factor(as.character(dat[, labelName]),
                             levels = c(negClassName, posClassName))
  dat <- dat[, c(subjectIdName, labelName, covNames, featNamesC)]
  dat <- na.omit(dat)
  list(dat = dat, featNamesC = featNamesC)
}


#' Runs correlation and partial correlation tests.
CorTests <- function(dat, labelName, confName, scoreName) {
  cors <- matrix(NA, 5, 1)
  colnames(cors) <- "estimate"
  rownames(cors) <- c(paste("cor(", scoreName, " , ", labelName, ")", sep = ""),
                      paste("cor(", scoreName, " , ", confName, ")", sep = ""),
                      paste("cor(", labelName, " , ", confName, ")", sep = ""),
                      paste("cor(", paste(scoreName, labelName, sep = " , "), " | ", confName, ")", sep = ""),
                      paste("cor(", paste(scoreName, confName, sep = " , "), " | ", labelName, ")", sep = ""))
  pvals <- cors
  colnames(pvals) <- "pval"

  aux1 <- cor.test(dat[, scoreName], dat[, labelName], method = "spearman", exact = FALSE)
  aux2 <- cor.test(dat[, scoreName], dat[, confName], method = "spearman", exact = FALSE)
  aux3 <- cor.test(dat[, labelName], dat[, confName], method = "spearman", exact = FALSE)
  aux4 <- ppcor::pcor.test(dat[, scoreName], dat[, labelName], dat[, confName], method = "spearman")
  aux5 <- ppcor::pcor.test(dat[, scoreName], dat[, confName], dat[, labelName], method = "spearman")

  cors[1, 1] <- aux1$estimate
  cors[2, 1] <- aux2$estimate
  cors[3, 1] <- aux3$estimate
  cors[4, 1] <- aux4$estimate
  cors[5, 1] <- aux5$estimate

  pvals[1, 1] <- aux1$p.value
  pvals[2, 1] <- aux2$p.value
  pvals[3, 1] <- aux3$p.value
  pvals[4, 1] <- aux4$p.value
  pvals[5, 1] <- aux5$p.value

  list(cors = cors, pvals = pvals)
}


#' Run Causality Tests based on Correlation (Random Forest)
#'
#' Runs the causality-based test for confounding
#' for the random forest classifiers
#'
#' @param dat unmatched data (dataframe, tibble)
#' @param datM matched data (dataframe, tibble)
#' @param nruns number of different runs (train/test splitting)
#' @param splitseeds random seed for data splitting
#' @param labelName name of prediction label
#' @param featNames name of desired feature names
#' @param confName name of confounding variable
#' @param negClassName boolean negative for diagnosis
#' @param posClassName boolean positive for diagnosis
#'
#' @examples
#' RunCondIndepTestsCorRf(
#' dat = tap.data,
#' datM = matched.tap.data,
#' nruns = 100, 123456,
#' labelName = "PD",
#' featNames = featNamesC,
#' confName = "age",
#' negClassName = FALSE,
#' posClassName = TRUE)
RunCondIndepTestsCorRf <- function(dat,
                                   datM,
                                   nruns,
                                   splitSeeds,
                                   labelName,
                                   featNames,
                                   confName,
                                   negClassName,
                                   posClassName) {
  GetRfAUC <- function(dat,
                       idxTrain,
                       idxTest,
                       labelName,
                       featNames,
                       negClassName,
                       posClassName) {
    dat <- dat[, c(labelName, featNames)]
    dat[, labelName] <- factor(as.character(dat[, labelName]),
                               levels = c(negClassName, posClassName))
    myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
    fit <- randomForest::randomForest(myFormula, data = dat[idxTrain,], ntree = 1000)
    predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "prob")
    rocObj <- pROC::roc(dat[idxTest, 1], predProbs[, posClassName], direction = "<",
                  levels = c(negClassName, posClassName))
    aucObs <- pROC::auc(rocObj)[1]
    list(aucObs = aucObs, predProbs = predProbs[, posClassName], rocObj = rocObj)
  }
  TrainTestRandomSplitIndex <- function(dat,
                                        nSplits = 2,
                                        respName) {
    aux <- levels(dat[, respName])
    negativeClassName <- aux[1]
    positiveClassName <- aux[2]
    neg <- which(dat[, respName] == negativeClassName)
    pos <- which(dat[, respName] == positiveClassName)
    nNeg <- length(neg)
    nPos <- length(pos)
    testNeg <- sample(neg, round(nNeg/nSplits), replace = FALSE)
    trainNeg <- setdiff(neg, testNeg)
    testPos <- sample(pos, round(nPos/nSplits), replace = FALSE)
    trainPos <- setdiff(pos, testPos)
    train <- c(trainNeg, trainPos)
    test <- c(testNeg, testPos)
    list(itrain = train, itest = test)
  }

  aucs <- matrix(NA, nruns, 2)
  colnames(aucs) <- c("orig", "matching")

  corRf <- matrix(NA, nruns, 10)
  colnames(corRf) <- c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)",
                       "pval(R,Y)", "pval(R,A)", "pval(A,Y)", "pval(R,Y|A)", "pval(R,A|Y)")
  corRfM <- corRf

  for (i in seq(nruns)) {
    cat(i, "\n")
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplitIndex(dat, nSplits = 2, respName = labelName)
    itrain <- splitDat$itrain
    itest <- splitDat$itest

    ## original
    cat("original", "\n")
    auxRf <- GetRfAUC(dat,
                      idxTrain = itrain,
                      idxTest = itest,
                      labelName,
                      featNames = featNames,
                      negClassName,
                      posClassName)
    aucs[i, "orig"] <- auxRf$aucObs
    R <- auxRf$predProbs
    Y <- dat[itest, labelName]
    Y <- as.numeric(Y)-1
    A <- dat[itest, confName]
    sdat <- data.frame(Y, A, R)
    auxTests <- CorTests(sdat, labelName = "Y", confName = "A", scoreName = "R")
    corRf[i, 1:5] <- auxTests$cors[, 1]
    corRf[i, 6:10] <- auxTests$pvals[, 1]

    ## matched
    cat("matching", "\n")
    splitDatM <- TrainTestRandomSplitIndex(datM, nSplits = 2, respName = labelName)
    itrainM <- splitDatM$itrain
    itestM <- splitDatM$itest
    auxRfM <- GetRfAUC(datM,
                       idxTrain = itrainM,
                       idxTest = itestM,
                       labelName,
                       featNames = featNames,
                       negClassName,
                       posClassName)
    aucs[i, "matching"] <- auxRfM$aucObs
    R <- auxRfM$predProbs
    Y <- datM[itestM, labelName]
    Y <- as.numeric(Y)-1
    A <- datM[itestM, confName]
    sdatM <- data.frame(Y, A, R)
    auxTestsM <- CorTests(sdatM, labelName = "Y", confName = "A", scoreName = "R")
    corRfM[i, 1:5] <- auxTestsM$cors[, 1]
    corRfM[i, 6:10] <- auxTestsM$pvals[, 1]
  }
  return(list(corRf = data.frame(corRf, check.names = FALSE), 
              corRfM = data.frame(corRfM, check.names = FALSE), 
              aucRf = data.frame(aucs, check.names = FALSE)))
}


#' Run Causality Tests based on Correlation (Ridge Regression)
#'
#' Runs the causality-based test for confounding
#' for the random forest classifiers
#'
#' @param dat unmatched data (dataframe, tibble)
#' @param datM matched data (dataframe, tibble)
#' @param nruns number of different runs (train/test splitting)
#' @param splitseeds random seed for data splitting
#' @param labelName name of prediction label
#' @param featNames name of desired feature names
#' @param confName name of confounding variable
#' @param negClassName boolean negative for diagnosis
#' @param posClassName boolean positive for diagnosis
#'
#' @examples
#' RunCondIndepTestsCorRr(
#' dat = tap.data,
#' datM = matched.tap.data,
#' nruns = 100, 123456,
#' labelName = "PD",
#' featNames = featNamesC,
#' confName = "age",
#' negClassName = FALSE,
#' posClassName = TRUE)
RunCondIndepTestsCorRr <- function(dat,
                                   datM,
                                   nruns,
                                   splitSeeds,
                                   labelName,
                                   featNames,
                                   confName,
                                   negClassName,
                                   posClassName) {
  GetRrAuc <- function(dat,
                       idxTrain,
                       idxTest,
                       respName,
                       featNames,
                       negClassName,
                       posClassName) {
    dat[, respName] <- factor(as.character(dat[, respName]),
                              levels = c(negClassName, posClassName))
    myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
    xtrain <- scale(dat[idxTrain, featNames])
    xtest <- scale(dat[idxTest, featNames])
    ytrain <- dat[idxTrain, respName]
    ytest <- dat[idxTest, respName]
    fit <- glmnet::cv.glmnet(x = xtrain, y = ytrain, family = "binomial",
                     alpha = 0)
    predProbs <- predict(fit, newx = xtest, type = "response", s = "lambda.min")
    rocObj <- pROC::roc(ytest, predProbs[, 1], direction = "<",
                  levels = c(negClassName, posClassName))
    aucObs <- pROC::auc(rocObj)[1]
    list(aucObs = aucObs, predProbs = predProbs[, 1], rocObj = rocObj)
  }
  TrainTestRandomSplitIndex <- function(dat,
                                        nSplits = 2,
                                        respName) {
    aux <- levels(dat[, respName])
    negativeClassName <- aux[1]
    positiveClassName <- aux[2]
    neg <- which(dat[, respName] == negativeClassName)
    pos <- which(dat[, respName] == positiveClassName)
    nNeg <- length(neg)
    nPos <- length(pos)
    testNeg <- sample(neg, round(nNeg/nSplits), replace = FALSE)
    trainNeg <- setdiff(neg, testNeg)
    testPos <- sample(pos, round(nPos/nSplits), replace = FALSE)
    trainPos <- setdiff(pos, testPos)
    train <- c(trainNeg, trainPos)
    test <- c(testNeg, testPos)
    list(itrain = train, itest = test)
  }

  aucs <- matrix(NA, nruns, 2)
  colnames(aucs) <- c("orig", "matching")

  corRr <- matrix(NA, nruns, 10)
  colnames(corRr) <- c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)",
                       "pval(R,Y)", "pval(R,A)", "pval(A,Y)", "pval(R,Y|A)", "pval(R,A|Y)")
  corRrM <- corRr

  for (i in seq(nruns)) {
    cat(i, "\n")
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplitIndex(dat, nSplits = 2, respName = labelName)
    itrain <- splitDat$itrain
    itest <- splitDat$itest

    ## original
    cat("original", "\n")
    auxRr <- GetRrAuc(dat,
                      idxTrain = itrain,
                      idxTest = itest,
                      labelName,
                      featNames,
                      negClassName,
                      posClassName)
    aucs[i, "orig"] <- auxRr$aucObs
    R <- auxRr$predProbs
    Y <- dat[itest, labelName]
    Y <- as.numeric(Y)-1
    A <- dat[itest, confName]
    sdat <- data.frame(Y, A, R)
    auxTests <- CorTests(sdat, labelName = "Y", confName = "A", scoreName = "R")
    corRr[i, 1:5] <- auxTests$cors[, 1]
    corRr[i, 6:10] <- auxTests$pvals[, 1]

    ## matched
    cat("matching", "\n")
    splitDatM <- TrainTestRandomSplitIndex(datM, nSplits = 2, respName = labelName)
    itrainM <- splitDatM$itrain
    itestM <- splitDatM$itest
    auxRrM <- GetRrAuc(datM,
                       idxTrain = itrainM,
                       idxTest = itestM,
                       labelName,
                       featNames,
                       negClassName,
                       posClassName)
    aucs[i, "matching"] <- auxRrM$aucObs
    R <- auxRrM$predProbs
    Y <- datM[itestM, labelName]
    Y <- as.numeric(Y)-1
    A <- datM[itestM, confName]
    sdatM <- data.frame(Y, A, R)
    auxTestsM <- CorTests(sdatM, labelName = "Y", confName = "A", scoreName = "R")
    corRrM[i, 1:5] <- auxTestsM$cors[, 1]
    corRrM[i, 6:10] <- auxTestsM$pvals[, 1]
  }
  return(list(corRr = data.frame(corRr, check.names = FALSE), 
              corRrM = data.frame(corRrM, check.names = FALSE), 
              aucRr = data.frame(aucs, check.names = FALSE)))
}


#' Runs distance correlation and partial distance correlation tests.
dCorTests <- function(nperm = 1000, dat, labelName, confName, scoreName) {
  dcors <- matrix(NA, 5, 1)
  colnames(dcors) <- "estimate"
  rownames(dcors) <- c(paste("dcor(", scoreName, " , ", labelName, ")", sep = ""),
                       paste("dcor(", scoreName, " , ", confName, ")", sep = ""),
                       paste("dcor(", labelName, " , ", confName, ")", sep = ""),
                       paste("dcor(", paste(scoreName, labelName, sep = " , "), " | ", confName, ")", sep = ""),
                       paste("dcor(", paste(scoreName, confName, sep = " , "), " | ", labelName, ")", sep = ""))
  pvals <- dcors
  colnames(pvals) <- "pval"
  permNulls <- matrix(NA, nperm, 5)
  colnames(permNulls) <- rownames(dcors)

  cat("running perm test 1", "\n")
  aux1 <- energy::dcor.test(dat[, scoreName], dat[, labelName], R = nperm)
  cat("running perm test 2", "\n")
  aux2 <- energy::dcor.test(dat[, scoreName], dat[, confName], R = nperm)
  cat("running perm test 3", "\n")
  aux3 <- energy::dcor.test(dat[, labelName], dat[, confName], R = nperm)
  cat("running perm test 4", "\n")
  aux4 <- energy::pdcor.test(dat[, scoreName], dat[, labelName], dat[, confName], R = nperm)
  cat("running perm test 5", "\n")
  aux5 <- energy::pdcor.test(dat[, scoreName], dat[, confName], dat[, labelName], R = nperm)

  dcors[1, 1] <- aux1$statistic
  dcors[2, 1] <- aux2$statistic
  dcors[3, 1] <- aux3$statistic
  dcors[4, 1] <- aux4$statistic
  dcors[5, 1] <- aux5$statistic

  pvals[1, 1] <- aux1$p.value
  pvals[2, 1] <- aux2$p.value
  pvals[3, 1] <- aux3$p.value
  pvals[4, 1] <- aux4$p.value
  pvals[5, 1] <- aux5$p.value

  permNulls[, 1] <- aux1$replicates
  permNulls[, 2] <- aux2$replicates
  permNulls[, 3] <- aux3$replicates
  permNulls[, 4] <- aux4$replicates
  permNulls[, 5] <- aux5$replicates
  return(list(dcors = dcors, pvals = pvals, permNulls = permNulls))
}


#' Run Causality Tests based on Distance Correlation (Random Forest)
#'
#' Runs the distance correaltion causality-based test for confounding for the random forest classifiers
#'
#' @param dat unmatched data (dataframe, tibble)
#' @param datM matched data (dataframe, tibble)
#' @param nruns number of different runs (train/test splitting)
#' @param splitseeds random seed for data splitting
#' @param labelName name of prediction label
#' @param featNames name of desired feature names
#' @param confName name of confounding variable
#' @param negClassName boolean negative for diagnosis
#' @param posClassName boolean positive for diagnosis
#'
#' @examples
#' RunCondIndepTestsdCorRf(
#' dat = tap.data,
#' datM = matched.tap.data,
#' nruns = 100, 123456,
#' labelName = "PD",
#' featNames = featNamesC,
#' confName = "age",
#' negClassName = FALSE,
#' posClassName = TRUE)
RunCondIndepTestsdCorRf <- function(dat,
                                    nruns,
                                    splitSeeds,
                                    labelName,
                                    featNames,
                                    confName,
                                    negClassName,
                                    posClassName,
                                    nperm = 1000) {
  GetRfAuc <- function(dat,
                       idxTrain,
                       idxTest,
                       labelName,
                       featNames,
                       negClassName,
                       posClassName) {
    dat <- dat[, c(labelName, featNames)]
    dat[, labelName] <- factor(as.character(dat[, labelName]),
                               levels = c(negClassName, posClassName))
    myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
    fit <- randomForest::randomForest(myFormula, data = dat[idxTrain,], ntree = 1000)
    predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "prob")
    rocObj <- pROC::roc(dat[idxTest, 1], predProbs[, posClassName], direction = "<",
                  levels = c(negClassName, posClassName))
    aucObs <- pROC::auc(rocObj)[1]
    list(aucObs = aucObs, predProbs = predProbs[, posClassName], rocObj = rocObj)
  }
  TrainTestRandomSplitIndex <- function(dat,
                                        nSplits = 2,
                                        respName) {
    aux <- levels(dat[, respName])
    negativeClassName <- aux[1]
    positiveClassName <- aux[2]
    neg <- which(dat[, respName] == negativeClassName)
    pos <- which(dat[, respName] == positiveClassName)
    nNeg <- length(neg)
    nPos <- length(pos)
    testNeg <- sample(neg, round(nNeg/nSplits), replace = FALSE)
    trainNeg <- setdiff(neg, testNeg)
    testPos <- sample(pos, round(nPos/nSplits), replace = FALSE)
    trainPos <- setdiff(pos, testPos)
    train <- c(trainNeg, trainPos)
    test <- c(testNeg, testPos)
    list(itrain = train, itest = test)
  }
  aucs <- matrix(NA, nruns, 2)
  colnames(aucs) <- c("orig", "matching")

  corRf <- matrix(NA, nruns, 10)
  colnames(corRf) <- c("cor(R,Y)", "cor(R,A)", "cor(A,Y)",
                       "cor(R,Y|A)", "cor(R,A|Y)",
                       "pval(R,Y)", "pval(R,A)", "pval(A,Y)",
                       "pval(R,Y|A)", "pval(R,A|Y)")
  corRfM <- corRf

  for (i in seq(nruns)) {
    cat(i, "\n")
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplitIndex(dat, nSplits = 2, respName = labelName)
    itrain <- splitDat$itrain
    itest <- splitDat$itest

    ## original
    cat("original", "\n")
    auxRf <- GetRfAuc(dat,
                      idxTrain = itrain,
                      idxTest = itest,
                      labelName,
                      featNames,
                      negClassName,
                      posClassName)
    aucs[i, "orig"] <- auxRf$aucObs
    R <- auxRf$predProbs
    Y <- dat[itest, labelName]
    Y <- as.numeric(Y)-1
    A <- dat[itest, confName]
    sdat <- data.frame(Y, A, R)
    auxTests <- dCorTests(nperm, sdat, labelName = "Y", confName = "A", scoreName = "R")
    corRf[i, 1:5] <- auxTests$dcors[, 1]
    corRf[i, 6:10] <- auxTests$pvals[, 1]
  }
  return(list(corRf = data.frame(corRf, check.names = FALSE)))
}


#' Run Causality Tests based on Distance Correlation (Ridge Regression)
#'
#' Runs the distance correaltion causality-based test for confounding for the random forest classifiers
#'
#' @param dat unmatched data (dataframe, tibble)
#' @param datM matched data (dataframe, tibble)
#' @param nruns number of different runs (train/test splitting)
#' @param splitseeds random seed for data splitting
#' @param labelName name of prediction label
#' @param featNames name of desired feature names
#' @param confName name of confounding variable
#' @param negClassName boolean negative for diagnosis
#' @param posClassName boolean positive for diagnosis
#'
#' @examples
#' RunCondIndepTestsdCorRr(
#' dat = tap.data,
#' datM = matched.tap.data,
#' nruns = 100, 123456,
#' labelName = "PD",
#' featNames = featNamesC,
#' confName = "age",
#' negClassName = FALSE,
#' posClassName = TRUE)
RunCondIndepTestsdCorRr <- function(dat,
                                    nruns,
                                    splitSeeds,
                                    labelName,
                                    featNames,
                                    confName,
                                    negClassName,
                                    posClassName,
                                    nperm = 1000) {
  GetRrAuc <- function(dat,
                       idxTrain,
                       idxTest,
                       respName,
                       featNames,
                       negClassName,
                       posClassName) {
    dat[, respName] <- factor(as.character(dat[, respName]),
                              levels = c(negClassName, posClassName))
    myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
    xtrain <- scale(dat[idxTrain, featNames])
    xtest <- scale(dat[idxTest, featNames])
    ytrain <- dat[idxTrain, respName]
    ytest <- dat[idxTest, respName]
    fit <- glmnet::cv.glmnet(x = xtrain, y = ytrain, family = "binomial",
                     alpha = 0)
    predProbs <- predict(fit, newx = xtest, type = "response", s = "lambda.min")
    rocObj <- pROC::roc(ytest, predProbs[, 1], direction = "<",
                  levels = c(negClassName, posClassName))
    aucObs <- pROC::auc(rocObj)[1]
    list(aucObs = aucObs, predProbs = predProbs[, 1], rocObj = rocObj)
  }
  TrainTestRandomSplitIndex <- function(dat,
                                        nSplits = 2,
                                        respName) {
    aux <- levels(dat[, respName])
    negativeClassName <- aux[1]
    positiveClassName <- aux[2]
    neg <- which(dat[, respName] == negativeClassName)
    pos <- which(dat[, respName] == positiveClassName)
    nNeg <- length(neg)
    nPos <- length(pos)
    testNeg <- sample(neg, round(nNeg/nSplits), replace = FALSE)
    trainNeg <- setdiff(neg, testNeg)
    testPos <- sample(pos, round(nPos/nSplits), replace = FALSE)
    trainPos <- setdiff(pos, testPos)
    train <- c(trainNeg, trainPos)
    test <- c(testNeg, testPos)
    list(itrain = train, itest = test)
  }

  aucs <- matrix(NA, nruns, 2)
  colnames(aucs) <- c("orig", "matching")

  corRr <- matrix(NA, nruns, 10)
  colnames(corRr) <- c("cor(R,Y)", "cor(R,A)", "cor(A,Y)",
                       "cor(R,Y|A)", "cor(R,A|Y)",
                       "pval(R,Y)", "pval(R,A)", "pval(A,Y)",
                       "pval(R,Y|A)", "pval(R,A|Y)")
  corRrM <- corRr

  for (i in seq(nruns)) {
    cat(i, "\n")
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplitIndex(dat, nSplits = 2, respName = labelName)
    itrain <- splitDat$itrain
    itest <- splitDat$itest

    ## original
    cat("original", "\n")
    auxRr <- GetRrAuc(dat,
                      idxTrain = itrain,
                      idxTest = itest,
                      labelName,
                      featNames,
                      negClassName,
                      posClassName)
    aucs[i, "orig"] <- auxRr$aucObs
    R <- auxRr$predProbs
    Y <- dat[itest, labelName]
    Y <- as.numeric(Y)-1
    A <- dat[itest, confName]
    sdat <- data.frame(Y, A, R)
    auxTests <- dCorTests(nperm, sdat, labelName = "Y", confName = "A", scoreName = "R")
    corRr[i, 1:5] <- auxTests$dcors[, 1]
    corRr[i, 6:10] <- auxTests$pvals[, 1]
  }
  return(list(corRr = data.frame(corRr, check.names = FALSE)))
}


#' Fit random forest and compute AUROC and balanced accuracy for the repeated measurement analyses.
GetAucBaccRf <- function(dat,
                         idxTrain,
                         idxTest,
                         labelName,
                         featNames,
                         negClassName,
                         posClassName) {
  ComputeBalancedAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    tp <- sum(ytest == 1 & yhat == 1)
    pos <- sum(ytest == 1)
    tn <- sum(ytest == 0 & yhat == 0)
    neg <- sum(ytest == 0)
    bacc <- (tp/pos + tn/neg)/2
    bacc
  }
  dat <- dat[, c(labelName, featNames)]
  dat[, labelName] <- factor(as.character(dat[, labelName]),
                             levels = c(negClassName, posClassName))
  myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
  fit <- randomForest::randomForest(myFormula, data = dat[idxTrain,])
  predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "prob")
  rocObj <- pROC::roc(dat[idxTest, 1], predProbs[, posClassName], direction = "<",
                levels = c(negClassName, posClassName))
  aucObs <- pROC::auc(rocObj)[1]
  bacc <- ComputeBalancedAccuracy(dat[idxTest, 1], predProbs[, posClassName], thr = 0.5)
  return(list(aucObs = aucObs, 
              predProbs = predProbs, 
              rocObj = rocObj, 
              bacc = bacc))
}


#' Fit ridge-regression and compute AUROC and balanced accuracy for the repeated measurement analyses.
GetAucBaccRr <- function(dat,
                         idxTrain,
                         idxTest,
                         labelName,
                         featNames,
                         negClassName,
                         posClassName) {
  ComputeBalancedAccuracy <- function(ytest, predProbs, thr = 0.5) {
    yhat <- ifelse(predProbs >= thr, 1, 0)
    ytest <- as.numeric(ytest) - 1
    tp <- sum(ytest == 1 & yhat == 1)
    pos <- sum(ytest == 1)
    tn <- sum(ytest == 0 & yhat == 0)
    neg <- sum(ytest == 0)
    bacc <- (tp/pos + tn/neg)/2
    bacc
  }
  dat <- dat[, c(labelName, featNames)]
  dat[, labelName] <- factor(as.character(dat[, labelName]),
                             levels = c(negClassName, posClassName))
  trainDat <- dat[idxTrain,]
  testDat <- dat[idxTest,]
  xtrain <- scale(trainDat[, featNames])
  xtest <- scale(testDat[, featNames])
  ytrain <- trainDat[, labelName]
  ytest <- testDat[, labelName]
  ## ridge regr fit
  fit <- glmnet::cv.glmnet(x = xtrain, y = ytrain, family = "binomial",
                   alpha = 0)
  predProbs <- predict(fit, newx = xtest, type = "response", s = "lambda.min")
  rocObj <- pROC::roc(ytest, predProbs[, 1], direction = "<",
                levels = c(negClassName, posClassName))
  aucObs <- pROC::auc(rocObj)[1]
  bacc <- ComputeBalancedAccuracy(dat[idxTest, 1], predProbs[, 1], thr = 0.5)
  return(list(aucObs = aucObs, 
              predProbs = predProbs, 
              rocObj = rocObj, 
              bacc = bacc))
}


#' Subjectwise train/test data splits.
#'
#' Run train test split given dataset with healthcode and
#' desired labels for repeated measurements
#'
#' @param dat featurized dataset
#' @param nSplits number of train/test splitting
#' @param subjectIdName subject index
#' @param labelName column for target variable
#' @param negClassName column name for negative label
#' @param posClassName column name for positve label
#'
#' @examples
#' GetIdxTrainTestSplitBySubject(dat = dat,
#' nSplits = 2,
#' subjectIdName = "healthCode",
#' labelName = "PD",
#' negClassName = "FALSE",
#' posClassName = "TRUE")
GetIdxTrainTestSplitBySubject <- function(dat,
                                          nSplits = 2,
                                          subjectIdName,
                                          labelName,
                                          negClassName = "FALSE",
                                          posClassName = "TRUE") {
  ids <- as.character(unique(dat[, subjectIdName]))
  labels <- dat[match(ids, dat[, subjectIdName]), labelName]
  caseIds <- ids[which(labels == posClassName)]
  controlIds <- ids[which(labels == negClassName)]
  nCase <- length(caseIds)
  nControl <- length(controlIds)
  testControlIds <- sample(controlIds, round(nControl/nSplits), replace = FALSE)
  trainControlIds <- setdiff(controlIds, testControlIds)
  testCaseIds <- sample(caseIds, round(nCase/nSplits), replace = FALSE)
  trainCaseIds <- setdiff(caseIds, testCaseIds)
  trainIds <- c(trainControlIds, trainCaseIds)
  testIds <- c(testControlIds, testCaseIds)
  idxTrain <- which(dat[, subjectIdName] %in% trainIds)
  idxTest <- which(dat[, subjectIdName] %in% testIds)
  return(list(idxTrain = idxTrain, idxTest = idxTest))
}

#' Subjectwise Label Shuffling
#'
#' Performs subjectwise shufflying of the data (for the shuffled
#' labels analyses)
#'
#' @param dat featurized dataframe (dataframe, tibble)
#' @param subjectIdName index of subject
#' @param labelName column of target variable
#'
#' @examples
#' SubjectWiseLabelShuffling0(tap.data, "healthCode", "PD")
SubjectWiseLabelShuffling0 <- function(dat, subjectIdName, labelName) {
  ids <- as.character(unique(dat[, subjectIdName]))
  nids <- length(ids)
  dat[, labelName] <- as.character(dat[, labelName])
  labels <- dat[match(ids, dat[, subjectIdName]), labelName]
  slabels <- labels[sample(nids)]
  for (i in seq(nids)) {
    dat[which(dat[, subjectIdName] == ids[i]), labelName] <- slabels[i]
  }
  dat[, labelName] <- as.factor(dat[, labelName])
  return(dat)
}



###############################################
## Functions for the variability comparisons
## PD versus non-PD participants
###############################################


#' compute IQR values given feature
GetIQR <- function(dat, featNames) {
  participants <- unique(dat$healthCode)
  npar <- length(participants)
  nfeat <- length(featNames)
  iqrs <- matrix(NA, npar, nfeat)
  rownames(iqrs) <- participants
  colnames(iqrs) <- featNames
  for (i in seq(npar)) {
    sdat <- dat[which(dat$healthCode == participants[i]), featNames]
    iqrs[i,] <- apply(sdat, 2, IQR, na.rm = TRUE)
  }
  return(iqrs)
}


#' Detrend data given features
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

#' Run Variabiltiy Comparison Permutation Test
#'
#' Run the permutation p-values for the variability comparison.
#' for each activity and given list of features
#'
#' @param nperm number of permutation
#' @param dat featurized data
#' @param selectedHelathCodes choice of subsampled healthcodes
#' @param featNames list of features
#'
#' @examples
#' RunVariabilityComparisonPermTest(nperm = 1000,
#' dat = tap.data,
#' selectedHealthCodes = <choice of healthCode>,
#' featNames = tapFeatures)
RunVariabilityComparisonPermTest <- function(nperm, dat, selectedHealthCodes, featNames) {
  set.seed(1234567)
  dat <- dat[dat$healthCode %in% selectedHealthCodes,]
  controls <- dat[which(dat$PD == "FALSE"),]
  cases <- dat[which(dat$PD == "TRUE"),]

  cat("detrending controls", "\n")
  controlsD <- LoessDetrendedFeatures(dat = controls, featNames = featNames)
  cat("detrending cases", "\n")
  casesD <- LoessDetrendedFeatures(dat = cases, featNames = featNames)

  cat("observed median IQR difference", "\n")
  casesV <- GetIQR(dat = casesD, featNames = featNames)
  controlsV <- GetIQR(dat = controlsD, featNames = featNames)
  obs <- apply(casesV, 2, median) - apply(controlsV, 2, median)

  cat("run permutations", "\n")
  nfeat <- length(featNames)
  permM <- matrix(NA, nperm, nfeat)
  colnames(permM) <- featNames
  datD <- rbind(casesD, controlsD)
  pdatD <- datD
  n <- nrow(datD)
  for (i in seq(nperm)) {
    cat("perm ", i, "\n")
    pdatD[, featNames] <- datD[sample(n, replace = FALSE), featNames]
    pcontrolsD <- pdatD[which(pdatD$PD == "FALSE"),]
    pcasesD <- pdatD[which(pdatD$PD == "TRUE"),]
    pcasesV <- GetIQR(dat = pcasesD, featNames = featNames)
    pcontrolsV <- GetIQR(dat = pcontrolsD, featNames = featNames)
    permM[i,] <- apply(pcasesV, 2, median) - apply(pcontrolsV, 2, median)
  }
  return(list(obs = obs, permM = permM))
}

#' Get Variability in P-Values(?)
GetVariabilityPval <- function(x) {
  obs <- x$obs
  nfeat <- length(obs)
  pvals <- matrix(NA, nfeat, 3)
  rownames(pvals) <- names(obs)
  colnames(pvals) <- c("pval", "obs", "approx")
  nperm <- nrow(x$permM)
  for (i in seq(nfeat)) {
    p1 <- (sum(x$permM[, i] > obs[i]) + 1)/(nperm + 1)
    p2 <- (sum(x$permM[, i] < obs[i]) + 1)/(nperm + 1)
    pvals[i, 2] <- obs[i]
    pvals[i, 1] <- 2 * min(c(p1, p2))
    psd <- sd(x$permM[, i])
    pm <- mean(x$permM[, i])
    pvals[i, 3] <- abs(obs[i] - pm)/psd
  }
  return(pvals)
}


#' Filter Close Interval Activities with a certain threshold
FilterOutActivitiesCloseInTime <- function(dat, secondsThr = 3600) {
  dat$createdOn <- as.POSIXct(dat$createdOn)
  hc <- as.character(unique(dat$healthCode))
  npar <- length(hc)

  sdat <- dat[dat$healthCode == hc[1],]
  ## sort by time
  sdat <- sdat[order(sdat$createdOn, decreasing = FALSE),]
  times <- sdat$createdOn
  ## get delta time in seconds
  deltas <- diff(times, lag = 1)
  ## keep indexes of activities that are at least secondsThr apart
  tokeep <- c(1, which(deltas >= secondsThr) + 1)
  sdat <- sdat[tokeep,]
  fdat <- sdat

  for (i in 2:npar) {
    sdat <- dat[dat$healthCode == hc[i],]
    sdat <- sdat[order(sdat$createdOn, decreasing = FALSE),]
    times <- sdat$createdOn
    deltas <- diff(times, lag = 1)
    tokeep <- c(1, which(deltas >= secondsThr) + 1)
    sdat <- sdat[tokeep,]
    fdat <- rbind(fdat, sdat)
  }
  return(fdat)
}


#' Get IQR distributions for PD and non-PD cases (produces the output for IQR boxplots figure)
#'
#' @param dat featurized data (tibble, dataframe)
#' @param selectedHealthCodes subsample of healthcodes being used
#' @param featNames list of features
VariabilityComparison <- function(dat, selectedHealthCodes, featNames) {
  set.seed(1234567)
  GetIQR <- function(dat, featNames) {
    participants <- unique(dat$healthCode)
    npar <- length(participants)
    nfeat <- length(featNames)
    iqrs <- matrix(NA, npar, nfeat)
    rownames(iqrs) <- participants
    colnames(iqrs) <- featNames
    for (i in seq(npar)) {
      sdat <- dat[which(dat$healthCode == participants[i]), featNames]
      iqrs[i,] <- apply(sdat, 2, IQR, na.rm = TRUE)
    }
    iqrs
  }
  dat <- dat[dat$healthCode %in% selectedHealthCodes,]
  controls <- dat[which(dat$PD == "FALSE"),]
  cases <- dat[which(dat$PD == "TRUE"),]
  controlsD <- LoessDetrendedFeatures(dat = controls, featNames = featNames)
  casesD <- LoessDetrendedFeatures(dat = cases, featNames = featNames)
  casesV <- GetIQR(dat = casesD, featNames = featNames)
  controlsV <- GetIQR(dat = controlsD, featNames = featNames)
  return(list(casesV = casesV, controlsV = controlsV))
}

## Performs exact age matching for males and females separately,
## and breaks ties by selecting the participants with the 
## largest number of records. (Further ties are broken by 
## randomly selecting tied participants.)
SeparateMaleFemaleAgeMatching <- function(dat) {
  GetPairs <- function(mdat) {
    aux <- table(mdat$subclass)
    subclasses <- names(aux)
    nclasses <- length(subclasses)
    cases <- c()
    controls <- c()
    ageCases <- c()
    ageControls <- c()
    for (i in seq(nclasses)) {
      idx0 <- which(mdat$subclass == subclasses[i] & mdat$pd == 0)
      idx1 <- which(mdat$subclass == subclasses[i] & mdat$pd == 1)
      mdat0 <- mdat[idx0,]
      mdat1 <- mdat[idx1,]
      n0 <- length(idx0) ## number of controls
      n1 <- length(idx1) ## number of cases
      
      ## when we have more controls (n0) than cases (n1),
      ## we select all cases, then sort the controls by
      ## number of records (nrecs) and select the top n1
      ## controls
      if (n0 > n1) {
        cases <- c(cases, mdat1$healthCode)
        ageCases <- c(ageCases, mdat1$age)
        mdat0 <- mdat0[order(mdat0$nrecs, decreasing = TRUE),]
        controls <- c(controls, mdat0$healthCode[seq(n1)])
        ageControls <- c(ageControls, mdat0$age[seq(n1)])
      }
      
      ## when we have more cases (n1) than controls (n0),
      ## we select all controls, then sort the cases by
      ## number of records (nrecs) and select the top n0
      ## cases
      if (n1 > n0) {
        controls <- c(controls, mdat0$healthCode)
        ageControls <- c(ageControls, mdat0$age)
        mdat0 <- mdat0[order(mdat0$nrecs, decreasing = TRUE),]
        cases <- c(cases, mdat1$healthCode[seq(n0)])
        ageCases <- c(ageCases, mdat1$age[seq(n0)])
      }
      
      ## when the number of cases is the same as controls
      ## we just select all cases and all controls
      if (n0 == n1) { 
        cases <- c(cases, mdat1$healthCode)
        ageCases <- c(ageCases, mdat1$age)
        controls <- c(controls, mdat0$healthCode)
        ageControls <- c(ageControls, mdat0$age)
      }
      
    }
    healthCode <- c(cases, controls)
    ncases <- length(cases)
    ncontrols <- length(controls)
    pd <- c(rep(1, ncases), rep(0, ncontrols))
    age <- as.numeric(c(ageCases, ageControls))
    idx <- c(seq(ncases), seq(ncontrols))
    data.frame(healthCode, pd, age, idx) 
  }
  dat$pd <- as.numeric(as.factor(dat$PD)) - 1
  dat <- dat[, c("healthCode", "pd", "gender", "age", "nrecs")]
  dat <- na.omit(dat)
  datM <- dat[dat$gender == "Male",]
  datF <- dat[dat$gender == "Female",]
  mM <- MatchIt::matchit(pd ~ age, data = datM, method = "exact")
  mdatM <- MatchIt::match.data(mM)
  hcM <- GetPairs(mdatM)
  hcM$gender <- rep("Male", nrow(hcM))
  mF <- MatchIt::matchit(pd ~ age, data = datF, method = "exact")
  mdatF <- MatchIt::match.data(mF)
  hcF <- GetPairs(mdatF)  
  hcF$gender <- rep("Female", nrow(hcF))
  return(rbind(hcM, hcF))
}


## get collapsed data 
##
GetCollapsedData <- function(x, labelName, covNames, subjectIdName, featNames) {
  ids <- as.character(unique(x[, subjectIdName]))
  nids <- length(ids)
  nfeat <- length(featNames)
  nvar <- length(c(subjectIdName, labelName, covNames))
  out <- data.frame(matrix(NA, nids, 2 * nfeat + nvar))
  cfeatNames <- c(paste(featNames, "med", sep = "."), paste(featNames, "iqr", sep = "."))
  colnames(out) <- c(subjectIdName, labelName, covNames, cfeatNames)
  rownames(out) <- ids
  for (i in seq(nids)) {
    sdat <- x[which(x[, subjectIdName] == ids[i]),]
    out[i, subjectIdName] <- ids[i]
    out[i, labelName] <- as.character(sdat[1, labelName])
    out[i, covNames] <- as.character(sdat[1, covNames])
    out[i, (nvar + 1):(nfeat + nvar)] <- apply(sdat[, featNames], 2, median, na.rm = TRUE)
    out[i, (nfeat + nvar + 1):(2 * nfeat + nvar)] <- apply(sdat[, featNames], 2, IQR, na.rm = TRUE)
  }
  for (j in seq(length(covNames))) {
    out[which(out[, covNames[j]] == "NA"), covNames[j]] <- NA
    out[, covNames[j]] <- out[, covNames[j]]
  }
  return(list(out = out, cfeatNames = cfeatNames))
}



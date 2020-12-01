## Install dependencies
##
library(install.load)
install_load("pROC", "randomForest", "plyr", "dplyr")
install_load("doMC")


####################################
## analysis functions
####################################

## Perform subject-wise label shuffling.
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
    
    dat
}


## This function makes sure that we have both labels in the training and test sets.
## It keeps generating new shuffles of the label data until both labels are present
## in the training and test sets.
## 
SubjectWiseLabelShuffling <- function(dat, subjectIdName, labelName, idxTrain, idxTest) {
    aux <- 0
    while (aux == 0) {
        datS <- SubjectWiseLabelShuffling0(dat, subjectIdName, labelName)
        labelsTrain <- unique(datS[idxTrain, labelName])
        labelsTest <- unique(datS[idxTest, labelName])
        if (length(labelsTrain) == 2 & length(labelsTest) == 2) {
            aux <- 1
        }
    }
    
    datS
}


## Generates the permutation null distributions
## using restricted permutations.
##
DRPermDistrAUC <- function(dat, 
                           idxTrain,
                           idxTest,
                           nperm, 
                           subjectIdName, 
                           labelName, 
                           featNames,
                           negClassName, 
                           posClassName,
                           verbose = FALSE,
                           parallel = TRUE) {
    dat <- dat[, c(labelName, featNames, subjectIdName)]
    dat[, labelName] <- factor(as.character(dat[, labelName]), 
                               levels = c(negClassName, posClassName)) 
    myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
    if (verbose) {
        res_auc <- plyr::llply(1:nperm, .parallel = parallel,  function(num){
            datPS <- SubjectWiseLabelShuffling(dat, subjectIdName, labelName, idxTrain, idxTest)
            fitPS <- randomForest(myFormula, data = datPS[idxTrain,])
            predProbsPS <- predict(fitPS, datPS[idxTest, -1, drop = FALSE], type = "prob")
            rocObjPS <- roc(datPS[idxTest, 1], predProbsPS[, posClassName], direction = "<", 
                            levels = c(negClassName, posClassName)) 
            pROC::auc(rocObjPS)[1]
        }, .progress = progress_text(char = "."))
    }
    else {
        res_auc <- plyr::llply(1:nperm, .parallel = parallel,  function(num){
            datPS <- SubjectWiseLabelShuffling(dat, subjectIdName, labelName, idxTrain, idxTest)
            fitPS <- randomForest(myFormula, data = datPS[idxTrain,])
            predProbsPS <- predict(fitPS, datPS[idxTest, -1, drop = FALSE], type = "prob")
            rocObjPS <- roc(datPS[idxTest, 1], predProbsPS[, posClassName], direction = "<", 
                            levels = c(negClassName, posClassName)) 
            pROC::auc(rocObjPS)[1]
        })
    }
    
    unlist(res_auc)
}


## Compute the AUC score, and pseudo p-value.
##
GetAUC <- function(dat,
                   idxTrain, 
                   idxTest, 
                   subjectIdName, 
                   labelName, 
                   featNames,
                   negClassName, 
                   posClassName) {
    GetNormApproxVarAUC <- function(ytest, predProbs, negClassName, posClassName) {
        GetTieStats <- function(x) {
            tj <- 0
            u <- unique(x)
            if (length(x) > length(u)) {
                idupli <- which(duplicated(x))
                ud <- unique(x[idupli])
                tau <- length(ud)
                tj <- rep(NA, tau)
                for (i in seq(tau)) {
                    tj[i] <- sum(x == ud[i])
                }
            }
            
            list(tj = tj, aux = sum(tj * (tj - 1) * (tj + 1)))
        }
        ytest <- factor(ytest)
        ylevels <- levels(ytest)
        n1 <- sum(ytest == negClassName)
        n2 <- sum(ytest == posClassName)
        ties <- GetTieStats(predProbs)
        n <- n1 + n2
        v <- (n + 1)/(12 * n1 * n2) - ties[[2]]/(12 * n1 * n2 * n * (n - 1))
        
        c(v = v, n = n, nNeg = n1, nPos = n2, statTies = ties[[2]])
    }
    dat <- dat[, c(labelName, featNames, subjectIdName)]
    dat[, labelName] <- factor(as.character(dat[, labelName]), 
                               levels = c(negClassName, posClassName)) 
    myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
    fit <- randomForest(myFormula, data = dat[idxTrain,])
    predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "prob")
    rocObj <- roc(dat[idxTest, 1], predProbs[, posClassName], direction = "<", 
                  levels = c(negClassName, posClassName))    
    aucObs <- pROC::auc(rocObj)[1]
    ## get variance to compute the pseudo p-value
    approxVar <- GetNormApproxVarAUC(dat[idxTest, 1], predProbs, negClassName, posClassName)
    
    list(aucObs = aucObs, rocObj = rocObj, approxVar = approxVar)
}



## Perform record-wise data split (where it is possible
## that part of the records of each subject is assigned to
## the training set, while the remaining records are assigned 
## to the test set).
##
GetIdxTrainTestSplitByRecord <- function(dat, nSplits = 2) {
    n <- nrow(dat)
    idxTrain <- sample(n, round(n/2), replace = FALSE)
    idxTest <- setdiff(seq(n), idxTrain)
    
    list(idxTrain = idxTrain, idxTest = idxTest)
}

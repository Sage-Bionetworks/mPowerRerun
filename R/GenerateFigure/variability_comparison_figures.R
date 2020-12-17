######################################################
# This script is used to run supplementary figure 5 and
# 6 of the mPower paper. It will create figure
# for feature variability comparison for each activity
#' @author_email: echaibub@synapse.org, 
#' aryton.tediarjo@sagebase.org
#######################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

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
SCRIPT_NAME <-  "variability_comparison_figures.R"
SYN_ID_REF <- list(figures = get_figure_ref(),
                   intermediate = get_intermediate_data_ref())
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
VPERM_PD_v_nonPD_SYN_ID <- SYN_ID_REF$intermediate$vperm_pd_vs_nonpd
FEATURE_LIST <- get_features()
SUPPL_FIGURE_5 <- paste0("mPower_",
                         gsub(" ", "_", get("metadata")$user_group), 
                         "_supplementary_figure_5",".png")
SUPPL_FIGURE_6 <- paste0("mPower_",
                         gsub(" ", "_", get("metadata")$user_group), 
                         "_supplementary_figure_6",".png")
ANNOTATIONS <-list(analysisType = "case vs controls",
                   analysisSubtype = "feature variability comparison",
                   userSubset = tolower(get("metadata")$user_group),
                   pipelineStep = "figures")

######################################################
## helper
#######################################################
## compute the permutation p-values
## (the third column of the output
## is used to rank the results)
GetPval <- function(x) {
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

## function to use use healthcode as row names
set_healthcode_as_index <- function(vComparison){
  map(vComparison, function(x){
    rownames(x) <- x$healthCode
    x <- x %>% 
      dplyr::select(-c(healthCode))
    return(x)
  })
}


main <- function(){
  
  #######################################################
  ## Instantiate Variables and Reference IDs
  #######################################################
  pd_vs_non_pd_variability_perm_test <- H5Fopen(synGet(VPERM_PD_v_nonPD_SYN_ID)$path)
  tapFeatures <- FEATURE_LIST$tap
  walkFeatures <- FEATURE_LIST$walk
  restFeatures <- FEATURE_LIST$rest
  voiceFeatures <- FEATURE_LIST$voice
  
  
  ##############################
  ## Supplementary Figure 5
  ##############################
  ptap <- h5read(pd_vs_non_pd_variability_perm_test, "/tap/permPvals")
  ptap$obs <- unlist(ptap$obs)
  pwal <- h5read(pd_vs_non_pd_variability_perm_test, "/walk/permPvals")
  pwal$obs <- unlist(pwal$obs)
  pres <- h5read(pd_vs_non_pd_variability_perm_test, "/rest/permPvals")
  pres$obs <- unlist(pres$obs)
  pvoi <- h5read(pd_vs_non_pd_variability_perm_test, "/voice/permPvals")
  pvoi$obs <- unlist(pvoi$obs)
  
  dtap <- GetPval(ptap)
  dwal <- GetPval(pwal)
  dres <- GetPval(pres)
  dvoi <- GetPval(pvoi)
  
  rownames(dtap) <- paste("tap", rownames(dtap), sep = "___")
  rownames(dwal) <- paste("wal", rownames(dwal), sep = "___")
  rownames(dres) <- paste("res", rownames(dres), sep = "___")
  rownames(dvoi) <- paste("voi", rownames(dvoi), sep = "___")
  
  dp <- data.frame(rbind(dtap, dwal, dres, dvoi))
  dp$apval <- p.adjust(dp$pval, method = "BH")
  odp <- dp[order(dp$approx, decreasing = TRUE),]
  odp <- odp[odp$apval <= 0.05,]
  blwd <- 2
  cm <- 0.95
  ca <- 1
  cl <- 1
  xlab <- "S statistic"
  ylab <- "Density"
  
  
  title <- SUPPL_FIGURE_5
  png(paste("", title, sep = ""), width = 2000, 
      height = 2000, res = 200)
  
  nms <- rownames(odp)
  par(mfrow = c(as.integer(nrow(odp) / 5) + 1, 5), 
      mar = c(4, 3.0, 1.2, 0.5), mgp = c(1.75, 0.5, 0))
  for (i in  seq(length(rownames(odp)))) {
    aux <- strsplit(nms[i], "___")[[1]]
    if (aux[1] == "tap") {
      idx <- match(aux[2], tapFeatures)
      hist(ptap$permM[, idx], 
           col = "white",
           probability = TRUE, 
           main = paste(tapFeatures[idx], "(tap)", sep = " "),
           xlim = c(min(c(ptap$obs[idx], min(ptap$permM[, idx]))), max(c(ptap$obs[idx], max(ptap$permM[,idx])))),
           cex.main = cm, cex.axis = ca, cex.lab = cl, xlab = xlab, ylab = ylab)
      abline(v = 0, lwd = blwd)
      abline(v = ptap$obs[idx], col = "red", lwd = blwd)
    }
    if (aux[1] == "wal") {
      idx <- match(aux[2], walkFeatures) 
      hist(pwal$permM[, idx], 
           col = "white",
           probability = TRUE, 
           main = paste(walkFeatures[idx], "(wal)", sep = " "),
           xlim = c(min(c(pwal$obs[idx], min(pwal$permM[, idx]))), max(c(pwal$obs[idx], max(pwal$permM[,idx])))),
           cex.main = cm, cex.axis = ca, cex.lab = cl, xlab = xlab, ylab = ylab)
      abline(v = 0, lwd = blwd)
      abline(v = pwal$obs[idx], col = "red", lwd = blwd)
    }  
    if (aux[1] == "res") {
      idx <- match(aux[2], restFeatures)
      hist(pres$permM[, idx], 
           col = "white",
           probability = TRUE, main = paste(restFeatures[idx], "(res)", sep = " "),
           xlim = c(min(c(pres$obs[idx], min(pres$permM[, idx]))), max(c(pres$obs[idx], max(pres$permM[,idx])))),
           cex.main = cm, cex.axis = ca, cex.lab = cl, xlab = xlab, ylab = ylab)
      abline(v = 0, lwd = blwd)
      abline(v = pres$obs[idx], col = "red", lwd = blwd)
    }    
    if (aux[1] == "voi") {
      idx <- match(aux[2], voiceFeatures)
      hist(pvoi$permM[, idx], 
           col = "white",
           probability = TRUE, main = paste(voiceFeatures[idx], "(voi)", sep = " "),
           xlim = c(min(c(pvoi$obs[idx], min(pvoi$permM[, idx]))), max(c(pvoi$obs[idx], max(pvoi$permM[,idx])))),
           cex.main = cm, cex.axis = ca, cex.lab = cl, xlab = xlab, ylab = ylab)
      abline(v = 0, lwd = blwd)
      abline(v = pvoi$obs[idx], col = "red", lwd = blwd)
    }     
  }
  dev.off()
  
  ##############################
  ## Store to Synapse
  ##############################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(executed = GIT_URL, 
                                  used = c(VPERM_PD_v_nonPD_SYN_ID)))
  unlink(title)
  
  
  ##############################
  ## Supplementary Figure 6
  ##############################
  vTap <- set_healthcode_as_index(
    h5read(pd_vs_non_pd_variability_perm_test, "/tap/varComp"))
  
  vWal <- set_healthcode_as_index(
    h5read(pd_vs_non_pd_variability_perm_test, "/walk/varComp"))
  
  vRes <- set_healthcode_as_index(
    h5read(pd_vs_non_pd_variability_perm_test, "/rest/varComp"))
  
  vVoi <- set_healthcode_as_index(
    h5read(pd_vs_non_pd_variability_perm_test, "/voice/varComp"))
  
  ## get the names of the features that showed 
  ## significant differences
  nms <- rownames(odp)
  aux <- unlist(strsplit(nms, "___"))
  task <- aux[seq(1, 2*length(nms), by = 2)]
  sigFeat <- aux[seq(2, 2*length(nms), by = 2)]
  
  colorPD <- "#DC0000FF"
  colorNonPD <- "#3C5488FF"
  
  cm <- 1.4
  ca <- 1.2
  cl <- 1.3
  cleg <- 0.75
  nfeat <- length(sigFeat)
  
  title <- SUPPL_FIGURE_6
  png(paste("", title, sep = ""), width = 2000, height = 2000, res = 200)
  
  par(mfrow = c(5, 5), mar = c(2.5, 1.5, 1, 0.5), mgp = c(2, 0.5, 0))
  for (i in seq(nfeat)) {
    if (task[i] == "tap") {
      boxplot(vTap$casesV[, sigFeat[i]], vTap$controlsV[, sigFeat[i]], outline = FALSE, 
              main = sigFeat[i], border = "black", 
              names = c("PD", "non-PD"),col = c(colorPD, colorNonPD), cex.main = cm, 
              cex.axis = ca, cex.lab = cl) 
    }
    if (task[i] == "res") {
      boxplot(vRes$casesV[, sigFeat[i]], vRes$controlsV[, sigFeat[i]], outline = FALSE, 
              main = sigFeat[i], border = "black", 
              names = c("PD", "non-PD"),col = c(colorPD, colorNonPD), cex.main = cm,
              cex.axis = ca, cex.lab = cl) 
    }
    if (task[i] == "wal") {
      boxplot(vWal$casesV[, sigFeat[i]], vWal$controlsV[, sigFeat[i]], outline = FALSE, 
              main = sigFeat[i], border = "black", 
              names = c("PD", "non-PD"),col = c(colorPD, colorNonPD), cex.main = cm,
              cex.axis = ca, cex.lab = cl)  
    }
    if (task[i] == "voi") {
      boxplot(vVoi$casesV[, sigFeat[i]], vVoi$controlsV[, sigFeat[i]], outline = FALSE, 
              main = sigFeat[i], border = "black", 
              names = c("PD", "non-PD"),col = c(colorPD, colorNonPD), cex.main = cm,
              cex.axis = ca, cex.lab = cl)
    }
  }
  dev.off()
  
  ##############################
  ## Store to Synapse
  ##############################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- list(study = tolower(get("metadata")$study),
                        analysisType = "case vs controls",
                        analysisSubtype = "feature variability comparison",
                        userSubset = tolower(get("metadata")$user_group),
                        pipelineStep = "figures")
  synStore(
    f, activity = Activity(
      "generate feature comparison",
      executed = GIT_URL, 
      used = c(VPERM_PD_v_nonPD_SYN_ID)))
  unlink(title)
  
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

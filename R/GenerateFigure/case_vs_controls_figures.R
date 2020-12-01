###################################################
#' This script is used to generate supplementary figure 
#' 3 and 4 of the mPower Paper
#' It will graph the AUC classification performance of 
#' each activity of PD cases vs controls
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
##################################################
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
SYN_ID_REF <- list(figures = get_figure_ref(),
                   intermediate = get_intermediate_data_ref())
SCRIPT_NAME <-  "case_vs_controls_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
UNREPEATED_PD_v_NONPD_SYN_ID <- SYN_ID_REF$intermediate$collapsed_pd_vs_nonpd
REPEATED_PD_v_NONPD_SYN_ID <- SYN_ID_REF$intermediate$repeated_pd_vs_nonpd
SUPPL_FIGURE_3 <- paste0("mPower_",
                         gsub(" ", "_", get("metadata")$user_group), 
                         "_supplementary_figure_3",".png")
SUPPL_FIGURE_4 <- paste0("mPower_",
                         gsub(" ", "_", get("metadata")$user_group), 
                         "_supplementary_figure_4",".png")
ANNOTATIONS_SUPPL_3 <- list(
  study = tolower(get("metadata")$study),
  analysisType = "case vs controls",
  analysisSubtype = "performance using random forests",
  userSubset = tolower(get("metadata")$user_group), 
  pipelineStep = "figures")
ANNOTATIONS_SUPPL_4 <- list(
  study = tolower(get("metadata")$study),
  analysisType = "case vs controls",
  analysisSubtype = "performance using ridge regression",
  userSubset = tolower(get("metadata")$user_group), 
  pipelineStep = "figures")

######################################################
## helpers
#######################################################

#' wrapper function to generate figure 3-4
get_supplementary_figures <- function(datS, datSP, datS2, datSP2,
                                      datC, datCP, datC2, datCP2){
        
        activity <- c("tapping", "walk", "balance", "voice")
        
        data.list <- map(list(datS = datS, 
                              datSP = datSP, 
                              datS2 = datS2, 
                              datSP2 = datSP2,
                              datC = datC, 
                              datCP = datCP, 
                              datC2 = datC2, 
                              datCP2 = datCP2), 
                         function(x){
                                 names(x) <- activity
                                 return(x)
                            })
        colorTap <- "#00A087FF"
        colorVoi <- "#7E6148FF"
        colorWal <- "#E64B36FF"
        colorRes <- "#8491B4FF"
        par(mfrow = c(2, 2), mar = c(4, 3, 3, 1), mgp = c(2, 0.75, 0))
        boxplot(data.list$datCP, border = rep("grey", 4), ylim = c(0.3, 0.9), at = seq(4), 
                ylab = "AUC", col = "white",
                main = "collapsed features (median and IQR)")
        boxplot(data.list$datC, border = c(colorTap, colorWal, colorRes, colorVoi), 
                ylim = c(0.3, 0.9), at = seq(4), ylab = "AUC",  col = "white",
                main = "collapsed features (median and IQR)", add = TRUE)
        abline(h = 0.5)
        mtext("(a)", side = 3, at = 4.3, cex = 1.5, line = -2)
       
         boxplot(data.list$datSP, border = rep("grey", 4), ylim = c(0.3, 0.9), 
                at = seq(4), ylab = "AUC", col = "white",
                main = "repeated measurements (subject-wise split)")
        boxplot(data.list$datS, border = c(colorTap, colorWal, colorRes, colorVoi), ylim = c(0.3, 0.9), 
                at = seq(4), ylab = "AUC", col = "white",
                main = "repeated measurements (subject-wise split)", add = TRUE)
        abline(h = 0.5)
        mtext("(b)", side = 3, at = 4.3, cex = 1.5, line = -2)
        ####
        
        boxplot(data.list$datCP2, border = rep("grey", 4), ylim = c(0.3, 0.9), 
                at = seq(4), ylab = "Accuracy", col = "white",
                main = "collapsed features (median and IQR)")
        boxplot(data.list$datC2, border = c(colorTap, colorWal, colorRes, colorVoi), 
                ylim = c(0.3, 0.9), at = seq(4), ylab = "Accuracy", col = "white",
                main = "collapsed features (median and IQR)", add = TRUE)
        abline(h = 0.5)
        mtext("(c)", side = 3, at = 4.3, cex = 1.5, line = -2)
        
        boxplot(data.list$datSP2, border = rep("grey", 4), ylim = c(0.3, 0.9), 
                at = seq(4), ylab = "Balanced accuracy", col = "white",
                main = "repeated measurements (subject-wise split)")
        boxplot(data.list$datS2, border = c(colorTap, colorWal, colorRes, colorVoi), ylim = c(0.3, 0.9), 
                at = seq(4), ylab = "Balanced accuracy", col = "white",
                main = "repeated measurements (subject-wise split)", add = TRUE)
        abline(h = 0.5)
        mtext("(d)", side = 3, at = 4.3, cex = 1.5, line = -2)
}



main <- function(){
  
  # retrieve model intermediary weights scorings
  unrepeated_pd_vs_non_pd <- H5Fopen(synGet(
    UNREPEATED_PD_v_NONPD_SYN_ID)$path)
  
  # retrieve model intermediary weights scorings
  repeated_pd_vs_non_pd <- H5Fopen(synGet(
    REPEATED_PD_v_NONPD_SYN_ID)$path)
  
  
  ###################################
  ## Supplementary Figure 3
  ###################################
  metric <- "Auc"
  datC <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                            "/tap/matched/rf2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/walk/matched/rf2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/rest/matched/rf2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/voice/matched/rf2")$metrics[[metric]])
  datCP <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                             "/tap/matched/rf2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/walk/matched/rf2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/rest/matched/rf2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/voice/matched/rf2")$metricsP[[metric]])
  
  metric <- "BalAcc"
  datC2 <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                             "/tap/matched/rf2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/walk/matched/rf2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/rest/matched/rf2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/voice/matched/rf2")$metrics[[metric]])
  datCP2 <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                              "/tap/matched/rf2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/walk/matched/rf2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/rest/matched/rf2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/voice/matched/rf2")$metricsP[[metric]])
  
  #### Retrieve AUC vs Balanced Accuracy ###
  metric <- "auc"
  datS <- data.frame(h5read(repeated_pd_vs_non_pd, 
                            "tap/rf_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "walk/rf_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "rest/rf_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "voice/rf_unshuffled")[[metric]])

  datSP <- data.frame(h5read(repeated_pd_vs_non_pd, 
                             "tap/rf_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "walk/rf_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "rest/rf_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "voice/rf_shuffled")[[metric]])
  metric <- "bacc"
  datS2 <- data.frame(h5read(repeated_pd_vs_non_pd, 
                             "tap/rf_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "walk/rf_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "rest/rf_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "voice/rf_unshuffled")[[metric]])
  
  datSP2 <- data.frame(h5read(repeated_pd_vs_non_pd, 
                              "tap/rf_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "walk/rf_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "rest/rf_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "voice/rf_shuffled")[[metric]])
  
  ### save figure ###
  title <- SUPPL_FIGURE_3
  png(paste(title), width = 2000, height = 2000, res = 200)
  get_supplementary_figures(datS, datSP, datS2, datSP2,
                            datC, datCP, datC2, datCP2)
  dev.off()
  
  #############################
  ## Store to Synapse
  #############################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS_SUPPL_3
  synStore(
    f, activity = Activity("generate case vs controls figure",
                           executed = GIT_URL, 
                           used = c(UNREPEATED_PD_v_NONPD_SYN_ID, 
                                    REPEATED_PD_v_NONPD_SYN_ID)))
  unlink(title)

  
  
  ###################################
  ## Supplementary Figure 4
  ###################################
  metric <- "Auc"
  datC <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                            "/tap/matched/rr2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/walk/matched/rr2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/rest/matched/rr2")$metrics[[metric]], 
                     h5read(unrepeated_pd_vs_non_pd, 
                            "/voice/matched/rr2")$metrics[[metric]])
  
  datCP <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                             "/tap/matched/rr2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/walk/matched/rr2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/rest/matched/rr2")$metricsP[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/voice/matched/rr2")$metricsP[[metric]])
  
  metric <- "BalAcc"
  datC2 <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                             "/tap/matched/rr2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/walk/matched/rr2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/rest/matched/rr2")$metrics[[metric]], 
                      h5read(unrepeated_pd_vs_non_pd, 
                             "/voice/matched/rr2")$metrics[[metric]])
  
  datCP2 <- data.frame(h5read(unrepeated_pd_vs_non_pd, 
                              "/tap/matched/rr2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/walk/matched/rr2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/rest/matched/rr2")$metricsP[[metric]], 
                       h5read(unrepeated_pd_vs_non_pd, 
                              "/voice/matched/rr2")$metricsP[[metric]])
  
  
  #### Retrieve AUC vs Balanced Accuracy ####
  metric <- "auc"
  datS <- data.frame(h5read(repeated_pd_vs_non_pd, 
                            "tap/rr_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "walk/rr_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "rest/rr_unshuffled")[[metric]], 
                     h5read(repeated_pd_vs_non_pd, 
                            "voice/rr_unshuffled")[[metric]])
  datSP <- data.frame(h5read(repeated_pd_vs_non_pd, 
                             "tap/rr_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "walk/rr_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "rest/rr_shuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "voice/rr_shuffled")[[metric]])
  
  metric <- "bacc"
  datS2 <- data.frame(h5read(repeated_pd_vs_non_pd, 
                             "tap/rr_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "walk/rr_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "rest/rr_unshuffled")[[metric]], 
                      h5read(repeated_pd_vs_non_pd, 
                             "voice/rr_unshuffled")[[metric]])
  datSP2 <- data.frame(h5read(repeated_pd_vs_non_pd, 
                              "tap/rr_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "walk/rr_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "rest/rr_shuffled")[[metric]], 
                       h5read(repeated_pd_vs_non_pd, 
                              "voice/rr_shuffled")[[metric]])
  
  ### save figure ###
  title <- SUPPL_FIGURE_4
  png(paste(title), width = 2000, height = 2000, res = 200)
  get_supplementary_figures(datS, datSP, datS2, datSP2,
                            datC, datCP, datC2, datCP2)
  dev.off()
  
  #############################
  ## Store to Synapse
  #############################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS_SUPPL_4
  synStore(f, activity = Activity(
    "generate case vs controls figure",
    executed = GIT_URL, 
    used = c(UNREPEATED_PD_v_NONPD_SYN_ID, 
             REPEATED_PD_v_NONPD_SYN_ID)))
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

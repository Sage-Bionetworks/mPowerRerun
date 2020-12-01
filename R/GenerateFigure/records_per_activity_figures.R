#########################################################
#' This script is used to generate supplementary figure 7
#' that plots the record per activity
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
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
                   intermediate = get_intermediate_data_ref(),
                   processed = get_processed_features_ref())
SCRIPT_NAME <-  "records_per_activity_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
FIGURE_TITLE <- paste0("mPower_",
                       gsub(" ", "_", get("metadata")$user_group), 
                       "_supplementary_figure_7",".png")
ANNOTATIONS <- list(study = tolower(get("metadata")$study),
                    analysisType = "records per mPower tasks",
                    userSubset = tolower(get("metadata")$user_group), 
                    pipelineStep = "figures")

######################################################
## Helpers
#######################################################
GetNumRecords <- function(dat, sortOutput = FALSE) {
  ids <- as.character(unique(dat$healthCode))
  nsubj <- length(ids)
  nrec <- matrix(NA, nsubj, 1)
  colnames(nrec) <- "nrec"
  rownames(nrec) <- ids
  for (i in seq(nsubj)) {
    nrec[i] <- length(which(dat$healthCode == ids[i]))
  }
  if (sortOutput) {
    o <- order(nrec[, 1], decreasing = TRUE)
    nrec <- nrec[o,, drop = FALSE]
  }
  return(nrec)
} 


main <- function(){
  ## get data activity
  datTap <- read.csv(synGet(SYN_ID_REF$processed$tapping)$path, 
                     sep = "\t", stringsAsFactors = TRUE)
  datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, 
                     sep = "\t", stringsAsFactors = TRUE)
  datRes <- read.csv(synGet(SYN_ID_REF$processed$resting)$path,
                     sep = "\t", stringsAsFactors = TRUE)
  datWal <- read.csv(synGet(SYN_ID_REF$processed$walking)$path, 
                     sep = "\t", stringsAsFactors = TRUE)

  ## filter n of 1  
  datTap0 <- GetDataForNof1(datTap, beforeThr = 15, afterThr = 15)
  datRes0 <- GetDataForNof1(datRes, beforeThr = 15, afterThr = 15)
  datWal0 <- GetDataForNof1(datWal, beforeThr = 15, afterThr = 15)
  datVoi0 <- GetDataForNof1(datVoi, beforeThr = 15, afterThr = 15)
  
  ## count the number of records per healthCode
  nrecTap <- GetNumRecords(datTap0, sortOutput = TRUE)
  nrecRes <- GetNumRecords(datRes0, sortOutput = TRUE)
  nrecWal <- GetNumRecords(datWal0, sortOutput = TRUE)
  nrecVoi <- GetNumRecords(datVoi0, sortOutput = TRUE)
  
  summary(nrecTap[, 1])
  summary(nrecWal[, 1])
  summary(nrecRes[, 1])
  summary(nrecVoi[, 1])
  
  ######################################################
  ## Generate Figures
  #######################################################
  colorTap <- "#00A087FF"
  colorVoi <- "#7E6148FF"
  colorWal <- "#E64B36FF"
  colorRes <- "#8491B4FF"
  cexleg <- 1.5
  cext <- 1.3
  cexn <- 0.3
  
  
  title <- FIGURE_TITLE
  fname <- paste("", title, sep = "")
  png(title, width = 1500, height = 2000, res = 200)
  par(mfrow = c(4, 1), mar = c(6, 4, 0.5, 0.5)+0.1, mgp = c(2.5, 0.75, 0))
  barplot(t(nrecTap), las = 2, cex.names = cexn, ylim = c(0, 450), col = colorTap,
          main = "", ylab = "Number of records per participant")
  abline(h = 30, col = "red")
  legend("right", legend = "(a)", cex = cexleg, bty = "n")
  mtext(side = 3, "Tapping activity", cex = cext, line = -3)
  ####
  barplot(t(nrecWal), las = 2, cex.names = cexn, ylim = c(0, 450), col = colorWal,
          main = "", ylab = "Number of records per participant")
  abline(h = 30, col = "red")
  legend("right", legend = "(b)", cex = cexleg, bty = "n")
  mtext(side = 3, "Walk activity", cex = cext, line = -3)
  ####
  barplot(t(nrecRes), las = 2, cex.names = cexn, ylim = c(0, 450), col = colorRes,
          main = "", ylab = "Number of records per participant")
  abline(h = 30, col = "red")
  legend("right", legend = "(c)", cex = cexleg, bty = "n")
  mtext(side = 3, "Balance activity", cex = cext, line = -3)
  ####
  barplot(t(nrecVoi), las = 2, cex.names = cexn, ylim = c(0, 450), col = colorVoi,
          main = "", ylab = "Number of records per participant")
  abline(h = 30, col = "red")
  legend("right", legend = "(d)", cex = cexleg, bty = "n")
  mtext(side = 3, "Voice activity", cex = cext, line = -3)
  dev.off()
  
  ##############################
  ## Store to Synapse
  ##############################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate records per activity fig",
    executed = GIT_URL, 
    used = c(SYN_ID_REF$processed$resting,
             SYN_ID_REF$processed$tapping,
             SYN_ID_REF$processed$walking,
             SYN_ID_REF$processed$voice)))
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


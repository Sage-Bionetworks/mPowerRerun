######################################################
# This script will plot supplementary figure 11,
# which plot user medication activity vs time-of-day
#######################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
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
BEFORE_MEDICATION <- 15
AFTER_MEDICATION <- 15
thr <- 0.05
SYN_ID_REF <- list(figures = get_figure_ref(),
                   processed = get_processed_features_ref(),
                   healthcode = get_healthcode_ref())
SCRIPT_NAME <-  "user_medication_activity_vs_tod_figures.R"
FEATURE_LIST <- get_features()
GIT_URL <- getPermlink(getRepo(get("git")$repo), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
FIGURE_TITLE <- paste0("mPower_",
                       gsub(" ", "_", get("metadata")$user_group), 
                       "_supplementary_figure_11",".png")
ANNOTATIONS <- list(study = tolower(get("metadata")$study),
                    analysisType = "n of 1 Analysis",
                    analysisSubtype = "user medication activity vs tod",
                    userSubset = tolower(get("metadata")$user_group),
                    pipelineStep = "figures")

######################################################
## Helper Functions
#######################################################

# get difference between before and after time of day
BeforeAfterTodDifferences <- function(dat, mtm = "BH") {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids) 
  out <- data.frame(matrix(NA, nids, 3))
  colnames(out) <- c("after-before", "t-statistic", "pval")
  rownames(out) <- ids
  for (i in seq(nids)) {
    sdat <- dat[dat$healthCode == ids[i],]
    after <- sdat$tod[sdat$medTimepoint == "Just after Parkinson medication (at your best)"]
    before <- sdat$tod[sdat$medTimepoint == "Immediately before Parkinson medication"]
    aux <- t.test(x = after, y = before, alternative = "two.sided")
    out[i, 1] <- mean(after) - mean(before)
    out[i, 2] <- aux$statistic
    out[i, 3] <- aux$p.value
  }
  out <- out[order(out[, 2]),] ## order by t-statistic
  out$adjPval <- p.adjust(out[, 3], method = mtm)
  return(out)
}

# generate plots
GenerateAfterMinusBeforePlot <- function(out, 
                                         thr = 0.05, 
                                         main = "",
                                         cex.names = 0.5,
                                         cex.lab = 1,
                                         cex.main = 1) {
  sigIdx <- which(out$adjPval <= thr)
  cols <- rep("grey", nrow(out))
  cols[sigIdx] <- "darkred"
  barplot(out[, 2], col = cols, main = main, 
          ylab = "(average ToD after med.) - (average ToD before med.)",
          names.arg = rownames(out), las = 2, cex.names = cex.names,
          cex.lab = cex.lab)
}



main <- function(){
  datTap <- read.csv(synGet(SYN_ID_REF$processed$tapping)$path, 
                     stringsAsFactors = TRUE, sep = "\t") 
  datVoi <- read.csv(synGet(SYN_ID_REF$processed$voice)$path, 
                     stringsAsFactors = TRUE, sep = "\t")
  datRes <- read.csv(synGet(SYN_ID_REF$processed$resting)$path, 
                     stringsAsFactors = TRUE, sep = "\t")
  datWal <- read.csv(synGet(SYN_ID_REF$processed$walking)$path, 
                     stringsAsFactors = TRUE, sep = "\t")
  
  tapFeatures <- FEATURE_LIST$tapping
  walkFeatures <- FEATURE_LIST$walking
  restFeatures <- FEATURE_LIST$resting
  voiceFeatures <- FEATURE_LIST$voice
  
  healthCode <- read.delim(synGet(SYN_ID_REF$healthcode$n_of_one)$path, 
                           header = TRUE, stringsAsFactors = TRUE) 
  tzTapData <- healthCode %>% dplyr::filter(activity == "tapping")
  tzWalData <- healthCode %>% dplyr::filter(activity == "walking")
  tzResData <- healthCode %>% dplyr::filter(activity == "resting")
  tzVoiData <- healthCode %>% dplyr::filter(activity == "voice")
  
  #########################################################
  ## shape tapping data
  #########################################################
  datTap0 <- GetDataForNof1(datTap, BEFORE_MEDICATION, AFTER_MEDICATION)
  datTap0 <- IncludeUTCandLocalTimeVariables(datTap0, tzTapData) %>%
    filter(!is.na(tod))
  
  #########################################################
  ## shape walk data
  #########################################################
  datWal0 <- GetDataForNof1(datWal, BEFORE_MEDICATION, AFTER_MEDICATION)
  datWal0 <- IncludeUTCandLocalTimeVariables(datWal0, tzWalData) %>%
    filter(!is.na(tod))
  
  
  #########################################################
  ## shape rest data
  #########################################################
  datRes0 <- GetDataForNof1(datRes, BEFORE_MEDICATION, AFTER_MEDICATION)
  datRes0 <- IncludeUTCandLocalTimeVariables(datRes0, tzResData) %>%
    filter(!is.na(tod))
  
  
  
  #########################################################
  ## shape voice data
  #########################################################
  datVoi0 <- GetDataForNof1(datVoi, BEFORE_MEDICATION, AFTER_MEDICATION)
  datVoi0 <- IncludeUTCandLocalTimeVariables(datVoi0, tzVoiData) %>%
    filter(!is.na(tod))
  
  outTap <- BeforeAfterTodDifferences(dat = datTap0, mtm = "BH")
  outWal <- BeforeAfterTodDifferences(dat = datWal0, mtm = "BH")
  outRes <- BeforeAfterTodDifferences(dat = datRes0, mtm = "BH")
  outVoi <- BeforeAfterTodDifferences(dat = datVoi0, mtm = "BH")
  
  thr <- 0.05
  cex.names <- 0.5
  cex.lab <- 1
  cex.main <- 1.5
  cexleg <- 1.5
  
  title <- FIGURE_TITLE
  png(title, width = 2000, height = 2000, res = 200)
  
  par(mfrow = c(2, 2), mar = c(10, 4, 3, 0.5), mgp = c(2, 0.75, 0))
  GenerateAfterMinusBeforePlot(outTap, thr, main = "tapping", cex.names, cex.lab, cex.main)
  legend("topleft", legend = "(a)", bty = "n", cex = cexleg)
  GenerateAfterMinusBeforePlot(outWal, thr, main = "walk", cex.names, cex.lab, cex.main)
  legend("topleft", legend = "(b)", bty = "n", cex = cexleg)
  GenerateAfterMinusBeforePlot(outRes, thr, main = "balance", cex.names, cex.lab, cex.main)
  legend("topleft", legend = "(c)", bty = "n", cex = cexleg)
  GenerateAfterMinusBeforePlot(outVoi, thr, main = "voice", cex.names, cex.lab, cex.main)
  legend("topleft", legend = "(d)", bty = "n", cex = cexleg)
  dev.off()
  
  #########################################################
  ## store to synapse
  #########################################################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate user medication vs tod fig",
    executed = GIT_URL, 
    used = c(SYN_ID_REF$processed$tapping,
             SYN_ID_REF$processed$walking,
             SYN_ID_REF$processed$voice,
             SYN_ID_REF$processed$resting)))
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



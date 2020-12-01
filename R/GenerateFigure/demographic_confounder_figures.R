###################################################
#' This script is used to generate supplementary figure 
#' 14 to 31 of the mPower Paper
#' It will graph demographic confounders 
#' for age, gender, education
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
SCRIPT_NAME <-  "demographic_confounder_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
CORR_TBL_SYN_ID <- SYN_ID_REF$intermediate$confounder_corr
DCORR_TBL_SYN_ID <- SYN_ID_REF$intermediate$confounder_dcorr
ANNOTATIONS <- list(study = tolower(get("metadata")$study),
                    analysisType = "demographics confounders",
                    userSubset = tolower(get("metadata")$user_group), 
                    pipelineStep = "figures")


AdjustPvals <- function(x, method = "bonf") {
  -log10(apply(x, 2, p.adjust, method))
}

GeneratePatternFigure <- function(tapD, walD, resD, voiD, 
                                  plottitle = "",
                                  cexleg = 1.2,
                                  cext1 = 1.2, 
                                  cext2 = 1.5, 
                                  titlePos = -0.5,
                                  mylim = c(0.3, 0.7)) {
  # figure parameter
  colorTap <- "#00A087FF"
  colorVoi <- "#7E6148FF"
  colorWal <- "#E64B36FF"
  colorRes <- "#8491B4FF"
  cexleg <- 1.2
  mylim <- c(-0.3, 0.7)
  figpath <- ""
  idx1 <- 1:5
  idx2 <- 6:10
  colnames(tapD) <- rep(c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)"), 2)
  colnames(walD) <- rep(c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)"), 2)
  colnames(resD) <- rep(c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)"), 2)
  colnames(voiD) <- rep(c("cor(R,Y)", "cor(R,A)", "cor(A,Y)", "cor(R,Y|A)", "cor(R,A|Y)"), 2)
  par(mfrow = c(2, 4), mar = c(0.1, 3.0, 5.5, 0.5), mgp = c(2.0, 0.75, 0))
  
  #### a
  boxplot(tapD[, idx1], las = 2, ylim = mylim, main = "", 
          ylab = "correlation", 
          col = "white",
          border = colorTap, 
          xaxt = "n")
  mtext(side = 3, "tapping", 
        cex = cext1, 
        line = 0.1)
  abline(h = 0, col = "grey")
  legend("topright", legend = "(a)", cex = cexleg, bty = "n")
  
  #### b
  boxplot(walD[, idx1], las = 2, 
          ylim = mylim, main = "", 
          ylab = "correlation", 
          border = colorWal, 
          xaxt = "n",
          col = "white")
  mtext(side = 3, "walk", cex = cext1, line = 0.1)
  abline(h = 0, col = "grey")
  legend("topright", 
         legend = "(b)", 
         cex = cexleg, bty = "n")
  
  ### c
  boxplot(resD[, idx1], las = 2, 
          ylim = mylim, main = "", 
          ylab = "correlation", 
          col = "white",
          border = colorRes, xaxt = "n")
  mtext(side = 3, "balance", cex = cext1, line = 0.1)
  mtext(side = 3, plottitle, cex = cext2, line = 3, at = titlePos)
  abline(h = 0, col = "grey")
  legend("topright", legend = "(c)", cex = cexleg, bty = "n")
  
  ### d
  boxplot(voiD[, idx1], las = 2, 
          col = "white",
          ylim = mylim, main = "", 
          ylab = "correlation", 
          border = colorVoi, xaxt = "n")
  mtext(side = 3, "voice", cex = cext1, line = 0.1)
  abline(h = 0, col = "grey")
  legend("topright", legend = "(d)", cex = cexleg, bty = "n")
  ####
  par(mar = c(5.5, 3.0, 0.1, 0.5), mgp = c(2.0, 0.75, 0))
  ####
  mtm <- "BH"
  aux <- tapD[, idx2]
  aux[aux == 0] <- 1e-323
  apvals <- AdjustPvals(aux, mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  
  ### e
  boxplot(apvals, 
          las = 2, 
          col = "white",
          ylim = c(0, maxlim+1), 
          main = "", 
          ylab = "-log10(adjusted p-value)", 
          border = colorTap)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(e)", cex = cexleg, bty = "n")
  #
  aux <- walD[, idx2]
  aux[aux == 0] <- 1e-323
  apvals <- AdjustPvals(aux, mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  
  ### e
  boxplot(apvals, las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          col = "white", 
          ylab = "-log10(adjusted p-value)", 
          border = colorWal)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(f)", cex = cexleg, bty = "n")
  #
  aux <- resD[, idx2]
  aux[aux == 0] <- 1e-323
  apvals <- AdjustPvals(aux, mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  
  ### f
  boxplot(apvals, las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          col = "white", 
          ylab = "-log10(adjusted p-value)", 
          border = colorRes)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(g)", cex = cexleg, bty = "n")
  #
  aux <- voiD[, idx2]
  aux[aux == 0] <- 1e-323
  apvals <- AdjustPvals(aux, mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  
  ### g
  boxplot(apvals, las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          col = "white", 
          ylab = "-log10(adjusted p-value)", 
          border = colorVoi)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(h)", cex = cexleg, bty = "n")
}

CompareCorVsdCorFigure <- function(tapCor, 
                                   walCor, 
                                   resCor, 
                                   voiCor, 
                                   tapdCor, 
                                   waldCor, 
                                   resdCor, 
                                   voidCor,
                                   plottitle = "",
                                   cexleg = 1.2,
                                   cext1 = 1.2, 
                                   cext2 = 1.5, 
                                   titlePos = -0.5) {
  # figure parameter
  colorTap <- "#00A087FF"
  colorVoi <- "#7E6148FF"
  colorWal <- "#E64B36FF"
  colorRes <- "#8491B4FF"
  cexleg <- 1.2
  mylim <- c(-0.3, 0.7)
  figpath <- ""
  idx2 <- 6:10
  mtm <- "BH"
  par(mfrow = c(2, 4), mar = c(0.1, 3.0, 5.5, 0.5), mgp = c(2.0, 0.75, 0))
  ####
  apvals <- AdjustPvals(tapCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, las = 2, 
          col = "white",
          ylim = c(0, maxlim+1),
          main = "", 
          ylab = "-log10(adjusted p-value)", 
          border = colorTap, 
          xaxt = "n")
  mtext(side = 3, "tapping", cex = cext1, line = 0.1)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(a)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(walCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          ylab = "-log10(adjusted p-value)", 
          border = colorWal, xaxt = "n")
  mtext(side = 3, "walk", cex = cext1, line = 0.1)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(b)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(resCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          ylab = "-log10(adjusted p-value)", 
          border = colorRes, xaxt = "n")
  mtext(side = 3, "balance", cex = cext1, line = 0.1)
  mtext(side = 3, plottitle, cex = cext2, line = 3, at = titlePos)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(c)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(voiCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          ylab = "-log10(adjusted p-value)", 
          border = colorVoi, xaxt = "n")
  mtext(side = 3, "voice", cex = cext1, line = 0.1)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(d)", cex = cexleg, bty = "n")  
  ####
  par(mar = c(5.5, 3.0, 0.1, 0.5), mgp = c(2.0, 0.75, 0))
  ####
  apvals <- AdjustPvals(tapdCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "", 
          ylab = "-log10(adjusted p-value)", 
          border = colorTap)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(e)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(waldCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1),
          main = "",
          ylab = "-log10(adjusted p-value)", 
          border = colorWal)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(f)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(resdCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          ylab = "-log10(adjusted p-value)", 
          border = colorRes)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(g)", cex = cexleg, bty = "n")
  #
  apvals <- AdjustPvals(voidCor[, idx2], mtm)
  maxlim <- max(apvals, na.rm = TRUE)
  boxplot(apvals, 
          col = "white",
          las = 2, 
          ylim = c(0, maxlim+1), 
          main = "",
          ylab = "-log10(adjusted p-value)",
          border = colorVoi)
  abline(h = -log10(0.05), col = "red")
  legend("topright", legend = "(h)", cex = cexleg, bty = "n")
}

main <- function(){
  
  #######################################################
  ## read table
  #######################################################
  corr_tbl  <- H5Fopen(synGet(CORR_TBL_SYN_ID)$path)
  dcorr_tbl    <- H5Fopen(synGet(DCORR_TBL_SYN_ID)$path)
  dcorr_tbl1   <- H5Fopen(synGet(DCORR_TBL_SYN_ID)$path)
  
  #######################################################
  ## Random Forest Figures
  #######################################################
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_14",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/age_rf/corRf"), 
                        walD = h5read(corr_tbl, "walk/age_rf/corRf"), 
                        resD = h5read(corr_tbl, "rest/age_rf/corRf"), 
                        voiD = h5read(corr_tbl, "voice/age_rf/corRf"),
                        plottitle = "age - unmatched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## umatched gender ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_15",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/gender_rf/corRf"), 
                        walD = h5read(corr_tbl, "walk/gender_rf/corRf"), 
                        resD = h5read(corr_tbl, "rest/gender_rf/corRf"), 
                        voiD = h5read(corr_tbl, "voice/gender_rf/corRf"),
                        plottitle = "gender - unmatched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## umatched education ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_16",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/education_rf/corRf"), 
                        walD = h5read(corr_tbl, "walk/education_rf/corRf"), 
                        resD = h5read(corr_tbl, "rest/education_rf/corRf"), 
                        voiD = h5read(corr_tbl, "voice/education_rf/corRf"),
                        plottitle = "education - unmatched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  ##############################
  ## Ridge Regression Figures
  ##############################
  ## unmatched age ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_17",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/age_rr/corRr"), 
                        walD = h5read(corr_tbl, "walk/age_rr/corRr"), 
                        resD = h5read(corr_tbl, "rest/age_rr/corRr"), 
                        voiD = h5read(corr_tbl, "voice/age_rr/corRr"),
                        plottitle = "age - unmatched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## umatched gender ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_18",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/gender_rr/corRr"), 
                        walD = h5read(corr_tbl, "walk/gender_rr/corRr"), 
                        resD = h5read(corr_tbl, "rest/gender_rr/corRr"), 
                        voiD = h5read(corr_tbl, "voice/gender_rr/corRr"),
                        plottitle = "gender - unmatched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## umatched education ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_19",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/education_rr/corRr"), 
                        walD = h5read(corr_tbl, "walk/education_rr/corRr"), 
                        resD = h5read(corr_tbl, "rest/education_rr/corRr"), 
                        voiD = h5read(corr_tbl, "voice/education_rr/corRr"),
                        plottitle = "education - unmatched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  ## matched age ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_20",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/age_rf/corRfM"), 
                        walD = h5read(corr_tbl, "walk/age_rf/corRfM"), 
                        resD = h5read(corr_tbl, "rest/age_rf/corRfM"), 
                        voiD = h5read(corr_tbl, "voice/age_rf/corRfM"), 
                        plottitle = "age - matched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  ## matched gender ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_21",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/gender_rf/corRfM"),
                        walD = h5read(corr_tbl, "walk/gender_rf/corRfM"),
                        resD = h5read(corr_tbl, "rest/gender_rf/corRfM"),
                        voiD = h5read(corr_tbl, "voice/gender_rf/corRfM"),
                        plottitle = "gender - matched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate demo confounding figure",
    executed = GIT_URL, used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  ## matched education ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_22",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/education_rf/corRfM"),
                        walD = h5read(corr_tbl, "walk/education_rf/corRfM"),
                        resD = h5read(corr_tbl, "rest/education_rf/corRfM"),
                        voiD = h5read(corr_tbl, "voice/education_rf/corRfM"),
                        plottitle = "education - matched - random forest",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## matched age ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_23",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/age_rr/corRrM"), 
                        walD = h5read(corr_tbl, "walk/age_rr/corRrM"), 
                        resD = h5read(corr_tbl, "rest/age_rr/corRrM"), 
                        voiD = h5read(corr_tbl, "voice/age_rr/corRrM"), 
                        plottitle = "age - matched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## matched gender ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_24",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/gender_rr/corRrM"),
                        walD = h5read(corr_tbl, "walk/gender_rr/corRrM"),
                        resD = h5read(corr_tbl, "rest/gender_rr/corRrM"),
                        voiD = h5read(corr_tbl, "voice/gender_rr/corRrM"),
                        plottitle = "gender - matched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  ## matched education ##
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_25",".png")
  png(title, width = 2100, height = 1300, res = 200)
  GeneratePatternFigure(tapD = h5read(corr_tbl, "tap/education_rr/corRrM"),
                        walD = h5read(corr_tbl, "walk/education_rr/corRrM"),
                        resD = h5read(corr_tbl, "rest/education_rr/corRrM"),
                        voiD = h5read(corr_tbl, "voice/education_rr/corRrM"),
                        plottitle = "education - matched - ridge regression",
                        cext2 = 1.7,
                        mylim = mylim)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  #################################################
  ## Compare correlation vs distance correlation
  #################################################
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_26",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/age_rf/corRfM"),
                         walCor = h5read(corr_tbl, "walk/age_rf/corRfM"),
                         resCor = h5read(corr_tbl, "rest/age_rf/corRfM"),
                         voiCor = h5read(corr_tbl, "voice/age_rf/corRfM"),
                         tapdCor = h5read(dcorr_tbl, "tap/age_rf/corRf"),
                         waldCor = h5read(dcorr_tbl, "walk/age_rf/corRf"),
                         resdCor = h5read(dcorr_tbl, "rest/age_rf/corRf"),
                         voidCor = h5read(dcorr_tbl, "voice/age_rf/corRf"),
                         plottitle = "age - cor. vs distance cor. tests - random forest",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_27",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/gender_rf/corRfM"),
                         walCor = h5read(corr_tbl, "walk/gender_rf/corRfM"),
                         resCor = h5read(corr_tbl, "rest/gender_rf/corRfM"),
                         voiCor = h5read(corr_tbl, "voice/gender_rf/corRfM"),
                         tapdCor = h5read(dcorr_tbl, "tap/gender_rf/corRf"),
                         waldCor = h5read(dcorr_tbl, "walk/gender_rf/corRf"),
                         resdCor = h5read(dcorr_tbl, "rest/gender_rf/corRf"),
                         voidCor = h5read(dcorr_tbl, "voice/gender_rf/corRf"),
                         plottitle = "gender - cor. vs distance cor. tests - random forest",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_28",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/education_rf/corRfM"),
                         walCor = h5read(corr_tbl, "walk/education_rf/corRfM"),
                         resCor = h5read(corr_tbl, "rest/education_rf/corRfM"),
                         voiCor = h5read(corr_tbl, "voice/education_rf/corRfM"),
                         tapdCor = h5read(dcorr_tbl, "tap/education_rf/corRf"),
                         waldCor = h5read(dcorr_tbl, "walk/education_rf/corRf"),
                         resdCor = h5read(dcorr_tbl, "rest/education_rf/corRf"),
                         voidCor = h5read(dcorr_tbl, "voice/education_rf/corRf"),
                         plottitle = "education - cor. vs distance cor. tests - random forest",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_29",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/age_rr/corRrM"),
                         walCor = h5read(corr_tbl, "walk/age_rr/corRrM"),
                         resCor = h5read(corr_tbl, "rest/age_rr/corRrM"),
                         voiCor = h5read(corr_tbl, "voice/age_rr/corRrM"),
                         tapdCor = h5read(dcorr_tbl, "tap/age_rr/corRr"),
                         waldCor = h5read(dcorr_tbl, "walk/age_rr/corRr"),
                         resdCor = h5read(dcorr_tbl, "rest/age_rr/corRr"),
                         voidCor = h5read(dcorr_tbl, "voice/age_rr/corRr"),
                         plottitle = "age - cor. vs distance cor. tests - ridge regression",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_30",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/gender_rr/corRrM"),
                         walCor = h5read(corr_tbl, "walk/gender_rr/corRrM"),
                         resCor = h5read(corr_tbl, "rest/gender_rr/corRrM"),
                         voiCor = h5read(corr_tbl, "voice/gender_rr/corRrM"),
                         tapdCor = h5read(dcorr_tbl, "tap/gender_rr/corRr"),
                         waldCor = h5read(dcorr_tbl, "walk/gender_rr/corRr"),
                         resdCor = h5read(dcorr_tbl, "rest/gender_rr/corRr"),
                         voidCor = h5read(dcorr_tbl, "voice/gender_rr/corRr"),
                         plottitle = "gender - cor. vs distance cor. tests - ridge regression",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
  unlink(title)
  
  
  
  title <- paste0("mPower_",
                  gsub(" ", "_", get("metadata")$user_group), 
                  "_supplementary_figure_31",".png")
  png(title, width = 2100, height = 1300, res = 200)
  CompareCorVsdCorFigure(tapCor = h5read(corr_tbl, "tap/education_rr/corRrM"),
                         walCor = h5read(corr_tbl, "walk/education_rr/corRrM"),
                         resCor = h5read(corr_tbl, "rest/education_rr/corRrM"),
                         voiCor = h5read(corr_tbl, "voice/education_rr/corRrM"),
                         tapdCor = h5read(dcorr_tbl, "tap/education_rr/corRr"),
                         waldCor = h5read(dcorr_tbl, "walk/education_rr/corRr"),
                         resdCor = h5read(dcorr_tbl, "rest/education_rr/corRr"),
                         voidCor = h5read(dcorr_tbl, "voice/education_rr/corRr"),
                         plottitle = "education - cor. vs distance cor. tests - ridge regression",
                         cexleg = 1.2,
                         cext1 = 1.2, 
                         cext2 = 1.7, 
                         titlePos = -0.5)
  dev.off()
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "generate confounder figure",
    executed = GIT_URL, 
    used = c(CORR_TBL_SYN_ID, DCORR_TBL_SYN_ID)))
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








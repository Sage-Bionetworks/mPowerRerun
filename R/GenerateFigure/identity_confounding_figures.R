###################################################
#' This script is used to generate supplementary figure 2
#' of the mPower paper. It will graph 
#' identity confounding of record-wise analysis
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
SCRIPT_NAME <-  "identity_confounding_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
IDENTITY_CONFOUNDING_SYN_ID <- SYN_ID_REF$intermediate$identity_confounding
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
FIGURE_TITLE <- FIGURE_TITLE <- paste0("mPower_",
                                       gsub(" ", "_", get("metadata")$user_group), 
                                       "_supplementary_figure_2",".png")
ANNOTATIONS <- list(study = tolower(get("metadata")$study),
                    analysisType = "identity confounding",
                    userSubset = tolower(get("metadata")$user_group), 
                    pipelineStep = "figures")


######################################################
## Helper
#######################################################

#' function to retrieve disease recognition pvalue
#' by computing how many numbers of AUC is larger than threshold AUC
#'  @param DrRWS: AUC for each training-test split
#'  @param StatsRWS: metrics for each splits
#'  @param furthest.split: split that is chosen by computing which
#'  train-test split has the furthest AUC from collective median splits
#'  @return pvalue whether disease can be recognized
get.disease.recog.pval <- function(DrRWS, StatsRWS, furthest.split){
    auc.thresh <- StatsRWS %>% 
        slice(furthest.split) %>% 
        .$auc %>% .[[1]]
    split.auc.scores <- DrRWS %>% 
        data.frame(.) %>% 
        dplyr::select(
            glue::glue("split.", as.character(furthest.split))) %>% 
        unlist(., use.names = FALSE)
    pval <- sum(split.auc.scores > auc.thresh)/nrow(DrRWS)
    return(pval)
}


main <- function(){
    
    ### read permutation results
    identity.confounding.tbl <- H5Fopen(synGet(IDENTITY_CONFOUNDING_SYN_ID)$path)
    
    ### read data
    tapDrRWS <- h5read(identity.confounding.tbl, "tap/drRWS")
    tapStatsRWS <- h5read(identity.confounding.tbl, "tap/statsRWS")
    voiDrRWS <- h5read(identity.confounding.tbl, "voice/drRWS")
    voiStatsRWS <- h5read(identity.confounding.tbl, "voice/statsRWS")
    resDrRWS <- h5read(identity.confounding.tbl, "rest/drRWS")
    resStatsRWS <- h5read(identity.confounding.tbl, "rest/statsRWS")
    walDrRWS <- h5read(identity.confounding.tbl, "walk/drRWS")
    walStatsRWS <- h5read(identity.confounding.tbl, "walk/statsRWS")
    
    
    colnames(tapDrRWS) <- paste("split", 1:30, sep = " ")
    colnames(walDrRWS) <- paste("split", 1:30, sep = " ")
    colnames(resDrRWS) <- paste("split", 1:30, sep = " ")
    colnames(voiDrRWS) <- paste("split", 1:30, sep = " ")
    
    
    jtap <- which.min(abs(tapStatsRWS[, "auc"] - median(tapStatsRWS[, "auc"])))
    jwal <- which.min(abs(walStatsRWS[, "auc"] - median(walStatsRWS[, "auc"])))
    jres <- which.min(abs(resStatsRWS[, "auc"] - median(resStatsRWS[, "auc"])))
    jvoi <- which.min(abs(voiStatsRWS[, "auc"] - median(voiStatsRWS[, "auc"])))
    
    xaxis <- seq(0, 1, length.out = 1000)
    dtap <- dnorm(xaxis, 0.5, sqrt(tapStatsRWS[jtap, "approxVar"]))
    dwal <- dnorm(xaxis, 0.5, sqrt(walStatsRWS[jwal, "approxVar"]))
    dres <- dnorm(xaxis, 0.5, sqrt(resStatsRWS[jres, "approxVar"]))
    dvoi <- dnorm(xaxis, 0.5, sqrt(voiStatsRWS[jvoi, "approxVar"]))
    
    pptap <- pnorm(median(tapDrRWS[, jtap]), 0.5, sqrt(tapStatsRWS[jtap, "approxVar"]), lower.tail = FALSE)
    ppwal <- pnorm(median(walDrRWS[, jwal]), 0.5, sqrt(walStatsRWS[jwal, "approxVar"]), lower.tail = FALSE)
    ppres <- pnorm(median(resDrRWS[, jres]), 0.5, sqrt(resStatsRWS[jres, "approxVar"]), lower.tail = FALSE)
    ppvoi <- pnorm(median(voiDrRWS[, jvoi]), 0.5, sqrt(voiStatsRWS[jvoi, "approxVar"]), lower.tail = FALSE)
    
    disease.recog.pval <- list()
    disease.recog.pval$tapping <- get.disease.recog.pval(tapDrRWS, tapStatsRWS, jtap)
    disease.recog.pval$resting <- get.disease.recog.pval(resDrRWS, resStatsRWS, jres)
    disease.recog.pval$walking <- get.disease.recog.pval(walDrRWS, walStatsRWS, jwal)
    disease.recog.pval$voice <- get.disease.recog.pval(voiDrRWS, voiStatsRWS, jvoi)
    
    
    
    xlim1 <- c(0, 1)
    ylim1 <- c(0, 90)
    ylim2 <- c(0.6, 0.9)
    cd <- 1
    alpha <- 0.5
    lwd1 <- 1
    ca <- 1.2
    cl <- 1.2
    cm <- 1.75
    ca2 <- 0.6
    
    mat <- matrix(c(1, 1, 1, 2,
                    3, 3, 3, 4,
                    5, 5, 5, 6,
                    7, 7, 7, 8), 
                  4, 4, byrow = TRUE)
    
    title <- FIGURE_TITLE
    png(title, width = 2000, height = 1800, res = 200)
    layout(mat)
    par(mar = c(4.5, 3, 1.5, 1), mgp = c(2, 0.75, 0))
    hist(tapDrRWS[, jtap], probability = TRUE, 
         col = rgb(0, 0, 1, alpha), xlim = xlim1, xlab = "AUC", 
         main = "tapping", ylim = ylim1, cex.axis  = ca, cex.lab = cl, cex.main = cm)
    abline(v = tapStatsRWS[jtap, "auc"], col = "brown", lwd = lwd1)
    abline(v = median(tapDrRWS[,jtap]), col = "black", lwd = lwd1)
    lines(xaxis, dtap, lty = 1, col = "black", lwd = lwd1)
    mtext("(a)", side = 3, at = 0, cex = ca)
    legend("topleft", 
           legend = c(glue::glue("disease recognition permutation p-value: ~{pvalue}", 
                                 pvalue = as.character(disease.recog.pval$tapping)), 
                      "pseudo p-value: ~0"), 
           text.col = c("brown", "black"), bty = "n", cex = ca)
    boxplot(tapDrRWS, las = 3, border = rgb(0, 0, 1, alpha), ylab = "AUC", xaxt = "n",
            ylim = c(0.75, 0.87), main = "tapping", cex.axis  = ca, cex.lab = cl)
    points(tapStatsRWS[, "auc"], pch = 20, col = "brown", cex = cd)
    axis(side = 1, at = seq(30), labels = paste("split", 1:30, sep = " "), las = 2, cex.axis = ca2)
    mtext("(e)", side = 3, at = 1, cex = ca)
    ####
    hist(walDrRWS[, jwal], probability = TRUE, 
         col = rgb(0, 0, 1, alpha), xlim = xlim1, xlab = "AUC", 
         main = "walk", ylim = ylim1, cex.axis  = ca, cex.lab = cl, cex.main = cm)
    abline(v = walStatsRWS[jwal, "auc"], col = "brown", lwd = lwd1)
    abline(v = median(walDrRWS[,jwal]), col = "black", lwd = lwd1)
    lines(xaxis, dwal, lty = 1, col = "black", lwd = lwd1)
    mtext("(b)", side = 3, at = 0, cex = ca)
    legend("topleft", legend = c(glue::glue("disease recognition permutation p-value: ~{pvalue}", 
                                            pvalue = as.character(disease.recog.pval$walking)), 
                                 "pseudo p-value: ~0"), 
           text.col = c("brown", "black"), bty = "n", cex = ca)
    boxplot(walDrRWS, las = 3, border = rgb(0, 0, 1, alpha), ylab = "AUC", xaxt = "n", 
            ylim = c(0.9, 0.95), main = "walk", cex.axis  = ca, cex.lab = cl)
    points(walStatsRWS[, "auc"], pch = 20, col = "brown", cex = cd)
    axis(side = 1, at = seq(30), labels = paste("split", 1:30, sep = " "), las = 2, cex.axis = ca2)
    mtext("(f)", side = 3, at = 1, cex = ca)
    ####
    hist(resDrRWS[, jres], probability = TRUE, 
         col = rgb(0, 0, 1, alpha), xlim = xlim1, xlab = "AUC", 
         main = "balance", ylim = ylim1, cex.axis  = ca, cex.lab = cl, cex.main = cm)
    abline(v = resStatsRWS[jres, "auc"], col = "brown", lwd = lwd1)
    abline(v = median(resDrRWS[,jres]), col = "black", lwd = lwd1)
    lines(xaxis, dres, lty = 1, col = "black", lwd = lwd1)
    mtext("(c)", side = 3, at = 0, cex = ca)
    legend("topleft", legend = c(glue::glue("disease recognition permutation p-value: ~{pvalue}", 
                                            pvalue = as.character(disease.recog.pval$resting)), 
                                 "pseudo p-value: ~0"), 
           text.col = c("brown", "black"), bty = "n", cex = ca)
    boxplot(resDrRWS, las = 3, border = rgb(0, 0, 1, alpha), ylab = "AUC", xaxt = "n", 
            ylim = c(0.65, 0.77), main = "rest", cex.axis  = ca, cex.lab = cl)
    points(resStatsRWS[, "auc"], pch = 20, col = "brown", cex = cd)
    axis(side = 1, at = seq(30), labels = paste("split", 1:30, sep = " "), las = 2, cex.axis = ca2)
    mtext("(g)", side = 3, at = 1, cex = ca)
    ####
    hist(voiDrRWS[, jvoi], probability = TRUE, 
         col = rgb(0, 0, 1, alpha), xlim = xlim1, xlab = "AUC", 
         main = "voice", ylim = ylim1, cex.axis  = ca, cex.lab = cl, cex.main = cm)
    abline(v = voiStatsRWS[jvoi, "auc"], col = "brown", lwd = lwd1)
    abline(v = median(voiDrRWS[,jvoi]), col = "black", lwd = lwd1)
    lines(xaxis, dvoi, lty = 1, col = "black", lwd = lwd1)
    mtext("(d)", side = 3, at = 0, cex = ca)
    legend("topleft", legend = c(glue::glue("disease recognition permutation p-value: ~{pvalue}", 
                                            pvalue = as.character(disease.recog.pval$voice)), 
                                 "pseudo p-value: ~0"), 
           text.col = c("brown", "black"), bty = "n", cex = ca)
    boxplot(voiDrRWS, las = 3, border = rgb(0, 0, 1, alpha), ylab = "AUC", xaxt = "n",
            ylim = c(0.82, 0.88), main = "voice", cex.axis  = ca, cex.lab = cl)
    points(voiStatsRWS[, "auc"], pch = 20, col = "brown", cex = cd)
    axis(side = 1, at = seq(30), labels = paste("split", 1:30, sep = " "), las = 2, cex.axis = ca2)
    mtext("(h)", side = 3, at = 1, cex = ca)
    dev.off()
    
    #############################
    ## Store to Synapse
    #############################
    f <- synapser::File(title, 
                        parent=FIGURE_OUTPUT_SYN_ID)
    f$annotations <- ANNOTATIONS
    synStore(f, activity = Activity(
        "generate identity confounding figures",
        executed = GIT_URL, 
        used = c(IDENTITY_CONFOUNDING_SYN_ID)))
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









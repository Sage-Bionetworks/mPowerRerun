######################################################
#' script to generate supplementary figure 10 
#' of the mPower paper, it will generate
#' whether user are adhering to the walking tasks
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
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
SYN_ID_REF <- list(figures = get_figure_ref(),
                   intermediate = get_intermediate_data_ref())
SCRIPT_NAME <-  "walking_adherence_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
ADHERENCE_SYN_ID <- SYN_ID_REF$intermediate$check_walking_adherence
FIGURE_TITLE <- paste0("mPower_",
  gsub(" ", "_", get("metadata")$user_group), 
  "_supplementary_figure_10",".png")
ANNOTATIONS <- list(
  analysisType = "walking task adherence",
  userSubset = tolower(get("metadata")$user_group), 
  pipelineStep = "figures")

main <- function(){
  ######################################################
  ## read data for adherence checking
  #######################################################
  dscores <- read.delim(synGet(ADHERENCE_SYN_ID)$path, 
                        sep = "\t", stringsAsFactors = FALSE)
  
  ######################################################
  ## Generate Figure
  #######################################################
  tb <- table(na.omit(dscores %>% filter(vertical != "NA"))$vertical)
  pc <- 100*tb/sum(tb)
  figpath <- ""
  myylim <- c(0, 20000)
  xlab <- "gravity acceleration"
  cexleg <- 2
  cm <- 1.75
  ca <- 1.5
  cl <- 1.5
  mat <- matrix(c(1, 1, 2,
                  3, 4, 5), 
                nr = 2, 
                nc = 3, 
                byrow = TRUE)
  
  title <- FIGURE_TITLE
  png(paste(figpath, title, sep = ""), width = 2400, height = 1800, res = 200)
  
  layout(mat)
  par(mar = c(4, 3.75, 1, 0.5), mgp = c(2.5, 0.75, 0))
  plot(seq(3), seq(3), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  barplot(pc, col = c("darkorange", "darkorange", "darkblue", "darkblue", "darkred", "darkred"),
          ylab = "percentage (%)", cex.names = 1.8, names.arg = c("-x", "+x", "-y", "+y", "-z", "+z"),
          cex.axis = ca, cex.lab = cl)
  legend("topleft", legend = "(b)", bty = "n", cex = cexleg)
  hist(dscores$median.x, ylim = myylim, xlab = xlab, cex.main = cm,
       main = "x-axis", col = "darkorange", cex.lab = cl, cex.axis = ca)
  legend("topleft", legend = "(c)", bty = "n", cex = cexleg)
  hist(dscores$median.y, ylim = myylim, xlab = xlab, cex.main = cm,
       main = "y-axis", col = "darkblue", cex.lab = cl, cex.axis = ca)
  legend("topleft", legend = "(d)", bty = "n", cex = cexleg)
  hist(dscores$median.z, ylim = myylim, xlab = xlab, cex.main = cm,
       main = "z-axis", col = "darkred", cex.lab = cl, cex.axis = ca)
  legend("topleft", legend = "(e)", bty = "n", cex = cexleg)
  dev.off()
  
  ######################################################
  ## Store Features
  #######################################################
  f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
  f$annotations <- ANNOTATIONS
  synStore(
    f, activity = Activity(
      "check walking adherence",
      executed = GIT_URL, used = c(ADHERENCE_SYN_ID)))
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










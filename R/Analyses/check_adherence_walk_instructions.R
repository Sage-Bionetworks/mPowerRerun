#################################################
#' This script is used for checking the adherence of
#' mPower users towards the walking test, which checks
#' whether they are correctly orienting their phones
#' in their respective pockets
#' 
#' @author: Elias Chaibub Neto, Aryton Tediarjo
#' @author_email: echaibub@synapse.org, aryton.tediarjo@sagebase.org
#################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(rhdf5)
library(doMC)
source("R/utils/populationAnalysisUtils.R")
source("R/utils/projectUtils.R")
source("R/utils/initializeVariables.R")

#######################################################
## Configuration
#######################################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))
registerDoMC(detectCores())


#######################################################
## Instantiate Variables and Reference IDs
#######################################################
SCRIPT_NAME <- "check_adherence_walk_instructions.R"
SYN_ID_REF <- list(intermediate = get_intermediate_data_ref(),
                   table = get_synapse_table_ref())
OUTPUT_DATA <- MODEL_OUTPUT <- paste0(
  "walking_task_adherence_", 
  gsub(" ", "_", get("metadata")$user_group), ".tsv")
WALK_ACTIVITY_TABLE_SYN_ID <-SYN_ID_REF$table$walking
OUTPUT_FOLDER_ID <- SYN_ID_REF$intermediate$output_folder
GIT_URL <- getPermlink(
  getRepo(get("git")$repo,
          ref="branch", 
          refName=get("git")$branch), 
  repositoryPath = file.path('R/Analyses',  SCRIPT_NAME))
ANNOTATIONS <- list(
  analysisType = "walking task adherence",
  dataSubtype = "processed",
  userSubset = get("metadata")$user_group,
  study = get("metadata")$study,
  pipelineStep= "intermediary data")


#' Get Direction based on mPower Walking Data
GetMainDirection <- function(metaData) {
  MainDirection <- function(dat) {
    mse1 <- mean((dat$x - 1)^2, na.rm = TRUE)
    mse2 <- mean((dat$x + 1)^2, na.rm = TRUE)
    mse3 <- mean((dat$y - 1)^2, na.rm = TRUE)
    mse4 <- mean((dat$y + 1)^2, na.rm = TRUE)
    mse5 <- mean((dat$z - 1)^2, na.rm = TRUE)
    mse6 <- mean((dat$z + 1)^2, na.rm = TRUE)
    aux <- which.min(c(mse1, mse2, mse3, mse4, mse5, mse6))
    medianAccel <- apply(dat, 2, median, na.rm = TRUE)
    if (aux == 1) {
      vertical <- "x+"
    }
    if (aux == 3) {
      vertical <- "y+"
    }
    if (aux == 5) {
      vertical <- "z+"
    }
    if (aux == 2) {
      vertical <- "x-"
    }
    if (aux == 4) {
      vertical <- "y-"
    }
    if (aux == 6) {
      vertical <- "z-"
    }
    list(vertical = vertical, 
         medianAccel = medianAccel)
  }
  nrec <- nrow(metaData)
  out <- data.frame(matrix(NA, nrec, 6))
  colnames(out) <- c("recordId", "healthCode", "vertical", 
                     "median.x", "median.y", "median.z")
  out$recordId <- metaData$recordId
  out$healthCode <- metaData$healthCode
  for (i in seq(nrec)) {
    cat(i, "\n")
    if (!is.na(metaData$deviceMotion_walking_outbound.fileLocation.items[i])) {
      dat <- jsonlite::fromJSON(as.character(
        metaData$deviceMotion_walking_outbound.fileLocation.items[i]))
      gdat <- dat$gravity
      aux <- MainDirection(gdat)
      out[i, "vertical"] <- aux$vertical
      out[i, 4:6] <- aux$medianAccel
    }
  }
  return(out)
}

main <- function(){
  
  #######################################################
  ## Data Processing
  #######################################################
  walk.tbl.syn <- synTableQuery(sprintf("SELECT * FROM %s", WALK_ACTIVITY_TABLE_SYN_ID))
  walk.tbl <- as.data.frame(walk.tbl.syn)
  columnsToDownload <- c("deviceMotion_walking_outbound.json.items") 
  walk.json.loc <- lapply(columnsToDownload, function(col.name){
    tbl.files = synapser::synDownloadTableColumns(walk.tbl.syn, col.name) %>%
      lapply(function(x) data.frame(V1 = x)) %>% 
      data.table::rbindlist(idcol = col.name) %>% 
      plyr::rename(c('V1' = gsub('.json','.fileLocation', col.name)))
  })
  
  walk.tbl$deviceMotion_walking_outbound.json.items <- as.character(
    walk.tbl$deviceMotion_walking_outbound.json.items)
  walk.tbl.meta <- data.table::rbindlist(list(walk.tbl %>%
                                                dplyr::left_join(
                                                  do.call(cbind, walk.json.loc))),
                                         use.names = T, 
                                         fill = T) %>% as.data.frame
  walk.tbl.meta$deviceMotion_walking_outbound.fileLocation.items <- as.character(
    walk.tbl.meta$deviceMotion_walking_outbound.fileLocation.items)
  
  dscores <- GetMainDirection(metaData = walk.tbl.meta) %>%
    write.table(., OUTPUT_DATA, sep="\t", row.names=F, quote=F)
  
  
  #######################################################
  ## Store Results to Synapse
  #######################################################
  f <- synapser::File(OUTPUT_DATA, OUTPUT_FOLDER_ID)
  f$annotations <- ANNOTATIONS
  synStore(f, activity = Activity(
    "Run Walk Adherence Information",
    executed = GIT_URL, 
    used = c(WALK_ACTIVITY_TABLE_SYN_ID)))
  unlink(OUTPUT_DATA)
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




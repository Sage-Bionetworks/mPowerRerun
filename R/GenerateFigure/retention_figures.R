######################################################
# This script will plot supplementary figure 1
# which plots user retentions 
#######################################################
library(install.load)
library(data.table)
library(tidyverse)
library(plyr)
library(survival)
library(survminer)
library(lubridate)
library(githubr)
install_load("gdata", "synapser", "jsonlite", "stringr")
install_load("gridExtra", "pheatmap", "printr", "ggthemes", "anytime")
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")
options(stringsAsFactors = F)


###################################
### configure credentials
####################################
synLogin()
config::get()
setGithubToken(
  readLines(get("git")$path))


##########################################################
### Get Reference and Final Cutoff Date for all activities
#########################################################
ENGAGEMENT_DATA <- 'syn17140076'
SCRIPT_NAME <- "retention_figures.R"
SYN_ID_REF <- list(raw = get_raw_features_ref(),
                   figures = get_figure_ref())
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
GIT_URL <- getPermlink(
  getRepo(get("git")$repo,
          ref="branch", 
          refName=get("git")$branch), 
  repositoryPath = file.path('R/GenerateFigure',  SCRIPT_NAME))
TAP_TBL_SYN_ID <- SYN_ID_REF$raw$tapping
VOICE_TBL_SYN_ID <- SYN_ID_REF$raw$voice
REST_TBL_SYN_ID <- SYN_ID_REF$raw$resting
WALK_TBL_SYN_ID <- SYN_ID_REF$raw$walking
ANNOTATIONS <- list(analysisType = "retention analysis",
                    userSubset = tolower(get("metadata")$user_group),
                    pipelineStep = "figures")

#' create logger for pipeline
sink('pipeline.log', append = TRUE)
cat(paste0(
  "[",Sys.time(), "]", " Running ", SCRIPT_NAME), "\n\n")
sink()

#addTime
fixCols <- function(df){
  df <- df %>% dplyr::mutate(createdOn = lubridate::ymd_hms(createdOn),
                      dayOfWeek = weekdays(createdOn),
                      weekend = ifelse(dayOfWeek %in% c('Sunday', 'Saturday'), T, F),
                      createdOnLocalTime = lubridate::ymd_hms(createdOnLocalTime),
                      hod_local = lubridate::hour(createdOnLocalTime) + 
                        lubridate::minute(createdOnLocalTime)/60,
                      participantWeek = participantWeek + 1,
                      participantDay = participantDay + 1 )
  return(df)
}

### The following data is downloaded from the engagement analysis project 
### https://github.com/Sage-Bionetworks/mhealth-engagement-analysis/blob/master/dataLoaders/loadData.R

##################
### mPower
##################
get_mpower_engagement_data <- function(cutoff.date){
  engagement <- fread(synGet(ENGAGEMENT_DATA)$path, 
                      data.table = F) %>%
    fixCols(.) %>% 
    dplyr::mutate(localDate = as.Date(lubridate::ymd_hms(createdOnLocalTime))) %>%
    dplyr::filter(localDate <= cutoff.date)
  return(engagement)
}
get_mpower_mdata <- function() {
  demo <- as.data.frame(synTableQuery(sprintf("SELECT * FROM %s", 
                                              get("synapse_tables")$demo)))
  colnames(demo) <- gsub("_|-", ".", names(demo))
  if("inferred.diagnosis" %in% names(demo)){
    PD <- demo$inferred.diagnosis
    demo <- demo %>% 
      dplyr::select(-c(professional.diagnosis, inferred.diagnosis)) %>% 
      mutate(professional.diagnosis = PD) %>% 
      filter(dataGroups %in% c("parkinson", "control", NA))}
  demo <- demo %>%
    mutate(caseStatus = professional.diagnosis) %>% 
    dplyr::select("healthCode", "appVersion", "phoneInfo", 
                  "education", "employment", "gender", 
                  "maritalStatus", "race", "caseStatus")
  return(demo)
}


##################################################
### Process Engagement Data
##################################################
cutoff.date <- lubridate::as_datetime(max(map_dbl(list(
  tap = TAP_TBL_SYN_ID,
  voice = VOICE_TBL_SYN_ID,
  walk =  WALK_TBL_SYN_ID,
  rest =  REST_TBL_SYN_ID),
  function(x){
    last.date <- max((read.csv(synGet(x)$path, sep = "\t") %>%
                        dplyr::mutate(createdOn = lubridate::ymd_hms(.$createdOn)) %>%
                        dplyr::filter(createdOn < lubridate::now()))$createdOn)
    return(last.date)})))

mpower_mdata <- get_mpower_mdata()
mpower <- get_mpower_engagement_data(cutoff.date)
mpower_user_stats <- mpower %>% 
  dplyr::group_by(healthCode) %>%
  dplyr::summarise(lastWeekinStudy = max(participantWeek),
                   totalDaysActive = n_distinct(lubridate::date(createdOn)),
                   firstDay = min(lubridate::date(createdOn), na.rm = T),
                   lastDay = max(lubridate::date(createdOn), na.rm = T),
                   duration_in_study = as.numeric(
                     max(lubridate::date(createdOn)) - min(lubridate::date(createdOn)) + 1)) %>%
  dplyr::mutate(study = "mpower")
mpower_user_stats <- base::merge(mpower_user_stats, 
                                 mpower_mdata %>% 
                                   dplyr::select(healthCode, caseStatus), all.x = T)
censor <- rep(1, nrow(mpower_user_stats)) 
fit.test <- survdiff(
  survival::Surv(time=duration_in_study, 
                 event=censor, type = "right") ~ caseStatus, data = mpower_user_stats)

fit.plot <- survfit(
  survival::Surv(time=duration_in_study, 
                 event=censor, type = "right") ~ caseStatus, data = mpower_user_stats)
stats <- as.data.frame(summary(fit.plot)$table)
renameCols <- c(LCL='0.95LCL', UCL='0.95UCL')

##################################################
### Generate Plot
##################################################
title <- paste0("mPower_",
                gsub(" ", "_", get("metadata")$user_group), 
                "_supplementary_figure_1",".png")
stats <- stats %>% 
  tibble::rownames_to_column(var='strata') %>%
  dplyr::rename(!!renameCols) %>%
  dplyr::mutate(CI=paste0(LCL,'-',UCL)) %>%
  dplyr::rename(N=records) %>%
  dplyr::select(strata, N, median, CI) %>%
  dplyr::mutate(strata = gsub('study=','', strata))
write.table(stats, file=title, 
            sep="\t", quote=F)

p1 <- survminer::ggsurvplot(fit.plot, pval = F, conf.int = T, 
                 palette = c('#3C5488FF', '#DC0000FF'),
                 xlab = "Duration in study ",  
                 risk.table = F,
                 risk.table.height = 0.3,
                 risk.table.y.text = FALSE,
                 legend.labs = c('nonPD', 'PD'),
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(base_size = 15))
ggsave(title, height = 6, width = 6, units="in", dpi=250, device = "png")

######################################################
## Store Features
#######################################################
f <- synapser::File(title, parent=FIGURE_OUTPUT_SYN_ID)
f$annotations <- ANNOTATIONS
synStore(f, activity = Activity(
  "generate retention figures",
  executed = GIT_URL, 
  used = c(ENGAGEMENT_DATA,
          get("demographics")$table,
          WALK_TBL_SYN_ID,
          TAP_TBL_SYN_ID,
          VOICE_TBL_SYN_ID,
          REST_TBL_SYN_ID)))
unlink(title)

sink('pipeline.log', append = TRUE)
cat(paste0("[",Sys.time(), "]", 
           " Done Generating Retention Figures"), 
    "\n\n")
sink()



  







######################################################
## script for healthcode manual conversion from RData
## as matched cohort used for 6 months paper is unable 
## to be replicated in the rerunning, thus it will be
## manually replicated using this script (only used for the paper rerun)
## Note: only use to replicate the exact figure in the paper
######################################################
library(synapser)
library(tidyverse)
library(plyr)
library(purrr)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

######################################################
## Configuration
######################################################
synLogin()
config::get()
setGithubToken(
    readLines(get("git")$path))

######################################################
## global variable
######################################################
SCRIPT_NAME <-  "convertHCfromRData.R"
SYN_ID_REF <- list(healthcode = get_healthcode_ref(),
                   table = get_synapse_table_ref())
OUTPUT_SYN_ID <- SYN_ID_REF$healthcode$output_folder
DEMOGRAPHICS_TBL_SYN_ID <- SYN_ID_REF$table$demo
IDENTITY_CONFOUNDING_MATCHED_DATA <- "syn11461177"
CASE_VS_CONTROLS_MATCHED_DATA <- "syn15355694"
TAP_TZ_DATA <- "syn8113978"
REST_TZ_DATA <- "syn8113984"
WALK_TZ_DATA <- "syn8113982"
VOICE_TZ_DATA <- "syn8113980"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/FilterHealthCode", SCRIPT_NAME))
OUTPUT_FILE <- list()
OUTPUT_FILE$pd_case_v_control <- paste0(
    "PD_case_vs_controls_matched_cohort_",
    get("metadata")$time_freeze,".tsv")
OUTPUT_FILE$identity_confounding <- paste0(
    "identity_confounding_matched_cohort_",
    get("metadata")$time_freeze,".tsv")
OUTPUT_FILE$N_of_1 <- paste0(
    "Nof1_filtered_cohort_",
    get("metadata")$time_freeze,".tsv")

######################################################
## Helpers
######################################################
get.demo <- function(){
    demo <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s", SYN_ID_REF$table$demo)))
    colnames(demo) <- gsub("_|-", ".", names(demo))
    if("inferred.diagnosis" %in% names(demo)){
        demo <- demo %>% 
            mutate(PD = demo$inferred.diagnosis) %>%
            dplyr::select(-c(professional.diagnosis, inferred.diagnosis)) %>%
            filter(dataGroups %in% c("parkinson", "control", NA))
    }else{
        demo <- demo %>% 
            dplyr::rename("PD" = "professional.diagnosis")
    }
    ## clean demographics data
    demo <- demo %>%
        dplyr::select(healthCode, education, employment,
                      are.caretaker, deep.brain.stimulation,
                      education, gender, health.history,
                      healthcare.provider, home.usage, maritalStatus,
                      medical.usage, medical.usage.yesterday, 
                      past.participation, phone.usage, 
                      PD, race, smartphone, age,
                      smoked, surgery, video.usage) %>%
        dplyr::filter(!is.infinite(age) | age <= 110) %>%
        plyr::ddply(.(healthCode), .fun = function(x){
            x$age = mean(x$age, na.rm = TRUE)
            x$PD = ifelse(
                length(unique(x$PD)) > 1, 
                NA, unique(x$PD))
            return(x[1,])})
    return(demo)
}

get.identity.confounding.match.data <- function(){
    demo <- get.demo()
    matched.hc <- list()
    matched <- read.delim(
        synGet(IDENTITY_CONFOUNDING_MATCHED_DATA)$path, 
        stringsAsFactors = FALSE)
    matched.hc$tapping  <- matched %>% filter(Assay == "Tapping") %>% list(.$X1, .$X2)
    matched.hc$voice <- matched %>% filter(Assay == "Voice") %>% list(.$X1, .$X2)
    matched.hc$walking <- matched %>% filter(Assay == "Walking") %>% list(.$X1, .$X2)
    matched.hc$resting <- matched %>% filter(Assay == "Rest") %>% list(.$X1, .$X2)
    matched.hc$memory <- matched %>% filter(Assay == "Memory") %>% list(.$X1, .$X2)
    matched.hc <- purrr::map(matched.hc, function(x){
        list(x[[2]], x[[3]]) %>% unlist()})
    healthcode <- purrr::map(names(matched.hc), 
                             function(activity){
                                 return(demo %>% 
                                            dplyr::filter(healthCode %in% matched.hc[[activity]]) %>%
                                            dplyr::mutate(activity = activity) %>%
                                            dplyr::select(healthCode, age, gender, activity))}) %>% purrr::reduce(., rbind)
    return(healthcode)
}

get.case.vs.controls.match.data <- function(){
    load(synGet(CASE_VS_CONTROLS_MATCHED_DATA)$path)
    demo <- get.demo()
    matched.hc <- list()
    matched.hc$tapping <- tapHC 
    matched.hc$resting <- resHC 
    matched.hc$voice <- voiHC 
    matched.hc$walking <- walHC 
    healthcode <- purrr::map(names(matched.hc), 
               function(activity){
                   return(demo %>% 
                              dplyr::filter(healthCode %in% matched.hc[[activity]]) %>%
                              dplyr::mutate(activity = activity) %>%
                              dplyr::select(healthCode, age, gender, activity))}) %>% purrr::reduce(., rbind)
    return(healthcode)
}

get.N.of.one.hc <- function(){
    hc.list <- list()
    hc.list$tapping <- read.delim(
        synGet(TAP_TZ_DATA)$path, header = TRUE) %>%
        dplyr::mutate(activity = "tapping")
    hc.list$resting <- read.delim(
        synGet(REST_TZ_DATA)$path, header = TRUE) %>%
        dplyr::mutate(activity = "resting")
    hc.list$walking <- read.delim(
        synGet(WALK_TZ_DATA)$path, header = TRUE) %>%
        dplyr::mutate(activity = "walking")
    hc.list$voice <- read.delim(
        synGet(VOICE_TZ_DATA)$path, header = TRUE) %>%
        dplyr::mutate(activity = "voice")
    return(hc.list %>% 
               purrr::reduce(., rbind) %>%
               dplyr::select(-Enter_State)) 
}


main <- function(){
    ## get data
    pd.case.vs.controls <- get.case.vs.controls.match.data()
    identity.confounding <- get.identity.confounding.match.data()
    n.of.one <- get.N.of.one.hc()
    
    ## store matched data to tsv
    pd.case.vs.controls %>% write.table(
        ., OUTPUT_FILE$pd_case_v_control, sep="\t", row.names=F, quote=F)
    identity.confounding %>% write.table(
        ., OUTPUT_FILE$identity_confounding, sep="\t", row.names=F, quote=F)
    n.of.one %>% write.table(
        ., OUTPUT_FILE$N_of_1, sep="\t", row.names=F, quote=F)
    
    ## store to synapse with the correct annotation
    file.entity<- synapser::File(OUTPUT_FILE$pd_case_v_control, 
                                 parent = OUTPUT_SYN_ID)
    file.entity$annotations <- list(
        analysisType = "case vs controls",
        pipelineStep = "healthcode subsampling",
        userSubset = get("metadata")$user_group
    )
    synStore(
        file.entity, activity = Activity(
        'replicate matched RData',
        executed = GIT_URL,
        used = c(CASE_VS_CONTROLS_MATCHED_DATA, DEMOGRAPHICS_TBL_SYN_ID)))
    unlink(OUTPUT_FILE$pd_case_v_control)
    
    ## store to synapse with the correct annotation
    file.entity<- synapser::File(OUTPUT_FILE$identity_confounding, parent = OUTPUT_SYN_ID)
    file.entity$annotations <- list(
        analysisType = "identity confounding",
        pipelineStep = "healthcode subsampling",
        userSubset = get("metadata")$user_group
    )
    synStore(
        file.entity, activity = Activity(
            'replicate matched RData',
            executed = GIT_URL,
            used = c(IDENTITY_CONFOUNDING_MATCHED_DATA, DEMOGRAPHICS_TBL_SYN_ID)))
    unlink(OUTPUT_FILE$identity_confounding)
    
    
    ## store to synapse with the correct annotation
    file.entity<- synapser::File(OUTPUT_FILE$N_of_1, parent = OUTPUT_SYN_ID)
    file.entity$annotations <- list(
        analysisType = "n of 1 analysis",
        pipelineStep = "healthcode subsampling",
        userSubset = get("metadata")$user_group
    )
    synStore(
        file.entity, activity = Activity(
            'replicate matched RData',
            executed = GIT_URL,
            used = c(TAP_TZ_DATA, 
                     REST_TZ_DATA,
                     VOICE_TZ_DATA,
                     WALK_TZ_DATA,
                     DEMOGRAPHICS_TBL_SYN_ID)))
    unlink(OUTPUT_FILE$pd_case_v_control)
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



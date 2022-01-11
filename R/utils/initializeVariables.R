#################################################################
#' Script to initalize variables and features used for the study
#' @author: Aryton Tediarjo
#' @author_email: aryton.tediarjo@sagebase.org
#################################################################
library(purrr)
library(synapser)
library(config)
library(githubr)

#' Reference Function for initializing features used for pipeline workflow
#' 
#' @return 
#' A list of each activity
get_features <- function(){
    tapFeatures <- c("meanTapInter", "medianTapInter", "iqrTapInter", "minTapInter","maxTapInter",      
                     "skewTapInter", "kurTapInter",  "sdTapInter", "madTapInter","cvTapInter",       
                     "rangeTapInter", "tkeoTapInter",  "ar1TapInter",  "ar2TapInter","fatigue10TapInter",
                     "fatigue25TapInter", "fatigue50TapInter", "meanDriftLeft","medianDriftLeft",  "iqrDriftLeft",     
                     "minDriftLeft", "maxDriftLeft", "skewDriftLeft",  "kurDriftLeft", "sdDriftLeft",      
                     "madDriftLeft", "cvDriftLeft",  "rangeDriftLeft", "meanDriftRight", "medianDriftRight", 
                     "iqrDriftRight", "minDriftRight", "maxDriftRight",  "skewDriftRight", "kurDriftRight",    
                     "sdDriftRight", "madDriftRight",  "cvDriftRight",  "rangeDriftRight", "numberTaps",       
                     "buttonNoneFreq")
    walkFeatures <- c('meanX', 'sdX', 'modeX', 'skewX', 'kurX', 'q1X', 'medianX', 'q3X', 
                      'iqrX', 'rangeX', 'acfX', 'zcrX', 'dfaX', 'cvX', 'tkeoX', 'F0X', 
                      'P0X', 'F0FX', 'P0FX', 'medianF0FX', 'sdF0FX', 'tlagX', 'meanY', 
                      'sdY', 'modeY', 'skewY', 'kurY', 'q1Y', 'medianY', 'q3Y', 'iqrY', 
                      'rangeY', 'acfY', 'zcrY', 'dfaY', 'cvY', 'tkeoY', 'F0Y', 'P0Y', 'F0FY', 
                      'P0FY', 'medianF0FY', 'sdF0FY', 'tlagY', 'meanZ', 'sdZ', 'modeZ', 'skewZ', 
                      'kurZ', 'q1Z', 'medianZ', 'q3Z', 'iqrZ', 'rangeZ', 'acfZ', 'zcrZ', 'dfaZ', 'cvZ', 
                      'tkeoZ', 'F0Z', 'P0Z', 'F0FZ', 'P0FZ', 'medianF0FZ', 'sdF0FZ', 'tlagZ', 'meanAA', 
                      'sdAA', 'modeAA', 'skewAA', 'kurAA', 'q1AA', 'medianAA', 'q3AA', 'iqrAA', 'rangeAA', 
                      'acfAA', 'zcrAA', 'dfaAA', 'cvAA', 'tkeoAA', 'F0AA', 'P0AA', 'F0FAA', 'P0FAA', 
                      'medianF0FAA', 'sdF0FAA', 'tlagAA', 'meanAJ', 'sdAJ', 'modeAJ', 'skewAJ', 'kurAJ', 'q1AJ', 
                      'medianAJ', 'q3AJ', 'iqrAJ', 'rangeAJ', 'acfAJ', 'zcrAJ', 'dfaAJ', 'cvAJ', 'tkeoAJ', 'F0AJ', 
                      'P0AJ', 'F0FAJ', 'P0FAJ', 'medianF0FAJ', 'sdF0FAJ', 'tlagAJ', 'corXY', 'corXZ', 'corYZ')
    restFeatures <- c('meanAA', 'sdAA', 'modeAA', 'skewAA', 'kurAA', 'q1AA', 'medianAA', 
                      'q3AA', 'iqrAA', 'rangeAA', 'acfAA', 'zcrAA', 'dfaAA', 'turningTime', 
                      'postpeak', 'postpower', 'alpha', 'dVol', 'ddVol')
    voiceFeatures <- c('Median_F0', 'Mean_Jitter', 'Median_Jitter', 'Mean_Shimmer', 'Median_Shimmer', 
                       'MFCC_Band_1', 'MFCC_Band_2', 'MFCC_Band_3', 'MFCC_Band_4', 'MFCC_Jitter_Band_1_Positive', 
                       'MFCC_Jitter_Band_2_Positive', 'MFCC_Jitter_Band_3_Positive', 'MFCC_Jitter_Band_4_Positive')
    return(list(tapping = tapFeatures,
                walking = walkFeatures,
                resting = restFeatures,
                voice = voiceFeatures))
}


#' Function for calling the ID of Synapse File View
#' given projectID and name of the generated file view
#' @return Synapse ID of the synapse table file view
get_file_view_ref <- function(){
    file_view <- map_chr(
        as.list(synGetChildren(get("output")$project_id)), function(x){
            if(x$name == get("output")$file_view_name){
                return(x$id)}else{return("")}}) %>% reduce(paste0)
    return(file_view)
}

#' Reference Function for calling the ID of Synapse File View
#' 
#' @return 
#' Synapse ID of file view
get_synapse_table_ref <- function(){
    SYN_ID_REF <- list()
    SYN_ID_REF$demo  <- get("synapse_tables")$demo
    SYN_ID_REF$walking  <- get("synapse_tables")$gait
    SYN_ID_REF$tapping   <- get("synapse_tables")$tap
    SYN_ID_REF$voice <- get("synapse_tables")$voice
    return(SYN_ID_REF)
}

#' Reference Function for calling the I/O of rawfeatures
#' 
#' @return 
#' List of synapse id of each activity and its output folders
get_raw_features_ref <- function(){
    file_view <- get_file_view_ref()
    SYN_ID_REF <- list()
    SYN_ID_REF$output_folder <- as.data.frame(synTableQuery(
        sprintf(
            "SELECT * FROM %s where pipelineStep = 'raw' AND type = 'folder'", 
            file_view)))$id
    RAW_FEATURE_FILES <- as.data.frame(synTableQuery(
        sprintf(
            "SELECT * FROM %s where parentId = '%s'", 
                file_view, SYN_ID_REF$output_folder)))
    SYN_ID_REF$tapping <- (RAW_FEATURE_FILES%>%
                               filter(task == "tapping" &
                                       pipelineStep == "raw"))$id
    SYN_ID_REF$walking <- (RAW_FEATURE_FILES%>%
                            filter(
                                task == "walking" & pipelineStep == "raw"))$id
    SYN_ID_REF$resting <- (RAW_FEATURE_FILES%>%
                               filter(
                                   task == "resting" & pipelineStep == "raw"))$id
    SYN_ID_REF$voice <- (RAW_FEATURE_FILES%>%
                             filter(
                                 task == "voice" & pipelineStep == "raw"))$id
    return(SYN_ID_REF)
}

#' Reference Function for calling the I/O of processed features
#' 
#' @return 
#' List of synapse id of each activity data and its output folders
get_processed_features_ref <- function(){
    file_view <- get_file_view_ref()
    SYN_ID_REF <- list()
    SYN_ID_REF$output_folder <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where pipelineStep = 'processed' AND type = 'folder'", 
                file_view)))$id
    PROCESSED_FEATURE_FILES <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where parentId = '%s'", 
                file_view, 
                SYN_ID_REF$output_folder)))
    SYN_ID_REF$tapping <- (PROCESSED_FEATURE_FILES %>%
                               filter(task == "tapping" &
                                          pipelineStep == "processed"))$id
    SYN_ID_REF$walking <- (PROCESSED_FEATURE_FILES %>%
                               filter(
                                   task == "walking" & pipelineStep ==  "processed"))$id
    SYN_ID_REF$resting <- (PROCESSED_FEATURE_FILES %>%
                               filter(
                                   task == "resting" & pipelineStep ==  "processed"))$id
    SYN_ID_REF$voice <- (PROCESSED_FEATURE_FILES %>%
                             filter(
                                 task == "voice" & pipelineStep ==  "processed"))$id
    return(SYN_ID_REF)
}

#' Reference Function for calling the I/O of intermediate data as a result after analysis
#' 
#' @return 
#' List of synapse id of each intermediate results and its output folders
get_intermediate_data_ref <- function(){
    file_view <- get_file_view_ref()
    SYN_ID_REF <- list()
    SYN_ID_REF$output_folder <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where pipelineStep = 'intermediary data' AND type = 'folder'", 
                file_view)))$id
    INTERMEDIARY_DATA <-  as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where parentId = '%s'", 
                file_view, 
                SYN_ID_REF$output_folder)))
    SYN_ID_REF$check_walking_adherence <- (INTERMEDIARY_DATA%>%
                                               filter(analysisType == "walking task adherence"))$id
    SYN_ID_REF$n_of_1_analysis <- (INTERMEDIARY_DATA%>%
                                       filter(analysisType == "n of 1 analysis"))$id
    SYN_ID_REF$confounder_corr <- (INTERMEDIARY_DATA%>%
                                       filter(analysisType == "demographics confounders" & 
                                                  analysisSubtype == "correlation test"))$id
    SYN_ID_REF$confounder_dcorr  <- (INTERMEDIARY_DATA%>%
                                        filter(analysisType == "demographics confounders" & 
                                                   analysisSubtype == "distance correlation test"))$id
    SYN_ID_REF$collapsed_pd_vs_nonpd <- (INTERMEDIARY_DATA%>%
                                              filter(analysisType == "case vs controls" & 
                                                         analysisSubtype == "collapsed measurements"))$id
    SYN_ID_REF$repeated_pd_vs_nonpd <- (INTERMEDIARY_DATA%>%
                                            filter(analysisType == "case vs controls" & 
                                                       analysisSubtype == "repeated measurements"))$id
    SYN_ID_REF$vperm_pd_vs_nonpd <- (INTERMEDIARY_DATA%>%
                                         filter(analysisType == "case vs controls" & 
                                                    analysisSubtype == "feature variability comparison"))$id
    SYN_ID_REF$identity_confounding <- (INTERMEDIARY_DATA%>%
                                         filter(analysisType == "identity confounding"))$id
    SYN_ID_REF$obj_pd_conf_score <- (INTERMEDIARY_DATA%>%
                                            filter(analysisType == "combined model"))$id
    return(SYN_ID_REF)
}

#' Reference Function for calling the I/O of matched healthcodes
#' 
#' @return 
#' List of synapse id of each matched healthcodes and its output folders
get_healthcode_ref <- function(){
    file_view <- get_file_view_ref()
    SYN_ID_REF <- list()
    SYN_ID_REF$output_folder <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where pipelineStep = 'healthcode subsampling' AND type = 'folder'", 
                file_view)))$id
    healthcode_df <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where parentId = '%s'", 
                file_view, 
                SYN_ID_REF$output_folder)))
    SYN_ID_REF$case_vs_controls <- (healthcode_df%>%
                                        dplyr::filter(
                                            analysisType == "case vs controls"))$id
    SYN_ID_REF$identity_confounding <- (healthcode_df%>%
                                        dplyr::filter(
                                            analysisType == "identity confounding"))$id
    SYN_ID_REF$n_of_one <- (healthcode_df %>%
                                        dplyr::filter(
                                            analysisType == "n of 1 analysis"))$id
    
    return(SYN_ID_REF)
}


#' Reference Function for calling the I/O of supplementary figures
#' 
#' @return 
#' List of synapse id of each supplementary figures and its output folders
get_figure_ref <- function(){
    file_view <- get_file_view_ref()
    SYN_ID_REF <- list()
    SYN_ID_REF$output_folder <- as.data.frame(synTableQuery(
        sprintf("SELECT * FROM %s where pipelineStep = 'figures' AND type = 'folder'", 
                file_view)))$id
    return(SYN_ID_REF)
}

get_obj_pd_ref <- function(){
    view_id <- "syn23545224"
    group <- config::get("metadata")$user_group
    SYN_ID_REF <- list()
    SYN_ID_REF$clinical <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'clinical' AND userSubset = '{group}'")))$id
    SYN_ID_REF$mapping <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'mapping' AND userSubset = '{group}'")))$id
    SYN_ID_REF$updrs <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'updrs' AND userSubset = '{group}'")))$id
    SYN_ID_REF$tapping <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'tapping' AND userSubset = '{group}'")))$id
    SYN_ID_REF$walking <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'walking' AND userSubset = '{group}'")))$id
    SYN_ID_REF$resting <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'resting' AND userSubset = '{group}'")))$id
    SYN_ID_REF$voice <- as.data.frame(synTableQuery(
        glue::glue(
            "SELECT * FROM {view_id} where pipelineStep = 'test data' AND task = 'voice' AND userSubset = '{group}'")))$id
    return(SYN_ID_REF)
}

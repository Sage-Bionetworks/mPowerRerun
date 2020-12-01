#############################################################
# Script to create project structure and File View in Synapse 
# for rerunning mPower Analysis
#############################################################
library(synapser)
library(config)
library(tidyverse)
library(githubr)
library(purrr)

#############################################################
# Configure Synapse Credentials, Github, and Config
#############################################################  
synLogin()
config::get()

#############################################################
# Helpers
#############################################################  

create.project.annotation <- function(){
  project.entity <- synGet(get("output")$project_id)
  project.annotations <- list(
    study = get("metadata")$study,
    diagnosis = c("Parkinson's Disease", "Control"),
    consortium = "mHealth",
    dataCollectionMethod = "active",
    deviceType ="handheld",
    sensorType = c("accelerometer", "gyroscope", "microphone"),
    devicePlatform = c("iOS"),
    deviceLocation = c("pocket", "flat surface"),
    reportedOutcome = c("medication report", "MDS-UPDRS"),
    digitalAssessmentCategory = c("motor coordination", "gait"),
    digitalAssessmentDetails = c("walking", "tapping", "resting", "phonation"),
    dhPortalIndex = FALSE,
    isDHProject = FALSE,
    dataAccessInstructions = NA,
    studyDescription = NA)
  synSetAnnotations(project.entity, annotations = project.annotations)
}


#' function for creating parent folder
create.parent.folder <- function(){
  folder <- Folder(get('output')$folder_name, 
                   parent = get('output')$project_id)
  folder <- synStore(folder)
  return(folder$properties$id)
}

#' function for subfolders of each analysis steps
create.pipeline.template <- function(parentId){
  ## create basic skeleton
  folder.mapping <- list(
    "intermediary data" = "Intermediary Data",
    "figures" = "Figures", 
    "raw" = "Raw Features",
    "processed" = "Processed Features",
    "healthcode subsampling" = "Filtered Healthcodes")
  
  store.folder <- purrr::map(names(folder.mapping),
             function(subtype){
               folder <- Folder(folder.mapping[[subtype]], 
                                parent = parentId)
               folder$annotations <- list(
                 study = get("metadata")$study,
                 pipelineStep = tolower(subtype),
                 userGroup = get("metadata")$user_group)
               folder <- synStore(folder)
               return(folder$properties$id)})
  return(parentId)
}


#' function for creating file view with annotations
create.file.view <- function(scopeId){
  EntityViewSchema(name= get('output')$file_view_name,
                  columns = c(
                             Column(name = "study",
                                    columnType = "STRING"),
                             Column(name = "deviceType",
                                    columnType = "STRING"),
                             Column(name = "sensorType",
                                    columnType = "STRING_LIST"),
                             Column(name = "digitalAssessmentCategory",
                                    columnType = "STRING_LIST"),
                             Column(name = "digitalAssessmentDetails",
                                    columnType = "STRING_LIST"),
                             Column(name = "dataAccessInstructions",
                                    columnType = "STRING"),
                             Column(name = "dataType",
                                    columnType = "STRING"),
                             Column(name = "dataSubtype",
                                    columnType = "STRING"),
                             Column(name = "analysisType",
                                    columnType = "STRING"),
                             Column(name = "analysisSubtype",
                                    columnType = "STRING"),
                             Column(name = "pipelineStep",
                                    columnType = "STRING"),
                             Column(name = "userSubset",
                                    columnType = "STRING")),
                           parent = get('output')$project_id,
                           add_default_columns=F,
                           scopes = scopeId,
                           includeEntityTypes=c(EntityViewType$FILE, 
                                                EntityViewType$FOLDER)) %>% synStore(.)
}

#############################################################
# Main Function
#############################################################  
main <- function(){
  create.project.annotation()
  create.parent.folder() %>% 
    create.pipeline.template(.) %>% 
    create.file.view(.)
}

main()




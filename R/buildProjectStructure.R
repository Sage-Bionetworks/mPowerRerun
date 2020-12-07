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
                 pipelineStep = tolower(subtype),
                 userSubset = get("metadata")$user_group)
               folder <- synStore(folder)
               return(folder$properties$id)})
  return(parentId)
}


#' function for creating file view with annotations
create.file.view <- function(scopeId){
  EntityViewSchema(name= get('output')$file_view_name,
                  columns = c(
                             Column(name = "task",
                                    columnType = "STRING",
                                    maximumSize = 10),
                             Column(name = "analysisType",
                                    columnType = "STRING",
                                    maximumSize = 30),
                             Column(name = "analysisSubtype",
                                    columnType = "STRING"),
                             Column(name = "pipelineStep",
                                    columnType = "STRING",
                                    maximumSize = 22),
                             Column(name = "userSubset",
                                    columnType = "STRING",
                                    maximumSize = 15),
                             Column(name = "consortium",
                                    columnType = "STRING",
                                    maximumSize = 10),
                             Column(name = "study",
                                    columnType = "STRING",
                                    maximumSize = 50),
                             Column(name = "studyOrProject",
                                    columnType = "STRING"),
                             Column(name = "sensorType",
                                    columnType = "STRING_LIST",
                                    maximumSize = 20),
                             Column(name = "deviceType",
                                    columnType = "STRING_LIST",
                                    maximumSize = 10),
                             Column(name = "devicePlatform",
                                    columnType = "STRING_LIST",
                                    maximumSize = 20),
                             Column(name = "dataCollectionMethod",
                                    columnType = "STRING_LIST",
                                    maximumSize = 15),
                             Column(name = "deviceLocation",
                                    columnType = "STRING_LIST",
                                    maximumSize = 15),
                             Column(name = "diagnosis",
                                    columnType = "STRING_LIST",
                                    maximumSize = 20),
                             Column(name = "reportedOutcome",
                                    columnType = "STRING_LIST",
                                    maximumSize = 20),
                             Column(name = "digitalAssessmentCategory",
                                    columnType = "STRING_LIST",
                                    maximumSize = 15),
                             Column(name = "digitalAssessmentDetails",
                                    columnType = "STRING_LIST",
                                    maximumSize = 10),
                             Column(name = "dataType",
                                    columnType = "STRING",
                                    maximumSize = 15),
                             Column(name = "dataSubtype",
                                    columnType = "STRING",
                                    maximumSize = 15),
                             Column(name = "dhPortalIndex",
                                    columnType = "BOOLEAN"),
                             Column(name = "dataDescriptionLocation",
                                    columnType = "STRING",
                                    maximumSize = 30),
                             Column(name = "dataAccessInstructions",
                                    columnType = "STRING",
                                    maximumSize = 30)),
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
  create.parent.folder() %>% 
    create.pipeline.template(.) %>% 
    create.file.view(.)
}

main()




library(synapser)
library(config)
library(tidyverse)
library(data.table)
library(bridgeclient)
library(githubr)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

# login to synapse
synapser::synLogin()

# fetch bridge credentials
# create .bridge_creds to pass credentials
bridgeclient::bridge_login(
    study = "parkinson",
    credentials_file = "R/InternalUse/.bridge_creds")

# get config files
config::get()

# authenticate github
setGithubToken(readLines(
    file.path(path.expand("~"), "git_token.txt")))
SCRIPT_NAME <- "filter_public_sharing_scope_objPD.R"
GIT_URL <- getPermlink(
    getRepo(get("git")$repo,
            ref="branch", 
            refName="filterObjPD"), 
    repositoryPath = file.path(
        'R/InternalUse', SCRIPT_NAME))

# reference
OUTPUT_PARENT_ID <- "syn26715583"
MAPPING_IDENTIFIER <- "syn8533708"
CLINICAL_DATA <- "syn12178233"
UPDRS_SCORES <- "syn8232164"
FEATURES <- list()
FEATURES$tap <- "syn8075029"
FEATURES$voice <- "syn8262306"
FEATURES$walk <- "syn8089484"
FEATURES$rest <- "syn8225448"

store_to_synapse <- function(data, filename, parent_id, ...){
    data %>%
        readr::write_tsv(filename)
    file <- synapser::File(filename, parent = parent_id)
    store <- synapser::synStore(file, ...)
}

filter_identifier <- function(data){
    filter_hc <- purrr::map_dfr(data$healthCode, function(hc){
        content <- bridgeclient::get_participant(health_code = hc)
        tibble::tibble(
            healthCode = content$healthCode,
            sharingScope = content$sharingScope)}) %>%
        dplyr::filter(sharingScope == "all_qualified_researchers") %>%
        dplyr::select(healthCode)
    data %>%
        dplyr::inner_join(filter_hc, by = c("healthCode"))
}

identifier <- synGet(MAPPING_IDENTIFIER)$path %>% 
    fread() %>%
    tibble::as_tibble() %>%
    dplyr::select(healthCode = healthcode, everything())

clinical_data <- synGet(CLINICAL_DATA)$path %>% 
    fread() %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(identifier)

updrs_scores <- synGet(UPDRS_SCORES, version = 6)$path %>% 
    fread() %>%
    tibble::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        id = glue::glue("ROC", sprintf("%03d", record))) %>%
    dplyr::inner_join(identifier) %>%
    dplyr::ungroup()

feature_list <- purrr::map(FEATURES, function(id){
    data <- synapser::synGet(id)$path %>% fread()
    data %>%
        tibble::as_tibble() %>%
        dplyr::inner_join(identifier)})

identifier %>% 
    store_to_synapse(filename = "objPD_mapping.tsv",
                     parent_id = OUTPUT_PARENT_ID,
                     activityName = "get subset from public sharing scope",
                     used = MAPPING_IDENTIFIER,
                     executed = GIT_URL)

updrs_scores %>% 
    store_to_synapse(filename = "objPD_updrs_scores.tsv",
                     parent_id = OUTPUT_PARENT_ID,
                     activityName = "get subset from public sharing scope",
                     used = UPDRS_SCORES,
                     executed = GIT_URL)

clinical_data %>% 
    store_to_synapse(filename = "objPD_clinical_data.tsv",
                     parent_id = OUTPUT_PARENT_ID,
                     activityName = "filter clinical files",
                     used = CLINICAL_DATA,
                     executed = GIT_URL)

purrr::map(names(feature_list), function(activity){
    output_filename <- glue::glue("objPD_", 
                                  activity, 
                                  "_features.tsv")
    feature_list[[activity]] %>% 
        store_to_synapse(filename = output_filename,
                         parent_id = OUTPUT_PARENT_ID,
                         activityName = "filter activities",
                         used = FEATURES[[activity]],
                         executed = GIT_URL)
})

    



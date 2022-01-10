library(synapser)
library(config)
library(tidyverse)
library(data.table)
library(bridgeclient)
library(githubr)
source("R/utils/initializeVariables.R")
source("R/utils/projectUtils.R")

synapser::synLogin()
bridgeclient::bridge_login(
    study = "parkinson",
    credentials_file = "R/InternalUse/.bridge_creds")

config::get()
setGithubToken(readLines(
    file.path(path.expand("~"), "git_token.txt")))
SCRIPT_NAME <- "filter_public_sharing_scope_objPD.R"
GIT_URL <- getPermlink(
    getRepo(get("git")$repo,
            ref="branch", 
            refName=get("git")$branch), 
    repositoryPath = file.path(
        'R/InternalUse', SCRIPT_NAME))

OUTPUT_PARENT_ID <- "syn26715149"
MAPPING_IDENTIFIER <- "syn8533708"
CLINICAL_DATA <- "syn12178233"
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
    dplyr::select(healthCode = healthcode, everything()) %>%
    filter_identifier()

clinical_data <- synGet(CLINICAL_DATA)$path %>% 
    fread() %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(identifier)

feature_list <- purrr::map(FEATURES, function(id){
    data <- synapser::synGet(id)$path %>% fread()
    data %>%
        tibble::as_tibble() %>%
        dplyr::inner_join(identifier)})

identifier %>% 
    store_to_synapse(filename = "public_shared_objPD_mapping.tsv",
                     parent_id = OUTPUT_PARENT_ID,
                     activityName = "get subset from public sharing scope",
                     used = MAPPING_IDENTIFIER)
metadata %>% 
    store_to_synapse(filename = "public_shared_objPD_clinical_data.tsv",
                     parent_id = OUTPUT_PARENT_ID,
                     activityName = "filter clinical files",
                     used = CLINICAL_DATA)

purrr::map(names(feature_list), function(activity){
    output_filename <- glue::glue("public_shared_objPD_", 
                                  activity, 
                                  "_features.tsv")
    feature_list[[activity]] %>% 
        store_to_synapse(filename = output_filename,
                         parent_id = OUTPUT_PARENT_ID,
                         activityName = "filter activities",
                         used = FEATURES[[activity]])
})

    



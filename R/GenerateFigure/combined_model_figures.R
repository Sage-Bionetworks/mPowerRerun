#############################################################
#' This script is used to generate main text figure 3 
#' of the mPower Paper. It generates plot of
#' combined model vs chosen metrics (UPDRS, SE-ADL, Hoehn-Yahr)
#' author: aryton.tediarjo@sagebase.org, larsson.omberg@sagebase.org
############################################################
library(synapser)
library(config)
library(tidyverse)
library(dplyr)
library(jsonlite)
library(githubr)
library(ggpval)
library(ggExtra)
library(ggpubr)
library(rstatix)
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
## Global Variables
#######################################################

## SYNAPSE ID REFERENCE FOR OBJECTIVE PD
CLINICAL_FILE <- 'syn8232164'
SAMPLE_MAP <- 'syn8533708'

## GLOBAL VARIABLES
N_TESTS <- 3
FIGURE_NAME <- paste0("mPower_",
                      gsub(" ", "_", get("metadata")$user_group), 
                      "_main_text_figure_3",".png")
SYN_ID_REF <- list(figures = get_figure_ref(),
                   intermediate = get_intermediate_data_ref())
SCRIPT_NAME <-  "combined_model_figures.R"
GIT_URL <- getPermlink(getRepo(get("git")$repo,
                               ref="branch", 
                               refName=get("git")$branch), 
                       repositoryPath = file.path("R/GenerateFigure", SCRIPT_NAME))
FIGURE_OUTPUT_SYN_ID <- SYN_ID_REF$figures$output_folder
MODEL_SCORE <- SYN_ID_REF$intermediate$obj_pd_conf_score
ANNOTATIONS <- list(
    analysisType = "combined model",
    analysisSubtype = "confidence score on objective PD users",
    userSubset = tolower(get("metadata")$user_group), 
    pipelineStep = "figures")


######################################################
## Helper
#######################################################
#' Function to get all required data
#' 
#' @return a named list containing dataframes of scores from AT_HOME_PD
get_objectivePD_data <- function(){
    model.df <- read.csv(synGet(MODEL_SCORE)$path, sep = "\t") %>% 
        tibble::as_tibble(.) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate_all(~ifelse(is.na(.), median(c(tapping, voice, walk, rest), na.rm = TRUE), .)) %>%
        dplyr::mutate(combined.model = median(c_across(where(is.numeric))))
    sample.map <- read.csv(synGet(SAMPLE_MAP)$path) %>% as_tibble(.)
    clinical.df <- read.csv(synGet(CLINICAL_FILE, version = 6)$path) %>%
        as_tibble(.) %>%
        rowwise() %>%
        dplyr::mutate(
            externalId = glue::glue("ROC", sprintf("%03d", record)),
            time = dplyr::case_when(
                descrip == 'Baseline (enrollment)'~ 0,
                descrip == 'Month 3' ~ 90, descrip == 'Month 6' ~ 180, TRUE ~ NaN)) %>%
        dplyr::ungroup()
    pd.mapping <- clinical.df %>% 
        dplyr::select(externalId, Do.you.have.Parkinson.disease.) %>% 
        dplyr::mutate_if(is.character, list(~na_if(.,""))) %>% 
        dplyr::rename("PD" = "Do.you.have.Parkinson.disease.") %>%
        drop_na()
    scores <- clinical.df %>% 
        dplyr::group_by(externalId) %>% 
        dplyr::summarise_at(c("UPDRS.Total.Score",  
                              "Schwab.and.England.Activities.of.Daily.Living.Scale.Score"), 
                            .funs = c(mean, sd), na.rm = TRUE) %>%
        dplyr::rename('updrs.mean' = "UPDRS.Total.Score_fn1", 
                      'updrs.std' = "UPDRS.Total.Score_fn2",
                      'se.adl.mean' = "Schwab.and.England.Activities.of.Daily.Living.Scale.Score_fn1",
                      'se.adl.std' = "Schwab.and.England.Activities.of.Daily.Living.Scale.Score_fn2") %>%
        dplyr::inner_join(pd.mapping, on = c("externalId")) %>% 
        dplyr::inner_join(sample.map, by = c("externalId" = "id")) %>%
        dplyr::inner_join(model.df, by = c("healthcode" = "healthCode")) %>% 
        dplyr::select(healthCode = healthcode, everything())
    hoehn.yahr <- clinical.df %>% 
        dplyr::select(externalId, Estimated.Hoehn.and.Yahr.stage) %>%
        dplyr::select(externalId, HY = Estimated.Hoehn.and.Yahr.stage) %>% 
        dplyr::mutate(HY = as.factor(HY)) %>%
        tidyr::drop_na(c(HY)) %>%
        dplyr::inner_join(scores %>% dplyr::select(externalId, healthCode, combined.model), by = c("externalId"))
    return(list(scores = scores, hoehn.yahr = hoehn.yahr))
}

#' Function to get stat metrics for annotating graphs with p-value corrections
#' 
#' @param scores: Dataframe containing metrics of combined model score with
#' UPDRS score and SE-ADL score
#' @param hoehn.yahr: Dataframe containing hoehn.yahr score to combined model score
#' @return a named list containing dataframes of metrics with adjusted scores
get_stat_test_metrics <- function(scores, hoehn.yahr){
    metrics.list <- list()
    metrics.list$updrs <- scores %>% 
        tidyr::drop_na(c(updrs.mean, combined.model)) %>% 
        dplyr::do(cor.test(.$updrs.mean, .$combined.model) %>% tidy) %>% 
        dplyr::mutate(p.value = p.value * N_TESTS) %>%
        dplyr::select(estimate, p.value)
    metrics.list$se.adl <- scores %>% 
        tidyr::drop_na(c(se.adl.mean, combined.model)) %>% 
        dplyr::do(cor.test(.$se.adl.mean, .$combined.model) %>% tidy) %>% 
        dplyr::mutate(p.value = p.value * N_TESTS) %>%
        dplyr::select(estimate, p.value)
    metrics.list$HY <- hoehn.yahr %>% 
        rstatix::t_test(combined.model ~ HY, var.equal = TRUE, 
                        p.adjust.method = "bonferroni", 
                        comparisons = list(c(0,1), c(0,2), c(0,3))) 
    return(metrics.list)
}    

#' Function to generate plot for main text figure 3, 
#' a combined model performance (median) 
#' faceted to 3 different PD metrics (SE-ADL, HoehnYahr, UPDRS)
#' @param scores: Dataframe containing metrics of combined model s
#' core with UPDRS score and SE-ADL score
#' @param hoehn.yahr: Dataframe containing hoehn.yahr 
#' score to combined model score
#' @param metrics.list: Named list containing statistical 
#' tests metrics for each plot
#' @return 
generate_plot <- function(scores, hoehn.yahr, metrics.list){
    #updrs.plot
    updrs.plot <- (ggscatter(
        scores %>% 
            drop_na(c(updrs.mean, combined.model)), 
        x = "combined.model", 
        y = "updrs.mean", 
        add = "reg.line", 
        add.params = list(color = "red", fill = "lightgray"), 
        conf.int = TRUE, fullrange = T) + 
            annotate("label", x = 0.3, y = 85, size = 2.5,
                     label = glue::glue("R = ", sprintf(metrics.list$updrs$estimate[[1]], fmt = '%#.2f')," ", 
                                        "P = ", formatC(metrics.list$updrs$p.value[[1]], format = "e", digits = 2))) + 
            geom_errorbar(aes(ymin = updrs.mean - updrs.std, ymax = updrs.mean + updrs.std)) + 
            scale_x_continuous(name = "Combined Model", expand=c(0,0), limits=c(0,1.05)) + 
            scale_y_continuous(name = "UPDRS Total Score", expand=c(0,0), limits=c(0,120)) + 
            coord_cartesian(xlim = c(0,1.05), ylim = c(-1.5,101)) + labs(tag = "a).") +
            theme_classic()) %>% 
        ggExtra::ggMarginal(., type = "densigram", 
                            xparams = list(fill = "white"), 
                            yparams = list(fill = "white"))
    
    ## SE-ADL plot
    se.adl.plot <- (ggscatter(
        scores %>% 
            drop_na(c(se.adl.mean, combined.model)),
        x = "combined.model", 
        y = "se.adl.mean", 
        add = "reg.line", 
        add.params = list(color = "red", fill = "lightgray"), 
        conf.int = TRUE, fullrange = T) + 
            annotate("label", x = 0.3, y = 10, size = 2.5,
                     label = glue::glue("R = ", sprintf(metrics.list$se.adl$estimate[[1]], fmt = '%#.2f')," ", 
                                        "P = ", formatC(metrics.list$se.adl$p.value[[1]], format = "e", digits = 2))) + 
            geom_errorbar(aes(ymin = se.adl.mean - se.adl.std, ymax = se.adl.mean + se.adl.std)) + 
            scale_x_continuous(name = "Combined Model", expand=c(0,0), limits=c(0,1.1)) + 
            scale_y_continuous(name = "SE-ADL", expand=c(0,0), limits=c(0,120)) + 
            coord_cartesian(xlim = c(0,1.1), ylim = c(-1.5,101)) + labs(tag = "b).") +
            theme_classic()) %>% 
        ggExtra::ggMarginal(., type = "densigram", 
                            xparams = list(fill = "white"), 
                            yparams = list(fill = "white"))
    
    ## create H & Y plot
    hy.boxplot <- ggplot(hoehn.yahr %>% 
                             tidyr::drop_na(), aes(x = HY, y = combined.model)) + 
        stat_boxplot(geom = "errorbar", width = 0.5) +
        geom_boxplot() + 
        scale_x_discrete(name = "Hoehn and Yahr Stage") + 
        scale_y_continuous(name = "Combined Model", expand=c(0,0), limits=c(0.2,1.05)) + labs(tag = "c).") +
        theme_classic() +
        stat_pvalue_manual(metrics.list$HY %>% 
                               add_xy_position %>% 
                               dplyr::mutate(annot = paste("P = ", p.adj)), 
                           label = "annot", size = 2.5)
    res <- gridExtra::grid.arrange(updrs.plot, se.adl.plot, hy.boxplot, nrow = 1)
    png(FIGURE_NAME,width = 2500, height = 800, units = "px", res = 200) 
    grid::grid.draw(res) 
    dev.off()
}

main <- function(){
    data <- get_objectivePD_data()
    metrics.list <- get_stat_test_metrics(data$scores, data$hoehn.yahr)
    plot <- generate_plot(data$scores, data$hoehn.yahr, metrics.list)
    f <- synapser::File(FIGURE_NAME, parent=FIGURE_OUTPUT_SYN_ID)
    f$annotations <- ANNOTATIONS
    activity <- synapser::Activity(
        "generate main text figure 3",
        used = c(CLINICAL_FILE, SAMPLE_MAP, MODEL_SCORE),
        executed = GIT_URL)
    synapser::synStore(f, activity = activity)
    unlink(FIGURE_NAME)
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




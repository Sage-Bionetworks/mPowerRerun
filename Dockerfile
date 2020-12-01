# build base image
FROM rocker/tidyverse:4.0.0

# install preliminary requirements
RUN apt-get update -y\
    && apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev libglu1-mesa-dev\
    && apt-get install -y curl libcurl4-openssl-dev\
    && apt-get install -y git

## install R environment
RUN R -e "install.packages('remotes')"\
    && R -e "install.packages('rgl', dependencies = TRUE)"\
    && R -e "remotes::install_github('Sage-Bionetworks/mpowertools')"\
    && R -e "install.packages('PythonEmbedInR', repos=c('http://cran.fhcrc.org', 'http://ran.synapse.org'))"\
    && R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
    && R -e "install.packages('BiocManager')"\
    && R -e "BiocManager::install('rhdf5')"\
    && R -e "install.packages('pROC')"\
    && R -e "install.packages('ppcor')"\
    && R -e "install.packages('randomForest')"\
    && R -e "install.packages('glmnet', dependencies=TRUE)"\
    && R -e "install.packages('lubridate')"\
    && R -e "install.packages('energy')"\
    && R -e "install.packages('lmtest')"\
    && R -e "install.packages('forecast')"\
    && R -e "install.packages('sandwich')"\
    && R -e "install.packages('gplots')"\
    && R -e "install.packages('dplyr')"\
    && R -e "install.packages('doMC')"\
    && R -e "devtools::install_github('th1vairam/CovariateAnalysis@dev')"\
    && R -e "devtools::install_github('brian-bot/githubr')"\
    && R -e "install.packages('stringr')"\
    && R -e "install.packages('plyr')"\
    && R -e "install.packages('MatchIt')"\
    && R -e "install.packages('config')"\
    && R -e "install.packages('install.load')"\
    && R -e "install.packages('survminer')"\
    && R -e "install.packages('ggpubr')"\
    && R -e "install.packages('ggpval')"\
    && R -e "install.packages('ggExtra')"\
    && R -e "install.packages('rstatix')"\
    && R -e "install.packages('relaimpo')"

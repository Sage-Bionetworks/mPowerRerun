# build base image
FROM rocker/tidyverse:4.0.0

# install preliminary requirements
RUN apt-get update -y\
    && apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev libglu1-mesa-dev\
    && apt-get install -y curl libcurl4-openssl-dev\
    && apt-get install -y git
    
## run git cloning
RUN git clone --branch renv https://github.com/Sage-Bionetworks/mPowerRerun /root/mPowerRerun

## change work dir
WORKDIR /root/mPowerRerun

## install R environment using lockfile
ENV RENV_VERSION 0.13.2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e "renv::init(bare = TRUE)"
RUN R -e "renv::restore()"
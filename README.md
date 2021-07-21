# Sage Bionetworks - mPowerRerun

<img alt="GitHub pull requests" src="https://img.shields.io/github/issues-pr/Sage-Bionetworks/mPowerRerun">  <img alt="GitHub issues" src="https://img.shields.io/github/issues/Sage-Bionetworks/mPowerRerun">

## About
This code repository contains streamlined approach in rerunning mPower Nature Biotech's Publication named **Remote Smartphone Monitoring of Parkinson Disease and Individual Response to Therapy**. This repository streamlines all the analysis results from the paper starting from our data warehouse to publication results.

**Analysis being done in the Data Pipeline:**
1. PD Case vs Retention Analysis
2. Identity Confounding on Repeated Measures
3. PD Case vs Controls Analysis (Ridge Regression and Random Forest)
4. Variability Comparisons on Extracted Features based of Random Forest Model
5. N of 1 Analysis for Medication vs Time of Day (Arima, Newey-West) and Feature Relative Importances
6. Assessment on demographics confounders based of correlation and distance correlation
7. Random Forest Combined Model Performances to Standardized PD Metrics (UPDRS, SE-ADL, Hoehn Yahr)

Wiki showcasing the results and guide for getting figure results from the analysis [link to wiki](https://www.synapse.org/#!Synapse:syn23277418/wiki/606593)

## How-to-Run:

### 1) Credentials Requirements

- [Synapse account](https://docs.synapse.org/articles/getting_started.html) 
- [Github Personal Access Token](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token)

### 2) Clone this Github Repo
```
git clone https://github.com/Sage-Bionetworks/mPowerRerun
```

### 3) Setting up Environment & Dependencies
You have two options for setting up your environment, either through your own local machine or a Docker Container.
We use R `renv` package to manage all the dependencies and its versioning (**Note: Environment Requires R Version 4.0 and above**)

If restoring environment fails, our Docker documentation provide a Docker Image that can restore the libraries and its corresponding system dependencies.

#### Installing dependencies

##### Install `renv` library
```R
install.packages("renv") 
```
##### Restore environment using `renv`
```R
renv::init(bare = T) #create renv empty environment
renv::restore() #restore all libraries
```

#### (Optional) Docker
- [Install Docker.](https://docs.docker.com/v17.12/install/#supported-platforms)

##### Create Docker Image & Run Container
This Docker container is built on top of  [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse/) producing a debian stable work environment.

- Create Docker Image & Container:

##### Build Docker Image
```bash
docker build -t mpower_rerun . 
```
##### Run Docker Container
[Tutorial on many ways to run a Rstudio Docker container]("https://hub.docker.com/r/rocker/rstudio)


### 4) Set up config.yml and .Renviron

The data workflow will make use of config.yml to preserve all the required data input for both data 6-months-study and the public release. Before running the data workflow, we will require your Github Personal Access Token from Step 1 and the Synapse Project ID to store where your results will be stored. 

`.RProfile` will be used to control which `config.yml` option you want to use (default is set for public release).

#### i) Set up Git
```R
git:
  path: <PATH_TO_GIT_PERSONAL_ACCESS_TOKEN>
  repo: "Sage-Bionetworks/mPowerRerun"
  branch: "main"
```

#### ii) Set project Output
```R
 output:
    project_id: <SYNAPSE_PROJECT_ID>
    folder_name: "mPower Rerun Results - Public"
    file_view_name: "mPower Rerun Results - Public - File View"
```

#### iii) Configure .Renviron
a. Reproducing Public Release Results:
```R
R_CONFIG_ACTIVE = "default"
```
b. Reproducing 6 Months Study Results:
```R
R_CONFIG_ACTIVE = "6_months_study"
```

### 5) Running the Data Pipeline
To run the data pipeline, `Makefile` is used to maintaining the data workflow.

a. Reproducing 6 Months Study:
```R
make regenerate_paper
```
b. Reproducing Public Release:
```R
make all
```

If you are interested in doing part of the pipeline, refer to `Makefile` to run each of the steps individually


## Misc. Info:
#### a. Serialized Model
Serialized model of our end **Random Forest** (trained on sensor features only) into a folder called serializedModel/ in .RDS serialized file during the [objectivePD cohort prediction](https://github.com/arytontediarjo/mPowerRerun/blob/master/R/Analyses/trainOnMPower_predictObjPD.R).

#### b. Debugging & Logging
The pipelne process will be tracked by log files. 
- pipeline.log will track timestamps of each code execution.
- error.log will show which script is having an error.

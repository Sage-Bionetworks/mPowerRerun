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
You have two options for setting up your environment, either through your own local machine or a Docker Container

#### a. RStudio
We use R `renv` package to manage all the dependencies and its versioning, thus you would need `renv` library installed to your RStudio.

##### Installing dependencies
```R
install.packages("renv")
renv::init(bare = T)
renv::restore()
```

#### b. Docker
- [Install Docker.](https://docs.docker.com/v17.12/install/#supported-platforms)

#### Create Docker Image & Run Container
This Docker container is built on top of  [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse/) producing a debian stable work environment.

- Create Docker Image & Container:
```bash
docker build -t <IMAGE_NAME> . 
docker run -d <IMAGE_NAME> -p 8787:8787 -e PASSWORD="sage" <IMAGE_NAME>
```
- `-d` flag allows you the retain use of the terminal while the Docker image is running 
- `-p` specifies port of choice
- `-it` for interactive mode in the Docker container

### 4) Set up Config
Once all environment is set up, you will need to set your Git Personal Access Token and desired Synapse project ID for storing all the results of the analysis

#### i) Set up Git
```
git:
  path: <PATH_TO_GIT_PERSONAL_ACCESS_TOKEN>
  repo: "Sage-Bionetworks/mPowerRerun"
  branch: "main"
```

#### ii) Set project Output
```
 output:
    project_id: <SYNAPSE_PROJECT_ID>
    folder_name: "mPower Rerun Results - Public"
    file_view_name: "mPower Rerun Results - Public - File View"
```

#### iii) Configure .RProfile
a. Reproducing 6 Months Study:
```
Sys.setenv(R_CONFIG_ACTIVE = "default")
```
b. Reproducing Public Release Results:
```
Sys.setenv(R_CONFIG_ACTIVE = "public_release")
```

### 5) Running the Data Pipeline
To run the data pipeline, `Makefile` is used to maintaining the data workflow.

a. Reproducing 6 Months Study:
```
make regenerate_paper
```
b. Reproducing Public Release:
```
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

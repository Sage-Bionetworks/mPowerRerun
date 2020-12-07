# Sage Bionetworks - mPowerRerun

Author: Elias Chaibub Neto, Larsson Omberg, Aryton Tediarjo

CP: aryton.tediarjo@sagebase.org

## Introduction
This code repository contains streamlined approach in rerunning Sage Bionetworks Nature Biotech's Research Study named **"MPower - A Smartphone approach to Remotely Monitor PD and Individual Responses to Therapy"**. This repo will act as a pipeline for extracting data from Synapse into the intermediate data (analysis metrics, machine learning performance) that is used for the figure deliverables.

**Analysis being done in this Git Repository:**
1. PD Case vs Control Survival Analysis
2. Identity Confounding on Repeated Measures
3. PD Case vs Controls Analysis (Ridge Regression and Random Forest)
4. Variability Comparisons on Extracted Features based of Random Forest Model
5. N of 1 Analysis for Medication vs Time of Day (Arima, Newey-West) and Feature Relative Importances
6. Assessment on demographics confounders based of correlation and distance correlation
7. Random Forest Combined Model Performances to Standardized PD Metrics (UPDRS, SE-ADL, Hoehn Yahr)

We also have a wiki showcasing the results and guide for getting figure results from the analysis [link to wiki](https://www.synapse.org/#!Synapse:syn22151120/wiki/604781)

## Environment

### 1.) Credentials Requirements

- [Synapse account](https://docs.synapse.org/articles/getting_started.html) 
- [Github Personal Access Token](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token)

### 2.) Using Docker (Suggested)

#### Reference to Docker 
[Install Docker.](https://docs.docker.com/v17.12/install/#supported-platforms)

#### Create Docker Image & Run Container
This Docker container is built on top of  [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse/) producing a debian stable work environment.

- Create Docker Image:
```bash
docker build -t <IMAGE_NAME> . 
```
- Interactively in Bash:
```bash
docker run -it <IMAGE_NAME> /bin/bash
```
- `-d` flag allows you the retain use of the terminal while the Docker image is running 
- `-p` specifies port of choice
- `-it` for interactive mode in the Docker container

## Pipeline Steps

Once the environment is all set up, here are some quick steps you can take to run the project. 

1. [Create a new project in Synapse](https://docs.synapse.org/articles/making_a_project.html)
2. Create a txt file containing your Git Token Credentials (can save it anywhere and you can point it using the config file), set your Git repository path
```yml
git:
    path: "<path_to_git_token>/git_token.txt" #your path to git token
    repo: "<path_to_git_repo>/mPowerRerun" #your cloned mPowerRerun Github Repo
    branch: "main"
```
3. Set metadata that will be used to set the Synapse Annotations of the data
```yml
metadata:
    study: 'mPower' # the name of the study
    user_group: 'public data' # this will be used for naming convention and annotation
```
4. Set your output information
```yml
output:
    project_id: 'syn23277418' # refer to your desired output project id
    folder_name: "mPower Rerun Results" # the name of the output folder of your analysis results
    file_view_name: 'mPower Rerun File View' # the synapse file view used to store the data into Synapse Tables (SQL format)
```
5. Afterwards, the config file will contain information regarding the synapse tables, where you can freely change the Synapse id of each Synapse tables (respective activities).
```yml
synapse_tables:
    demo: "syn7222419"
    gait: "syn7222425" #walk and balance test 
    tap: "syn7222423"
    voice: "syn7222424"
additional:
    voice_features: 'syn22041873' #generated voice features from matlab
```

Note: Tap, Walk and Rest activities are fully reproducible, voice feature extraction will require Matlab, so we provided featurized dataset in Synapse. 

Once configurations are made, this project will be encapsulated into the usage of GNU Makefile, thus running `make all` to reproduce custom data or `make regenerate_paper` to reproduce the exact publications results in the project directory with your bash/terminal will streamline the whole process. You are also able to run it per stage of analysis (refer to the Makefile). 

## Miscellaneous
#### a. Serialized Model
We are storing the serialized model of our end **Random Forest** (trained on sensor features only) into a folder called serializedModel/ in .RDS serialized file during the [objectivePD cohort prediction](https://github.com/arytontediarjo/mPowerRerun/blob/master/R/Analyses/trainOnMPower_predictObjPD.R). This file can be used for your analytical purposes or making predictions based on the predefined sensor features of each activity.

#### b. Debugging & Logging
The pipelne process will be tracked by a logger; pipeline.log will track timestamps of each code execution, error.log will show which script is having an error.

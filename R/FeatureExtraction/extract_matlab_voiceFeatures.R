# extract_voiceFeatures <- function(HC){
# Wrapper code to extract voice features using Max Little's code in MATLAB

## It is assumed your working directory is where this file is

## REQUIREMENTS
#     libav-tools package in ubuntu for conversion of m4a to wav files
#     MATALAB with parallel computing and signal processing toolbox

## Clear R console screen output
cat("\014")

## Set input parameters
WORKING_DIR = '~/Documents/mpowervoice/mPowerAnalysis/voice_module'
SOURCE_TBL = 'syn7222424'; # Exported voice table id from internal project
FEATURES_PARENT_ID = 'syn7231638' # Destination directory in synapse
FEATURES_ID = 'syn11335364' # Destination file id where features are to be stored

## Set working and library directories
setwd(WORKING_DIR)

## Load libraries
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapser)
library(githubr) ## Needs dev branch

synapser::synLogin()

## Clear temp files in folder
system('rm temp*')

# ## Get commits from github
thisFileName <- 'extract_matlab_voiceFeatures.R'

thisRepo <- getRepo(repository = "itismeghasyam/mPowerAnalysis", ref="branch", refName='voice_module_dev') # Validate your github token before fetching repo
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('voice_module/',thisFileName))

activityName = 'Extract voice features';
activityDescription = "Extract voice features from audio m4a files using Max Little's features";

## Get all voice records from the source table
voice.tbl.syn = synTableQuery(paste0("SELECT * FROM ",SOURCE_TBL))
#," WHERE healthCode IN ('",HC,"')"))
voice.tbl <- voice.tbl.syn$asDataFrame()

## Get subset of healthCodes we need to download the datafrom
# Doing this because my local computer does not have space to download all the voice files
# and then convert them into wav files
# n = total number of records per healthCode
# classes: 
# 1. n>=400 [done]
# 2. n>=300 and n<400 [done]
# 3. n>=200 and n<300 [done]
# 4. n>=100 and n<200 [done]
# 5. n>=50 and n<100 [done]
# 6. n>=20 and n<50 [done]
# 7. n>=8 and n<20 [done]
# 8. n>=4 and n<8 [done]
# 9. n>=2 and n<4 [done]
# 10. n==1 [done]
hc.subset <- voice.tbl %>% 
  dplyr::select(recordId, healthCode) %>% 
  dplyr::group_by(healthCode) %>% 
  dplyr::count() %>% 
  unique() %>% 
  dplyr::filter(n>=2, n<4)

HC <- hc.subset[['healthCode']] %>% 
  unlist() %>% 
  paste0(collapse = '\',\'')

## Get subset class voice records from the source table
voice.tbl.syn = synTableQuery(paste0("SELECT * FROM ",SOURCE_TBL," WHERE healthCode IN ('",HC,"')"))
voice.tbl <- voice.tbl.syn$asDataFrame()

## Download all the voice audio and countdown files and attach file locations as columns
columnsToDownload = c("audio_audio.m4a",
                      "audio_countdown.m4a") 

voice.json.loc = lapply(columnsToDownload, function(col.name){
  tbl.files = synapser::synDownloadTableColumns(voice.tbl.syn, col.name) %>%
    lapply(function(x) data.frame(V1 = x)) %>% 
    data.table::rbindlist(idcol = col.name) %>% 
    plyr::rename(c('V1' = gsub('.m4a','.fileLocation', col.name)))
})


# tbl.file.ids = synapser::synDownloadTableColumns(voice.tbl, 'audio_audio.m4a') %>%
#   unlist %>% data.frame 
# colnames(tbl.file.ids) = 'audioM4aFileLocation'
# tbl.file.ids[,'audio_audio.m4a'] = rownames(tbl.file.ids)
# 
# cntd.file.ids = synapseClient::synDownloadTableColumns(voice.tbl, 'audio_countdown.m4a') %>%
#   unlist %>% data.frame 
# colnames(cntd.file.ids) = 'countdownM4aFileLocation'
# cntd.file.ids[,'audio_countdown.m4a'] = rownames(cntd.file.ids)

## Join voice file location to the voice.tbl data frame
# voice.tbl = voice.tbl@values
# voice.tbl$row.id = rownames(voice.tbl)
# voice.tbl = voice.tbl %>%
#   tidyr::separate(row.id, into = c("source.row.id","source.row.version"), sep = "_") %>%
#   dplyr::mutate(source.row.id = as.integer(source.row.id),
#                 source.row.version = as.integer(source.row.version),
#                 row.id = paste(source.row.id, source.row.version, sep = '_')) %>%
#   dplyr::full_join(tbl.file.ids) %>%
#   dplyr::full_join(cntd.file.ids) %>%
#   dplyr::filter(!is.na(audio_audio.m4a),
#                 !is.na(audio_countdown.m4a))


voice.tbl.meta <- data.table::rbindlist(list(voice.tbl %>%
                                               dplyr::left_join(do.call(cbind, voice.json.loc[1]))),
                                        use.names = T, fill = T) %>%
  as.data.frame() 

voice.tbl.meta <- data.table::rbindlist(list(voice.tbl.meta %>%
                                               dplyr::left_join(do.call(cbind, voice.json.loc[2]))),
                                        use.names = T, fill = T) %>%
  as.data.frame %>% 
  dplyr::filter(!is.na(audio_audio.m4a),
                !is.na(audio_countdown.m4a)) %>% 
  dplyr::mutate(row.id = paste(ROW_ID, ROW_VERSION, sep = '_')) %>% 
  dplyr::rename(audioM4aFileLocation = audio_audio.fileLocation,
                countdownM4aFileLocation = audio_countdown.fileLocation,
                source.row.id = ROW_ID,
                source.row.version = ROW_VERSION)

# ## Get feature extracted record ids from synapse
# features.tbl = read.table(synGet(FEATURES_ID)$path, header = T, quote='', sep = '\t', stringsAsFactors = F) %>%
#   dplyr::mutate(createdOn = as.POSIXct(createdOn))
# 
# ## Get diff of voice and feature table
# voice.tbl = voice.tbl %>%
#   dplyr::filter(recordId %in% setdiff(voice.tbl$recordId, features.tbl$recordId))

if (dim(voice.tbl.meta)[1] != 0){
  ## Convert m4a to wav files and store the location to a file for matlab
  wav.file.locations = ddply(voice.tbl.meta, .(row.id), .fun = function(tbl){
    ## Convert audio m4a to wave files
    bname = basename(as.character(tbl$audioM4aFileLocation))
    dname = dirname(as.character(tbl$audioM4aFileLocation))
    
    tbl$audioWAVFileLocation = tryCatch({
      system(paste('avconv','-y','-i',tbl$audioM4aFileLocation,'-b 64k', paste0(dname, '/audio_audio.wav')))
      paste0(dname, '/audio_audio.wav')
    }, error = function(e){
      NA
    })
    
    ## Convert countdown m4a to wave files
    bname = basename(as.character(tbl$countdownM4aFileLocation))
    dname = dirname(as.character(tbl$countdownM4aFileLocation))
    
    tbl$countdownWAVFileLocation = tryCatch({
      system(paste('avconv','-y','-i',tbl$countdownM4aFileLocation,'-b 64k', paste0(dname, '/countdown_audio.wav')))
      paste0(dname, '/countdown_audio.wav')
    }, error = function(e){
      NA
    })
    
    return(tbl)
  })
  
  ## Write wav file names to a temp file for matlab program to read
  wav.file.locations %>% 
    dplyr::select(audioWAVFileLocation, countdownWAVFileLocation) %>%
    write.table(file = 'tempWavFileLocations.txt', sep = '\t', row.names = F, quote= F)
  
  ## Extract features (using MATLAB externally)
  input.file = paste0("'",WORKING_DIR,"/tempWavFileLocations.txt","'")
  output.file = paste0("'", WORKING_DIR, "/tempExtractedFeatures", "'")
  system(paste0('/Applications/MATLAB_R2019b.app/bin/matlab -nodesktop -nosplash -r "extract_voiceFeatures(',
                input.file,',', output.file, '); quit;"'))
  
  ## Read extracted features in to R
  cols = c('row.id', 'source.row.id', 'source.row.version', 'recordId', 'healthCode',
           'audio_audio.m4a', 'audio_countdown.m4a', 'audioM4aFileLocation', 'audioWAVFileLocation', 
           'countdownM4aFileLocation', 'countdownWAVFileLocation', 'createdOn', 'appVersion', 'phoneInfo', 
           'medTimepoint', 'Median_F0', 'Mean_Jitter', 'Median_Jitter', 'Mean_Shimmer', 'Median_Shimmer', 
           'MFCC_Band_1', 'MFCC_Band_2', 'MFCC_Band_3', 'MFCC_Band_4', 'MFCC_Jitter_Band_1_Positive',
           'MFCC_Jitter_Band_2_Positive', 'MFCC_Jitter_Band_3_Positive', 'MFCC_Jitter_Band_4_Positive')
  extracted.features = list(
    notFiltered = read.csv('tempExtractedFeaturesnotFiltered.csv') %>%
      dplyr::left_join(wav.file.locations) %>%
      dplyr::select(one_of(cols)),
    kalmanFiltered = read.csv('tempExtractedFeatureskalmanFiltered.csv') %>%
      dplyr::left_join(wav.file.locations) %>%
      dplyr::select(one_of(cols)),
    frequencyFiltered = read.csv('tempExtractedFeaturesfrequencyFiltered.csv') %>%
      dplyr::left_join(wav.file.locations) %>%
      dplyr::select(one_of(cols))) %>%
    data.table::rbindlist(use.names = T, fill = T, idcol = 'Filtering') %>%
    dplyr::arrange(source.row.id)
  
  ## Store extracted features in to synapse
  write.table(extracted.features, file = 'VoiceFeaturesAudio9.tsv', quote=F, row.names = F, sep = '\t')
  save(extracted.features, file = 'voice9.rdata')
  obj = File('VoiceFeaturesAudio.tsv', name = 'Voice Features Audio (Filtered)', parentId = FEATURES_PARENT_ID)
  obj = synStore(obj,
                 used = SOURCE_TBL,
                 activityName = activityName,
                 activityDescription = activityDescription,
                 executed = thisFile)
}
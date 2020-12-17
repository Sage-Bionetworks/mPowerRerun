#################################################################################
# GLOBALS                                                                       #
#################################################################################
GENERATE_FIGURE_PATH=R/GenerateFigure/
GENERATE_RAW_FEATURES_PATH=R/FeatureExtraction/
GENERATE_PROCESSED_FEATURES_PATH=R/FeatureEngineering/
GENERATE_INTERMEDIATE_PATH=R/Analyses/
GENERATE_FILTERED_HEALTHCODE_PATH=R/FilterHealthCode/

#################################################################################
# TARGET COMMANDS                                                             #
#################################################################################

all: clean_log project raw_features processed_features filtered_healthcodes analyses figures
regenerate_paper: clean_log project raw_features processed_features replicate_healthcodes analyses figures

# create directory for I/O when running in makefile
clean_log:
	rm -f pipeline.log;
	rm -f error.log;

# create template for project in synapse
project:
	Rscript R/buildProjectStructure.R

# generate raw features
raw_features:
	Rscript $(GENERATE_RAW_FEATURES_PATH)/getVoiceSubset.R || exit 1
	Rscript $(GENERATE_RAW_FEATURES_PATH)/extractTapFeatures.R  || exit 1
	Rscript $(GENERATE_RAW_FEATURES_PATH)/extractWalkFeatures.R || exit 1
	Rscript $(GENERATE_RAW_FEATURES_PATH)/extractRestFeatures.R || exit 1
	echo "FEATURE EXTRACTION DONE"

# clean raw features, and add any logic for used features
processed_features:
	for f in $(GENERATE_PROCESSED_FEATURES_PATH)/*; \
		do Rscript $$f || exit 1; \
		done; \
	echo "FEATURE CLEANING DONE!"

# filter healthcodes for subsampling users based on study
filtered_healthcodes:
	Rscript $(GENERATE_FILTERED_HEALTHCODE_PATH)/pdCaseControlsMatching.R || exit 1
	Rscript $(GENERATE_FILTERED_HEALTHCODE_PATH)/Nof1FilterHealthcodes.R  || exit 1
	Rscript $(GENERATE_FILTERED_HEALTHCODE_PATH)/identityConfoundingMatching.R || exit 1
	echo "FILTERING HEALTHCODE DONE!"

# replicate healthcode from previous usage due to unreproducible samples
replicate_healthcodes:
	Rscript $(GENERATE_FILTERED_HEALTHCODE_PATH)/convertHCfromRData.R || exit 1

# run all analyses codes
analyses:
	for f in $(GENERATE_INTERMEDIATE_PATH)/*; \
		do Rscript $$f || exit 1; \
		done; \
	echo "ANALYSIS DONE!"

# get all figures
figures:
	for f in $(GENERATE_FIGURE_PATH)/*; \
		do Rscript $$f || exit 1; \
		done; \
	echo "GENERATE FIGURES DONE!"
	

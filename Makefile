# Configuration
VERSION := 0.9.0
BLACKOUT_FILE := opera-disp-s1-blackout-dates-2025-02-13.json
DATE := $(shell date +%Y-%m-%d)
# Verbosely echo commands
SHELL = sh -xv

# Find the latest CMR survey file
CMR_SURVEY_TAR := $(shell ls -t cmr_survey*.csv.tar.gz | head -n1)
# Extract date range from CMR_SURVEY_TAR filename
# YYYY-mm-dd_to_YYYY-mm-dd
DATE_RANGE := $(shell echo $(CMR_SURVEY_TAR) | sed -n 's/.*\.\([0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}_to_[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}\).*/\1/p')
CMR_SURVEY_CSV := cmr_survey_$(DATE_RANGE).csv
# E.g.:
# echo cmr_survey.2016-07-01_to_2024-12-10.csv.tar.gz  | sed -n 's/.*\.\(.*\)\.csv\.tar\.gz/\1/p'
# 2016-07-01_to_2024-12-10
# Define consistent_bursts filename with date range
CONSISTENT_BURSTS := opera-disp-s1-consistent-burst-ids-$(DATE)-$(DATE_RANGE).json

# Define reference dates filename with today's date
REFERENCE_DATES := opera-disp-s1-reference-dates-$(DATE).json

# Main target
all: opera-s1-disp-$(VERSION).gpkg $(CONSISTENT_BURSTS) $(REFERENCE_DATES)

# Create Opera DB
opera-s1-disp-$(VERSION).gpkg:
	opera-db create

# Extract CMR survey
$(CMR_SURVEY_CSV): $(CMR_SURVEY_TAR)
	tar -xzf $< -O > $@

# Make burst catalog
# WANT THIS TO BE LIKE: opera-disp-s1-consistent-burst-ids-2024-10-11-2016-07-01_to_2024-09-04.json
$(CONSISTENT_BURSTS): $(CMR_SURVEY_CSV) opera-s1-disp-$(VERSION).gpkg
	opera-db make-burst-catalog --blackout-file $(BLACKOUT_FILE) $^ > $@

# Make reference dates
$(REFERENCE_DATES): $(BLACKOUT_FILE)
	opera-db make-reference-dates --output $@ --blackout-file $(BLACKOUT_FILE)

# Clean up intermediate files
clean:
	rm -f $(CMR_SURVEY_CSV) $(CONSISTENT_BURSTS)

# Clean all generated files
cleanall: clean
	rm -f opera-s1-disp-$(VERSION).gpkg opera-disp-s1-consistent-bursts-*.json

.PHONY: all clean cleanall

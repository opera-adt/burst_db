# Configuration
# VERSION can be passed in: make VERSION=0.13.0
# Or defaults to reading from git tag
VERSION ?= $(shell git describe --tags --abbrev=0 | sed 's/^v//')
ifeq ($(VERSION),)
    $(error No VERSION specified and no git tags found. Use: make VERSION=0.13.0)
endif

SNOW_PARQUET := ../snow-analysis/opera-region4-snow-analysis.parquet
DATE := $(shell date +%Y-%m-%d)
# Global blackout period to add to all frames
GLOBAL_BLACKOUT_START := 2025-04-29T19:40:10
GLOBAL_BLACKOUT_END := 2025-05-01T19:33:34
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
	@echo "================================================"
	@echo "Build complete for version $(VERSION)"
	@echo "================================================"

# Create Opera DB
opera-s1-disp-$(VERSION).gpkg:
	opera-db create

# Extract CMR survey
$(CMR_SURVEY_CSV): $(CMR_SURVEY_TAR)
	tar -xzf $< -O > $@

# Create blackout dates from snow analysis, then add global blackout period
BLACKOUT_FILE := opera-disp-s1-blackout-dates-$(DATE).json
$(BLACKOUT_FILE): $(SNOW_PARQUET)
	opera-db create-blackout $(SNOW_PARQUET)
	opera-db add-global-blackout $@ \
		--start-date "$(GLOBAL_BLACKOUT_START)" \
		--end-date "$(GLOBAL_BLACKOUT_END)" \
		--output-file $@.tmp
	mv $@.tmp $@

# Make burst catalog
# E.g.: opera-disp-s1-consistent-burst-ids-2024-10-11-2016-07-01_to_2024-09-04.json
# Also we make one without blackout dates for comparison
$(CONSISTENT_BURSTS): $(CMR_SURVEY_CSV) opera-s1-disp-$(VERSION).gpkg $(BLACKOUT_FILE)
	opera-db make-burst-catalog $(CMR_SURVEY_CSV) opera-s1-disp-$(VERSION).gpkg
	mv $(CONSISTENT_BURSTS) opera-disp-s1-consistent-burst-ids-no-blackout.json
	opera-db make-burst-catalog --blackout-file $(BLACKOUT_FILE)  $(CMR_SURVEY_CSV) opera-s1-disp-$(VERSION).gpkg

# Make reference dates
$(REFERENCE_DATES): $(BLACKOUT_FILE)
	opera-db make-reference-dates --output $@ --blackout-file $(BLACKOUT_FILE)

# Clean up intermediate files
clean:
	rm -f $(CMR_SURVEY_CSV) $(CONSISTENT_BURSTS)

# Clean all generated files
cleanall: clean
	rm -f opera-s1-disp-$(VERSION).gpkg opera-disp-s1-consistent-bursts-*.json \
		opera-disp-s1-blackout-dates-*.json *.duckdb

# Show current version
show-version:
	@echo "Current version: $(VERSION)"

# Show current configuration
show-config:
	@echo "================================================"
	@echo "Build Configuration"
	@echo "================================================"
	@echo "VERSION: $(VERSION)"
	@echo "DATE: $(DATE)"
	@echo "CMR_SURVEY_CSV: $(CMR_SURVEY_CSV)"
	@echo "DATE_RANGE: $(DATE_RANGE)"
	@echo "GLOBAL_BLACKOUT_START: $(GLOBAL_BLACKOUT_START)"
	@echo "GLOBAL_BLACKOUT_END: $(GLOBAL_BLACKOUT_END)"
	@echo "BLACKOUT_FILE: $(BLACKOUT_FILE)"
	@echo "CONSISTENT_BURSTS: $(CONSISTENT_BURSTS)"
	@echo "REFERENCE_DATES: $(REFERENCE_DATES)"
	@echo "================================================"

.PHONY: all clean cleanall show-version show-config

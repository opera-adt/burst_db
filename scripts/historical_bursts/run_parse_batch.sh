#!/bin/bash
set -e

# AWS_BATCH_JOB_ARRAY_INDEX is assumed to be set
# For testing, you can uncomment the line below
# AWS_BATCH_JOB_ARRAY_INDEX=2

# Get the "START_INDEX" to add to the aws batch index
# make it 0 if it's not set
START_INDEX=${START_INDEX:-0}
DAYS_TO_ADD=$((AWS_BATCH_JOB_ARRAY_INDEX + START_INDEX))
# DATE_TO_DOWNLOAD=$(date -d "2015-01-01 + 10 days" "+%Y%m%d")
DATE_TO_DOWNLOAD=$(date -d "2015-01-01 + $DAYS_TO_ADD days" "+%Y%m%d")


# parse_bursts.py [-h] [--out-dir OUT_DIR] [--bucket BUCKET] [--in-folder IN_FOLDER] [--out-folder OUT_FOLDER] [--start-date START_DATE] [--end-date END_DATE]
#                        [--max-workers MAX_WORKERS] [--satellite {A,B}] [--single-or-double {*,S,D}] [--no-clean]

# Get the satellite (A or B) from the env variable
SATELLITE=${SATELLITE:-A}

# Set the output folder
out_folder="${OUT_FOLDER:-burst_csvs}"


echo /src/parse_bursts.py --bucket burst-database-scott --in-folder safe_folders_zipped \
    --out-folder "$out_folder" --start-date "$DATE_TO_DOWNLOAD" --end-date "$DATE_TO_DOWNLOAD" \
    --satellite "$SATELLITE" --max-workers 2
# Run the Python script
python /src/parse_bursts.py --bucket burst-database-scott --in-folder safe_folders_zipped \
    --out-folder "$out_folder" --start-date "$DATE_TO_DOWNLOAD" --end-date "$DATE_TO_DOWNLOAD" \
    --satellite "$SATELLITE" --max-workers 2

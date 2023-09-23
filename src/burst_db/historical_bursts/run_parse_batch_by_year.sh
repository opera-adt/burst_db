#!/bin/bash
set -e


# usage: parse_bursts.py [-h] [--out-dir OUT_DIR] [--bucket BUCKET] [--in-folder IN_FOLDER]
#                        [--out-folder OUT_FOLDER] [--start-date START_DATE] [--end-date END_DATE]
#                        [--max-workers MAX_WORKERS] [--satellite {A,B}] [--single-or-double {*,S,D}] [--no-clean]
#                        [--combine-by-date] [--full-safe-list FULL_SAFE_LIST] [--no-upload]

# AWS_BATCH_JOB_ARRAY_INDEX is assumed to be set
# For testing, you can uncomment the line below
# AWS_BATCH_JOB_ARRAY_INDEX=2

# Get the "DAYS_OFFSET" to add to the aws batch index
# make it 0 if it's not set
DAYS_OFFSET=${DAYS_OFFSET:-0}
# Number of days to run for each job
DAYS_PER_JOB=${DAYS_PER_JOB:-50}
days_to_add=$((AWS_BATCH_JOB_ARRAY_INDEX * DAYS_PER_JOB + DAYS_OFFSET))

# Get the satellite (A or B) from the env variable
SATELLITE=${SATELLITE:-A}
# Start the S1A at 2014-10-01, but S1B at 2016-08-20
if [ "$SATELLITE" = "A" ]; then
    day0="2014-10-01"
elif [ "$SATELLITE" = "B" ]; then
    day0="2016-08-20"
else
    echo "SATELLITE must be A or B"
    exit 1
fi

start_date=$(date -d "$day0 + $days_to_add days" "+%Y%m%d")
end_date=$(date -d "$day0 + $((days_to_add - 1)) days + $DAYS_PER_JOB days" "+%Y%m%d")

# years_to_add=$((AWS_BATCH_JOB_ARRAY_INDEX + START_INDEX))

# start_date=$(date -d "2015-01-01 + $years_to_add years" "+%Y%m%d")
# end_date=$(date -d "2015-01-01 + $years_to_add years + 1 year" "+%Y%m%d")



# Set the output folder
out_folder="${OUT_FOLDER:-bursts-by-date}"

# unzip /src/all_cmr_stac_safes_sorted.txt.zip # unzip not installed
python -m zipfile -e /src/all_cmr_stac_safes_sorted.txt.zip .

# Echo command and run
echo /src/parse_bursts.py --bucket burst-database-scott --in-folder safe_folders_zipped \
    --out-folder "$out_folder" --start-date "$start_date" --end-date "$end_date" \
    --combine-by-date --satellite "$SATELLITE" --max-workers 3 \
    --full-safe-list all_cmr_stac_safes_sorted.txt
python /src/parse_bursts.py --bucket burst-database-scott --in-folder safe_folders_zipped \
    --out-folder "$out_folder" --start-date "$start_date" --end-date "$end_date" \
    --combine-by-date --satellite "$SATELLITE" --max-workers 3 \
    --full-safe-list all_cmr_stac_safes_sorted.txt
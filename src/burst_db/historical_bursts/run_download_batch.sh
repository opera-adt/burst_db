#!/bin/bash
set -e

# AWS_BATCH_JOB_ARRAY_INDEX is assumed to be set
# For testing, you can uncomment the line below
# AWS_BATCH_JOB_ARRAY_INDEX=2

bucket="burst-database-scott"
out_folder="safe_folders_zipped"

# Get the "START_INDEX" to add to the aws batch index
# make it 0 if it's not set
START_INDEX=${START_INDEX:-0}

REMAINING_FOLDER="${REMAINING_FOLDER:-remaining}"
MAX_WORKERS=${MAX_WORKERS:-2}
BATCH_SIZE=${BATCH_SIZE:-20}
BATCH_INDEX_START=${BATCH_INDEX_START:-0}

# Convert index to zero-padded string
# zero_padded_index=$(printf "%04d" $AWS_BATCH_JOB_ARRAY_INDEX)
zero_padded_index=$(printf "%04d" $(($AWS_BATCH_JOB_ARRAY_INDEX + $START_INDEX)))

# Build the S3 path for the "safe_list"
# We have divided the total list into smaller text files named "part_0000.txt", "part_0001.txt", etc.
part_file_name="part_${zero_padded_index}.txt"
remote_filename="s3://${bucket}/${REMAINING_FOLDER}/${part_file_name}"

# Download safe_list from S3 to a local file (assuming you want to download it)
echo "Downloading $remote_filename to $part_file_name"
aws s3 cp "$remote_filename"  "$part_file_name"

# Run the Python script
echo /src/download.py --batch-size "$BATCH_SIZE" --start-idx "$BATCH_INDEX_START" --max-work "$MAX_WORKERS" --bucket ${bucket} --folder ${out_folder} --safe-list "$part_file_name"
python /src/download.py --batch-size "$BATCH_SIZE" --start-idx "$BATCH_INDEX_START" --max-work "$MAX_WORKERS" --bucket ${bucket} --folder ${out_folder} --safe-list "$part_file_name"

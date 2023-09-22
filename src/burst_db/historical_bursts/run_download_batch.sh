#!/bin/bash
set -e

# AWS_BATCH_JOB_ARRAY_INDEX is assumed to be set
# For testing, you can uncomment the line below
# AWS_BATCH_JOB_ARRAY_INDEX=2

# Get the "START_INDEX" to add to the aws batch index
# make it 0 if it's not set
START_INDEX=${START_INDEX:-0}

# Convert index to zero-padded string
# zero_padded_index=$(printf "%04d" $AWS_BATCH_JOB_ARRAY_INDEX)
zero_padded_index=$(printf "%04d" $(($AWS_BATCH_JOB_ARRAY_INDEX + $START_INDEX)))

# Build the S3 path for the "safe_list"
part_file_name="part_${zero_padded_index}.txt"
remote_filename="s3://burst-database-scott/remaining/${part_file_name}"

# Download safe_list from S3 to a local file (assuming you want to download it)
echo "Downloading $remote_filename to $part_file_name"
aws s3 cp "$remote_filename"  "$part_file_name"

# Run the Python script
python /src/download.py --batch 20 --max-work 5 --bucket burst-database-scott --folder safe_folders_zipped3 --safe-list "$part_file_name"

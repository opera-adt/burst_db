#!/bin/bash
set -e

# AWS_BATCH_JOB_ARRAY_INDEX is assumed to be set
# For testing, you can uncomment the line below
# AWS_BATCH_JOB_ARRAY_INDEX=2

# Convert index to zero-padded string
zero_padded_index=$(printf "%04d" $AWS_BATCH_JOB_ARRAY_INDEX)

# Build the S3 path for the "safe_list"
part_file_name="part_${zero_padded_index}.txt"
remove_filename="s3://burst-database-scott/remaining/${part_file_name}"

# Download safe_list from S3 to a local file (assuming you want to download it)
aws s3 cp "$remove_filename"  "$part_file_name"

# Run the Python script
python /src/download.py --batch 100 --max-work 10 --bucket burst-database-scott --folder safe_folders_zipped3 --safe-list "$part_file_name"

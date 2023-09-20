#!/bin/bash
set -e

increment=40
max_value=3176
start_idx=135
safe_list="$1"

for ((start=start_idx, end=start_idx+increment; end <= max_value; start+=increment, end+=increment)); do
    date
    echo "start: $start, end: $end"
    ASFSMD_CLIENT=s3fs python download.py --batch 100 --max-work 20 --bucket burst-database-scott --folder safe_folders_zipped2 --safe-list "$safe_list" --start $start --end $end
done
``
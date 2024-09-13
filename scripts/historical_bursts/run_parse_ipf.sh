#!/bin/bash
set -e

for year in $(seq 2014 2023); do
    for month in $(seq -w 1 12); do
        echo "Processing $year/$month..."
        python ~/repos/burst_db/src/burst_db/historical_bursts/parse_ipf.py "zipped-metadata/$year/$month/*/*.zip" --out "granule_ipfs_${year}_${month}.csv"
    done
done

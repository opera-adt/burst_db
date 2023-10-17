#!/usr/bin/env python
"""
A script to export the burst geogrid from the augmented burst map

"""
import argparse
import os

import burst_database_core as bd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extracts burst geogrid data from"
            "augmented burstmap. Writes out the data into "
            ".csv, .json, or .sqlite"
        )
    )
    parser.add_argument(
        "path_augmented_burstmap",
        type=str,
        help=(
            "Path to the augmented burst map from which "
            "the burst geogrid will be extracted"
        ),
    )
    parser.add_argument(
        "path_database_out",
        type=str,
        help="Path to the output database file\nSupported format: .csv, .json, .sqlite",
    )

    args = parser.parse_args()

    if not os.path.exists(args.path_augmented_burstmap):
        raise ValueError(
            f"Cannot find the augmented burst map: {args.path_augmented_burstmap}"
        )

    records_burst_data = bd.extract_burst_geogrid_data(args.path_augmented_burstmap)

    if args.path_database_out.endswith(".csv"):
        bd.export_to_csv(records_burst_data, args.path_database_out)
    elif args.path_database_out.endswith(".json"):
        bd.export_to_json(records_burst_data, args.path_database_out)
    elif (args.path_database_out.endswith(".sqlite")) or (
        args.path_database_out.endswith(".sqlite3")
    ):
        bd.export_to_sqlite(records_burst_data, args.path_database_out)
    else:
        raise ValueError(
            f"Output file format was not supported:{args.path_database_out}"
        )

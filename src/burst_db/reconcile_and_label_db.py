#!/usr/bin/env python3
r"""Standalone script to reconcile burst databases and add processing mode labels.

This script performs two main tasks:
1. Reconcile differences between old and new burst database JSON files
2. Add processing_mode labels (historical/forward) to sensing times

Usage:
    python reconcile_and_label_db.py --old-db old.json \\
        --new-db new.json --output output.json
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any


def load_burst_db(filepath: str | Path) -> dict[str, Any]:
    """Load a burst database JSON file.

    Handles both nested (metadata/data) and flat structures.

    Parameters
    ----------
    filepath : str | Path
        Path to the JSON file.

    Returns
    -------
    dict
        The loaded JSON data.

    """
    with open(filepath, "r") as f:
        loaded = json.load(f)
    return loaded


def get_data_section(db: dict[str, Any]) -> dict[str, Any]:
    """Extract the data section from a burst database.

    Parameters
    ----------
    db : dict
        The loaded burst database.

    Returns
    -------
    dict
        The data section containing frame information.

    """
    return db.get("data", db)


def parse_sensing_time(time_str: str) -> datetime:
    """Parse an ISO 8601 sensing time string.

    Parameters
    ----------
    time_str : str
        Time string in format YYYY-MM-DDTHH:MM:SS.

    Returns
    -------
    datetime
        Parsed datetime object.

    """
    return datetime.fromisoformat(time_str)


def get_date_only(time_str: str) -> str:
    """Extract the date portion from a sensing time string.

    Parameters
    ----------
    time_str : str
        Time string in format YYYY-MM-DDTHH:MM:SS.

    Returns
    -------
    str
        Date string in format YYYY-MM-DD.

    """
    return time_str.split("T")[0]


def format_sensing_time(dt: datetime) -> str:
    """Format a datetime to ISO 8601 sensing time string.

    Parameters
    ----------
    dt : datetime
        Datetime object to format.

    Returns
    -------
    str
        Formatted string in YYYY-MM-DDTHH:MM:SS format.

    """
    return dt.strftime("%Y-%m-%dT%H:%M:%S")


def find_frames_with_different_bursts(
    old_data: dict[str, Any], new_data: dict[str, Any]
) -> list[str]:
    """Find frame IDs where burst IDs differ between old and new databases.

    Parameters
    ----------
    old_data : dict
        Data section from old database.
    new_data : dict
        Data section from new database.

    Returns
    -------
    list[str]
        List of frame IDs with different burst IDs.

    """
    different_frames = []

    # Check all frames in old database
    for frame_id, old_frame_data in old_data.items():
        if frame_id not in new_data:
            continue

        old_bursts = set(old_frame_data.get("burst_id_list", []))
        new_bursts = set(new_data[frame_id].get("burst_id_list", []))

        if old_bursts != new_bursts:
            different_frames.append(frame_id)

    return different_frames


def reconcile_frame(
    old_frame: dict[str, Any], new_frame: dict[str, Any]
) -> tuple[dict[str, Any], bool]:
    """Reconcile a single frame between old and new databases.

    If new has more burst IDs than old AND sensing times have no overlap
    (indicating a complete restart after a gap), keep new data as-is.
    Otherwise:
    - If new has more burst IDs than old, remove the extras.
    - If new is missing sensing times from old, add them.

    Parameters
    ----------
    old_frame : dict
        Frame data from old database.
    new_frame : dict
        Frame data from new database.

    Returns
    -------
    tuple[dict, bool]
        Reconciled frame data and a boolean indicating if reconciliation was skipped
        (True if kept as-is due to no overlap).

    """
    old_bursts = old_frame.get("burst_id_list", [])
    new_bursts = new_frame.get("burst_id_list", [])
    old_times = old_frame.get("sensing_time_list", [])
    new_times = new_frame.get("sensing_time_list", [])

    old_times_set = set(old_times)
    new_times_set = set(new_times)

    # For overlap and missing time checks, compare only dates (ignore hours)
    old_dates_set = {get_date_only(t) for t in old_times}
    new_dates_set = {get_date_only(t) for t in new_times}

    # Check if there's any overlap in sensing times (by date only)
    dates_overlap = old_dates_set & new_dates_set

    # If new has more burst IDs AND no overlap in sensing times,
    # this is a complete restart after a gap - keep new data as-is
    if len(new_bursts) > len(old_bursts) and not dates_overlap:
        return new_frame.copy(), True

    reconciled = {}

    # Handle burst IDs: if new has more than old, use old's burst IDs
    if len(new_bursts) > len(old_bursts):
        reconciled["burst_id_list"] = old_bursts.copy()
    else:
        reconciled["burst_id_list"] = new_bursts.copy()

    # Handle sensing times: start with new, add missing from old
    # Find sensing times in old where the DATE is not present in new
    # (if the date matches but hours differ, we don't consider it missing)
    missing_times = {t for t in old_times if get_date_only(t) not in new_dates_set}

    # Combine and sort all sensing times
    all_times = list(new_times_set | missing_times)
    all_times_sorted = sorted(all_times, key=parse_sensing_time)
    reconciled["sensing_time_list"] = all_times_sorted

    return reconciled, False


def reconcile_databases(
    old_db: dict[str, Any], new_db: dict[str, Any]
) -> dict[str, Any]:
    """Reconcile old and new burst databases.

    For frames with different burst IDs:
    - If new has more burst IDs AND no overlap in sensing times (complete restart),
      keep new data as-is (just label it)
    - Otherwise, if new has more burst IDs, remove extras (use old's burst IDs)
    - Add missing sensing times from old to new

    Parameters
    ----------
    old_db : dict
        Old burst database.
    new_db : dict
        New burst database.

    Returns
    -------
    dict
        Reconciled database with metadata preserved from new.

    """
    old_data = get_data_section(old_db)
    new_data = get_data_section(new_db)

    # Start with a copy of new database structure
    reconciled_db = {}
    if "metadata" in new_db:
        reconciled_db["metadata"] = new_db["metadata"].copy()
        reconciled_db["metadata"]["reconciled_from"] = {
            "old_db": "provided",
            "new_db": "provided",
            "reconciliation_time": datetime.now().isoformat(),
        }

    reconciled_data = {}

    # Find frames that need reconciliation
    different_frames = find_frames_with_different_bursts(old_data, new_data)
    print(f"Found {len(different_frames)} frames with different burst IDs")

    skipped_count = 0
    reconciled_count = 0

    for frame_id in new_data:
        if frame_id in different_frames and frame_id in old_data:
            # Try to reconcile this frame
            reconciled_frame, was_skipped = reconcile_frame(
                old_data[frame_id], new_data[frame_id]
            )
            reconciled_data[frame_id] = reconciled_frame

            if was_skipped:
                skipped_count += 1
                print(
                    f"  Frame {frame_id}: kept as-is (no overlap, new burst IDs - "
                    f"likely restart after gap)"
                )
            else:
                reconciled_count += 1
                print(
                    f"  Frame {frame_id}: reconciled "
                    f"(bursts: {len(new_data[frame_id].get('burst_id_list', []))} -> "
                    f"{len(reconciled_frame['burst_id_list'])}, "
                    f"times: {len(new_data[frame_id].get('sensing_time_list', []))} -> "
                    f"{len(reconciled_frame['sensing_time_list'])})"
                )
        else:
            # Keep new frame data as-is
            reconciled_data[frame_id] = new_data[frame_id].copy()

    print(f"  Reconciled: {reconciled_count}, Kept as-is (no overlap): {skipped_count}")

    reconciled_db["data"] = reconciled_data
    return reconciled_db


def identify_time_groups(
    sensing_times: list[str], gap_threshold_years: float = 2.0
) -> list[list[str]]:
    """Split sensing times into groups based on temporal gaps.

    If there's a gap of >= gap_threshold_years between consecutive sensing times,
    split into separate groups.

    Parameters
    ----------
    sensing_times : list[str]
        List of sensing time strings in ISO 8601 format.
    gap_threshold_years : float, optional
        Gap threshold in years to split groups. Default is 2.0.

    Returns
    -------
    list[list[str]]
        List of groups, where each group is a list of sensing times.

    """
    if not sensing_times:
        return []

    # Sort sensing times
    sorted_times = sorted(sensing_times, key=parse_sensing_time)

    gap_threshold = timedelta(days=gap_threshold_years * 365)

    groups = []
    current_group = [sorted_times[0]]

    for i in range(1, len(sorted_times)):
        prev_time = parse_sensing_time(sorted_times[i - 1])
        curr_time = parse_sensing_time(sorted_times[i])
        gap = curr_time - prev_time

        if gap >= gap_threshold:
            # Start a new group
            groups.append(current_group)
            current_group = [sorted_times[i]]
        else:
            current_group.append(sorted_times[i])

    # Don't forget the last group
    if current_group:
        groups.append(current_group)

    return groups


def assign_processing_modes(
    sensing_times: list[str], batch_size: int = 15, gap_threshold_years: float = 2.0
) -> dict[str, str]:
    """Assign processing modes to sensing times.

    Logic:
    - If a group has fewer than batch_size sensing times total, all are "no_run"
      (not enough data to form a complete batch)
    - Otherwise, organize sensing times in batches of batch_size from the beginning
    - If a batch has full batch_size sensing times, all are "historical"
    - If a batch has less than batch_size (usually the last batch), all are "forward"
    - If there's a gap of >= gap_threshold_years, split into separate batch sets
      and restart the labeling
    - Labels are numbered (e.g., "historical_01", "forward_01") to indicate
      which group they belong to. The number increments after each gap.
      "no_run" labels do not have numeric suffixes.

    Parameters
    ----------
    sensing_times : list[str]
        List of sensing time strings in ISO 8601 format.
    batch_size : int, optional
        Number of sensing times per batch. Default is 15.
    gap_threshold_years : float, optional
        Gap threshold in years to restart batch counting. Default is 2.0.

    Returns
    -------
    dict[str, str]
        Dictionary mapping sensing times to their processing modes.
        E.g., {"2016-07-02T23:05:35": "historical_01", ...}

    """
    if not sensing_times:
        return {}

    # Sort sensing times
    sorted_times = sorted(sensing_times, key=parse_sensing_time)

    # Split into groups based on gaps
    groups = identify_time_groups(sorted_times, gap_threshold_years)

    # Create mapping from sensing time to mode
    time_to_mode = {}

    for group_num, group in enumerate(groups, start=1):
        num_times = len(group)
        # Format group number with zero-padding (e.g., 01, 02, ...)
        group_suffix = f"_{group_num:02d}"

        # If group has fewer than batch_size items, all are "no_run"
        # Not enough data to form a complete batch
        if num_times < batch_size:
            for time_str in group:
                time_to_mode[time_str] = "no_run"
            continue

        # Process this group with batch logic
        num_full_batches = num_times // batch_size

        for i, time_str in enumerate(group):
            batch_index = i // batch_size
            is_full_batch = batch_index < num_full_batches

            if is_full_batch:
                time_to_mode[time_str] = f"historical{group_suffix}"
            else:
                # Last partial batch
                time_to_mode[time_str] = f"forward{group_suffix}"

    # Return dictionary mapping sensing times to modes (sorted by time)
    return {t: time_to_mode[t] for t in sorted_times}


def add_processing_modes(
    db: dict[str, Any], batch_size: int = 15, gap_threshold_years: float = 2.0
) -> dict[str, Any]:
    """Add processing_mode to all frames in the database.

    Parameters
    ----------
    db : dict
        Burst database with data section.
    batch_size : int, optional
        Number of sensing times per batch. Default is 15.
    gap_threshold_years : float, optional
        Gap threshold in years to restart batch counting. Default is 2.0.

    Returns
    -------
    dict
        Database with processing_mode added to each frame.

    """
    result = {}

    # Copy metadata if present
    if "metadata" in db:
        result["metadata"] = db["metadata"].copy()
        result["metadata"]["processing_mode_params"] = {
            "batch_size": batch_size,
            "gap_threshold_years": gap_threshold_years,
            "labeling_time": datetime.now().isoformat(),
        }

    data = get_data_section(db)
    result_data = {}

    for frame_id, frame_data in data.items():
        sensing_times = frame_data.get("sensing_time_list", [])
        # assign_processing_modes now returns a dict {sensing_time: label}
        sensing_time_labels = assign_processing_modes(
            sensing_times, batch_size, gap_threshold_years
        )

        result_data[frame_id] = {
            "burst_id_list": frame_data.get("burst_id_list", []),
            "sensing_time_list": sensing_time_labels,
        }

    result["data"] = result_data
    return result


def get_processing_mode_summary(db: dict[str, Any]) -> dict[str, Any]:
    """Get a summary of processing mode distribution.

    Parameters
    ----------
    db : dict
        Database with sensing_time_list as dict {time: label}.

    Returns
    -------
    dict
        Summary statistics.

    """
    data = get_data_section(db)

    total_frames = len(data)
    total_historical = 0
    total_forward = 0
    total_no_run = 0
    frames_with_gaps = 0
    max_group_num = 0

    for _frame_id, frame_data in data.items():
        # sensing_time_list is now a dict {sensing_time: label}
        sensing_time_dict = frame_data.get("sensing_time_list", {})
        times = list(sensing_time_dict.keys())
        modes = list(sensing_time_dict.values())

        # Count modes that start with "historical", "forward", or "no_run"
        total_historical += sum(1 for m in modes if m.startswith("historical"))
        total_forward += sum(1 for m in modes if m.startswith("forward"))
        total_no_run += sum(1 for m in modes if m.startswith("no_run"))

        # Track max group number
        for m in modes:
            if "_" in m:
                try:
                    group_num = int(m.split("_")[-1])
                    max_group_num = max(max_group_num, group_num)
                except ValueError:
                    pass

        # Check for gaps
        groups = identify_time_groups(times)
        if len(groups) > 1:
            frames_with_gaps += 1

    return {
        "total_frames": total_frames,
        "total_sensing_times": total_historical + total_forward + total_no_run,
        "historical_count": total_historical,
        "forward_count": total_forward,
        "no_run_count": total_no_run,
        "frames_with_temporal_gaps": frames_with_gaps,
        "max_group_number": max_group_num,
    }


def save_database(db: dict[str, Any], filepath: str | Path) -> None:
    """Save a burst database to JSON file.

    Parameters
    ----------
    db : dict
        Database to save.
    filepath : str | Path
        Output file path.

    """
    with open(filepath, "w") as f:
        json.dump(db, f, indent=2, default=str)
    print(f"Saved database to {filepath}")


def main():
    """Run the reconciliation and labeling workflow."""
    parser = argparse.ArgumentParser(
        description="Reconcile burst databases and add processing mode labels.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Reconcile two databases and add processing modes
    python reconcile_and_label_db.py --old-db old.json --new-db new.json \\
           --output output.json

    # Reconcile and also update the input (new) database with corrections
    python reconcile_and_label_db.py --old-db old.json --new-db new.json \\
           --output output.json --update-input

    # Only add processing modes to a single database (no reconciliation)
    python reconcile_and_label_db.py --new-db input.json --output output.json \\
           --no-reconcile

    # Customize batch size and gap threshold
    python reconcile_and_label_db.py --old-db old.json --new-db new.json \\
        --output output.json \\
        --batch-size 20 --gap-threshold 1.5
        """,
    )

    parser.add_argument(
        "--old-db",
        type=str,
        help="Path to the old burst database JSON file (for reconciliation).",
    )
    parser.add_argument(
        "--new-db",
        type=str,
        required=True,
        help="Path to the new burst database JSON file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path for the output JSON file.",
    )
    parser.add_argument(
        "--no-reconcile",
        action="store_true",
        help="Skip reconciliation step, only add processing modes.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=15,
        help="Number of sensing times per batch for processing mode assignment. "
        "Default: 15",
    )
    parser.add_argument(
        "--gap-threshold",
        type=float,
        default=2.0,
        help="Gap threshold in years to restart batch counting. Default: 2.0",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print detailed information during processing.",
    )
    parser.add_argument(
        "--update-input",
        action="store_true",
        help="Update the input (new) database file with reconciled burst IDs and "
        "sensing times. This saves the corrections back to the original file "
        "before processing mode labels are added.",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.no_reconcile and not args.old_db:
        parser.error("--old-db is required unless --no-reconcile is specified")

    print("=" * 60)
    print("Burst Database Reconciliation and Labeling")
    print("=" * 60)

    # Load the new database
    print(f"\nLoading new database: {args.new_db}")
    new_db = load_burst_db(args.new_db)
    print(f"  Loaded {len(get_data_section(new_db))} frames")

    if args.no_reconcile:
        # Skip reconciliation, just use new database
        working_db = new_db
        print("\nSkipping reconciliation (--no-reconcile specified)")
    else:
        # Load old database and reconcile
        print(f"\nLoading old database: {args.old_db}")
        old_db = load_burst_db(args.old_db)
        print(f"  Loaded {len(get_data_section(old_db))} frames")

        print("\n" + "-" * 60)
        print("Task 1: Reconciling databases")
        print("-" * 60)
        working_db = reconcile_databases(old_db, new_db)

        # Update the input file if requested
        if args.update_input:
            print("\n" + "-" * 60)
            print("Updating input database with reconciled data")
            print("-" * 60)
            # Create updated version preserving original metadata, only updating data
            updated_input = {}
            if "metadata" in new_db:
                updated_input["metadata"] = new_db["metadata"]
            updated_input["data"] = {}
            working_data = get_data_section(working_db)
            for frame_id, frame_data in working_data.items():
                # Convert sensing_time_list back to a list if it's a dict
                sensing_times = frame_data.get("sensing_time_list", [])
                if isinstance(sensing_times, dict):
                    sensing_times = list(sensing_times.keys())
                updated_input["data"][frame_id] = {
                    "burst_id_list": frame_data.get("burst_id_list", []),
                    "sensing_time_list": sensing_times,
                }
            save_database(updated_input, args.new_db)
            print(f"  Updated input database: {args.new_db}")

    print("\n" + "-" * 60)
    print("Task 2: Adding processing mode labels")
    print("-" * 60)
    print(f"  Batch size: {args.batch_size}")
    print(f"  Gap threshold: {args.gap_threshold} years")

    result_db = add_processing_modes(
        working_db,
        batch_size=args.batch_size,
        gap_threshold_years=args.gap_threshold,
    )

    # Print summary
    summary = get_processing_mode_summary(result_db)
    print("\nProcessing mode summary:")
    print(f"  Total frames: {summary['total_frames']}")
    print(f"  Total sensing times: {summary['total_sensing_times']}")
    print(f"  Historical: {summary['historical_count']}")
    print(f"  Forward: {summary['forward_count']}")
    print(f"  No run: {summary['no_run_count']}")
    print(
        f"  Frames with temporal gaps (>= {args.gap_threshold} years): "
        f"{summary['frames_with_temporal_gaps']}"
    )
    print(f"  Max group number: {summary['max_group_number']}")

    # Verbose output
    if args.verbose:
        print("\n" + "-" * 60)
        print("Detailed frame information:")
        print("-" * 60)
        data = get_data_section(result_db)
        for frame_id, frame_data in list(data.items())[:5]:  # Show first 5
            # sensing_time_list is now a dict {sensing_time: label}
            sensing_time_dict = frame_data.get("sensing_time_list", {})
            times = list(sensing_time_dict.keys())
            modes = list(sensing_time_dict.values())
            groups = identify_time_groups(times)

            # Count historical, forward, and no_run (with any suffix)
            hist_count = sum(1 for m in modes if m.startswith("historical"))
            fwd_count = sum(1 for m in modes if m.startswith("forward"))
            no_run_count = sum(1 for m in modes if m.startswith("no_run"))

            # Get unique mode labels
            unique_modes = sorted(set(modes))

            print(f"\n  Frame {frame_id}:")
            print(f"    Burst IDs: {len(frame_data.get('burst_id_list', []))}")
            print(f"    Sensing times: {len(times)}")
            print(f"    Temporal groups: {len(groups)}")
            print(f"    Historical: {hist_count}, \
                      Forward: {fwd_count}, \
                      No run: {no_run_count}")
            print(f"    Labels: {unique_modes}")

            if len(groups) > 1:
                print("    Group sizes:", [len(g) for g in groups])
        if len(data) > 5:
            print(f"\n  ... and {len(data) - 5} more frames")

    # Save result
    print("\n" + "-" * 60)
    save_database(result_db, args.output)
    print("=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()

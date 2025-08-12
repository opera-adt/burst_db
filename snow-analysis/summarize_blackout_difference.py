"""Summarize of "consistent-burst-id" JSON before and after blackout-date filtering."""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rich.console import Console
from rich.table import Table

PRIORITY_F = Path("src/burst_db/data/NApriorityrollout_framebased_v8_13Mar2025.geojson")
REGION3B_F = Path("src/burst_db/data/region3b_v1_23Jul2025.json")


def _load_consistent_json(path: Path) -> pd.DataFrame:
    """Return a DataFrame with sensing/burst counts for each frame."""
    data = json.loads(path.read_text())["data"]
    df = pd.DataFrame(
        {
            "frame_id": list(map(int, data.keys())),
            "sensing_time_count": [len(v["sensing_time_list"]) for v in data.values()],
            "burst_id_count": [len(v["burst_id_list"]) for v in data.values()],
        }
    )
    df["total_bursts"] = df.sensing_time_count * df.burst_id_count
    df["ministack_count"] = (df.sensing_time_count // 15).astype(int)
    return df


def _rich_summary_table(df: pd.DataFrame, title: str) -> None:
    """Print a neat Rich table of aggregated stats by priority."""
    console = Console()
    table = Table(title=title, header_style="bold magenta")
    table.add_column("Priority")
    table.add_column("# Frames", justify="right")
    table.add_column("Acqs kept", justify="right")
    table.add_column("Acqs lost", justify="right")
    table.add_column("% lost", justify="right")
    table.add_column("Ministacks ≥15", justify="right")

    for _, row in (
        df.groupby("priority")
        .agg(
            frames=("frame_id", "count"),
            kept=("sensing_time_count_selected", "sum"),
            lost=("acqs_lost", "sum"),
            ministacks=("ministack_count_selected", "sum"),
        )
        .assign(pct_lost=lambda x: 100 * x.lost / (x.kept + x.lost))
        .sort_index()
        .reset_index()
    ).iterrows():
        table.add_row(
            row.priority,
            f"{row.frames:>6d}",
            f"{row.kept:>9,d}",
            f"{row.lost:>9,d}",
            f"{row.pct_lost:>6.1f}",
            f"{row.ministacks:>10,d}",
        )
    console.print(table)


if __name__ == "__main__":
    import sys

    path = Path(sys.argv[1])
    file_with_drops = next(path.glob("*consistent-burst-ids-2025*.json"))
    all_file = next(path.glob("*consistent*no-blackout*.json"))

    df_selected = _load_consistent_json(file_with_drops).rename(
        columns=lambda c: f"{c}_selected" if c != "frame_id" else c
    )
    df_all = _load_consistent_json(all_file).rename(
        columns=lambda c: f"{c}_all" if c != "frame_id" else c
    )

    df_cmp = df_all.merge(df_selected, on="frame_id", how="inner")

    # lost acquisitions per frame
    df_cmp["acqs_lost"] = (
        df_cmp.sensing_time_count_all - df_cmp.sensing_time_count_selected
    )
    df_cmp["pct_lost"] = 100 * df_cmp.acqs_lost / df_cmp.sensing_time_count_all

    # add priority & Region-3b flags
    gdf_pri = gpd.read_file(PRIORITY_F)[["frame_id", "priority"]]
    region3b = set(gpd.read_file(REGION3B_F).frame_id.unique())

    df_cmp = df_cmp.merge(gdf_pri, on="frame_id", how="left").fillna(
        {"priority": "unknown"}
    )
    df_cmp.loc[df_cmp.frame_id.isin(region3b), "priority"] = "3b"
    df_cmp.priority = df_cmp.priority.astype(str)

    _rich_summary_table(df_cmp, "DISP-S1 Consistent Burst IDs - after blackout filter")

    sns.set_theme(style="whitegrid")

    # Histogram of sensing-time counts *after filtering*
    plt.figure(figsize=(10, 5))
    bins = np.arange(0, df_cmp.sensing_time_count_selected.max() + 12, 12)
    sns.histplot(
        data=df_cmp,
        x="sensing_time_count_selected",
        hue="priority",
        bins=bins,
        multiple="stack",
        edgecolor=".2",
    )
    plt.xlabel("Sensing-time acquisitions per frame (after blackout)")
    plt.ylabel("# Frames")
    plt.title("Coverage per frame by priority (blackout-filtered)")
    plt.tight_layout()

    plt.figure(figsize=(6, 6))
    sns.histplot(
        data=df_cmp[df_cmp.acqs_lost > 0],
        x="acqs_lost",
        hue="priority",
        bins=bins,
        multiple="stack",
        edgecolor=".2",
    )
    plt.xlabel("Acquisitions lost per frame")
    plt.ylabel("# Frames")
    plt.title("Loss of acquisitions per frame from snow blackouts")
    plt.tight_layout()

    # Scatter: original vs filtered acquisitions
    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=df_cmp[df_cmp.sensing_time_count_selected >= 15],
        x="sensing_time_count_all",
        y="sensing_time_count_selected",
        hue="priority",
        alpha=0.6,
        edgecolor="none",
    )
    plt.plot(
        [0, df_cmp.sensing_time_count_all.max()],
        [0, df_cmp.sensing_time_count_all.max()],
        "--k",
        lw=1,
    )
    plt.xlabel("Acquisitions before blackout")
    plt.ylabel("Acquisitions after blackout")
    plt.title("Loss of acquisitions per frame")
    plt.tight_layout()

    plt.show()

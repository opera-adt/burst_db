"""Module for deciding which months to drop from DISP-S1 frames.

Uses GFS temperature & categorical-snow fields.

The module expects two inputs that you already have saved

* ``ds`` - an ``xarray.Dataset`` holding at least:
    - ``categorical_snow_surface``    - categorical/accumulated snow (0/1) per pixel.
    - ``temperature_2m``      - surface temperature (C).
    - ``time`` dim - 6 hour steps
    - ``latitude``/``longitude`` dims - regular 0.25° GFS grid.

* ``frames_gdf`` - ``geopandas.GeoDataFrame`` from the Frame priority listing.

"""

from __future__ import annotations

from enum import Enum
from typing import Literal

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon
from tqdm.auto import tqdm


def aggregate_weather(
    ds: xr.Dataset,
    *,
    snow_var: str = "categorical_snow_surface",
    temp_var: str = "temperature_2m",
    win: str | int = "1W",
    snow_agg: str = "sum",
) -> xr.Dataset:
    """Aggregate GFS to daily cadence.

    Parameters
    ----------
    ds : xr.Dataset
        Saved GFS dataset.
    snow_var : str
        Name of categorical snow variable.
    snow_var, temp_var : str
        Name of temperature variable to use.
    win : str | int
        Rolling window; e.g. ``"3D"`` or ``72"`` (hours).
    snow_agg : {"sum", "max"}
        Whether fresh-snow should be summed or max-ed across *win*.

    """
    ds_daily_snow = ds[snow_var].resample(time="1D").max()
    ds_daily_temp = ds[temp_var].resample(time="1D").mean()

    if snow_agg == "sum":
        snow_roll = ds_daily_snow.resample(time=win).sum()
    else:
        snow_roll = ds_daily_snow[snow_var].resample(time=win).max()

    tmin = ds_daily_temp.resample(time=win).min()
    tmax = ds_daily_temp.resample(time=win).max()

    out = xr.Dataset({"snow": snow_roll, "tmin": tmin, "tmax": tmax})
    out.snow.attrs["long_name"] = f"days with snow per {win}"
    return out


def bad_month_mask(
    agg: xr.Dataset,
    *,
    snow_threshold: float = 5.0,  # snow weekly sum (days with snow during the week)
    freezing_threshold: float = 1,  # (°C)
    temp_var: Literal["tmin", "tmax"] = "tmax",
    combine: Literal["or", "and"] = "or",
) -> xr.DataArray:
    """Return bool mask (time, latitude, longitude) flagged true for bad period.

    When `combine` is true, the period is bad if *either* condition holds.
    """
    if combine == "or":
        bad = (agg["snow"] >= snow_threshold) | (agg[temp_var] <= freezing_threshold)
    else:
        bad = (agg["snow"] >= snow_threshold) & (agg[temp_var] <= freezing_threshold)
    return bad


def _subset_mask_to_frame(mask: xr.DataArray, poly: Polygon) -> xr.DataArray:
    """Clip the *mask* to pixel centers inside *poly*."""
    mask["latitude"]
    mask["longitude"]
    sub = mask.sel(
        latitude=slice(poly.bounds[3], poly.bounds[1]),
        longitude=slice(poly.bounds[0], poly.bounds[2]),
    )
    yy, xx = np.meshgrid(sub["latitude"].values, sub["longitude"].values, indexing="ij")
    pts = gpd.GeoSeries(
        gpd.points_from_xy(xx.ravel(), yy.ravel()), index=pd.RangeIndex(xx.size)
    )
    inside = pts.within(poly)
    inside_mask = inside.values.reshape(xx.shape)
    return sub.where(inside_mask)


def daily_bad_fraction(mask: xr.DataArray, poly: Polygon) -> pd.Series:
    """Return Series indexed by date with fraction of pixels flagged bad."""
    sub = _subset_mask_to_frame(mask, poly)
    frac = sub.mean(dim=("latitude", "longitude"))
    return frac.to_series()


def plot_frame_timeline(
    frame_id: int,
    mask: xr.DataArray,
    frames_gdf: gpd.GeoDataFrame,
    *,
    mask_fraction: float = 0.5,
    ax: plt.Axes | None = None,
):
    """Plot daily good/bad flags for one frame as a heat-strip & monthly bars."""
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 1))
    poly = frames_gdf.loc[frames_gdf.frame_id == frame_id, "geometry"].iloc[0]
    sub = _subset_mask_to_frame(mask, poly)
    bad_by_day = sub.mean(dim=("latitude", "longitude")) > mask_fraction
    x = pd.to_datetime(bad_by_day["time"].values)
    y = bad_by_day.values.astype(int)
    ax.bar(
        x, np.ones_like(y), width=1, color=["#d62728" if b else "#2ca02c" for b in y]
    )
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title(f"Frame {frame_id}: bad (red) / good (green) periods")
    ax.set_xlim(x.min(), x.max())
    if hasattr(ax, "figure") and hasattr(ax.figure, "tight_layout"):
        ax.figure.tight_layout()
    return ax


class Mode(Enum):
    """Enumeration for blackout summarization modes."""

    CONSERVATIVE = "conservative"
    AGGRESSIVE = "aggressive"
    MEDIAN = "median"


def summarize_blackouts(
    runs: list[tuple[pd.Timestamp, pd.Timestamp]],
    *,
    mode: Mode | str = Mode.MEDIAN,
    pivot_month: int = 7,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Collapse winter *runs* into a single blackout window that respects the year wrap.

    We map every date onto a *pivot* water-year that starts on ``pivot_month``
    (default **July 1**).  Anything before July is treated as *next* year so a
    winter like Nov-2024 → Apr-2025 stays continuous.

    Returns a tuple of two timestamps whose **year is arbitrary** (2000 / 2001)
    but whose order is always start < end.
    """
    if not runs:
        raise ValueError("No blackout runs detected")
    mode = Mode(mode) if isinstance(mode, str) else mode

    def _to_pivot(ts: pd.Timestamp) -> pd.Timestamp:
        year = 2000 + (ts.month < pivot_month)
        return ts.replace(year=year)

    starts = [_to_pivot(s) for s, _ in runs]
    ends = [_to_pivot(e) for _, e in runs]

    if mode == Mode.CONSERVATIVE:
        return min(starts), max(ends)
    if mode == Mode.AGGRESSIVE:
        return max(starts), min(ends)
    if mode == Mode.MEDIAN:

        def med(lst):
            return sorted(lst)[len(lst) // 2]

        return med(starts), med(ends)
    raise ValueError("Invalid mode")


def get_blackout_windows(
    agg: xr.Dataset,
    gdf_priority: gpd.GeoDataFrame,
    snow_threshold: float = 3,
    freezing_threshold: float = -2,
    temp_var: Literal["tmin", "tmax"] = "tmax",
    mask_fraction: float = 0.5,
    frame_id: int | None = None,
    debug: bool = False,
) -> pd.DataFrame:
    """Get blackout windows for frames based on weather conditions.

    Parameters
    ----------
    agg : xr.DataArray
        Aggregated weather data.
    gdf_priority : gpd.GeoDataFrame
        GeoDataFrame containing frame geometries and IDs.
    snow_threshold : float
        Snow threshold for bad conditions.
    freezing_threshold : float
        Temperature threshold for bad conditions.
    temp_var : Literal["tmin", "tmax"]
        Temperature variable to use.
    mask_fraction : float
        Fraction of pixels that must be bad to consider a day bad.
    frame_id : int | None
        Specific frame ID to process, or None for all frames.
    debug : bool
        Whether to print debug information.

    Returns
    -------
    pd.DataFrame
        DataFrame with blackout windows for each frame.

    """
    mask = bad_month_mask(
        agg,
        snow_threshold=snow_threshold,
        freezing_threshold=freezing_threshold,
        temp_var=temp_var,
        combine="or",
    )

    frames = [frame_id] if frame_id is not None else gdf_priority.frame_id.values

    rows = []
    for frame_id in tqdm(frames, desc="Processing frames"):
        poly = gdf_priority.loc[gdf_priority.frame_id == frame_id, "geometry"].iloc[0]
        if _subset_mask_to_frame(mask, poly).size == 0:
            if debug:
                print(f"No mask found for frame {frame_id}; skipping")
            continue

        frac = daily_bad_fraction(mask, poly)  # Series indexed by date
        # runs = blackout_runs(frac, mask_fraction=mask_fraction, max_gap=9)

        start_aggressive = end_aggressive = start_median = end_median = (
            start_conservative
        ) = end_conservative = pd.NaT
        runs = None
        try:
            runs = get_annual_seasons(frac, mask_fraction=mask_fraction)
        except ValueError:
            print(f"No runs found for frame {frame_id}; skipping")
        if runs:
            # runs = runs[1:]  # skip year starting in jan 1
            start_aggressive, end_aggressive = summarize_blackouts(
                runs, mode="aggressive"
            )
            start_median, end_median = summarize_blackouts(runs, mode="median")
            start_conservative, end_conservative = summarize_blackouts(
                runs, mode="conservative"
            )
        rows.append(
            {
                "frame_id": frame_id,
                "mask_fraction": mask_fraction,
                "snow_threshold": snow_threshold,
                "freezing_threshold": freezing_threshold,
                "temp_var": temp_var,
                "start_aggressive": start_aggressive,
                "end_aggressive": end_aggressive,
                "start_median": start_median,
                "end_median": end_median,
                "start_conservative": start_conservative,
                "end_conservative": end_conservative,
            }
        )
    df = pd.DataFrame(rows)
    df["blackout_duration_aggressive"] = np.abs(
        df["start_aggressive"] - df["end_aggressive"]
    ).dt.days
    df["blackout_duration_median"] = np.abs(
        df["start_median"] - df["end_median"]
    ).dt.days
    df["blackout_duration_conservative"] = np.abs(
        df["start_conservative"] - df["end_conservative"]
    ).dt.days
    return df


def get_annual_seasons(
    frac_series: pd.Series,
    *,
    mask_fraction: float = 0.5,
    pivot_month: int = 8,
    min_total_winter_periods: int = 2,
    min_run_len: int = 1,
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Identify the primary freeze-start and thaw-end date for each year.

    This function operates on a "water year" basis (default: July-June)
    to correctly capture winter seasons that span the new year. For each
    water year, it finds the start of the first significant freeze and the
    end of the last significant thaw.

    Parameters
    ----------
    frac_series : pd.Series
        Daily fraction of bad pixels for a single frame.
    mask_fraction : float
        Fraction of bad pixels to consider a day "bad".
    pivot_month : int
        The month starting the "water year" (e.g., 7 for July).
    min_total_winter_periods : int
        The minimum total number of "bad" weeks required in a water year to
        qualify as a winter season.
    min_run_len : int
        A consecutive block of "bad" days must be at least this long
        to be considered a significant freeze event.

    Returns
    -------
    list[tuple[pd.Timestamp, pd.Timestamp]]
        A list of (freeze_start, thaw_end) tuples for each valid winter
        season found in the data.

    """
    # Define the water year for each entry in the time series
    water_year = frac_series.index.map(
        lambda ts: ts.year if ts.month >= pivot_month else ts.year - 1
    )

    seasons = []
    # Group the data by water year and process each one
    for _year_idx, group in frac_series.groupby(water_year):
        # Identify all "bad" days in this water year
        is_bad_series = group >= mask_fraction
        bad_days = group[is_bad_series]

        # 1. Skip years with insufficient winter conditions
        if len(bad_days) < min_total_winter_periods:
            continue

        # 2. Find all consecutive runs of "bad" days
        # Create a grouper for consecutive blocks
        block_grouper = (is_bad_series != is_bad_series.shift()).cumsum()
        bad_blocks = is_bad_series[is_bad_series]

        # Calculate the size (length) of each block of bad days
        block_sizes = bad_blocks.groupby(block_grouper).size()

        # 3. Filter for blocks that are long enough to be significant
        significant_blocks = block_sizes[block_sizes >= min_run_len]
        if significant_blocks.empty:
            continue  # No significant freeze events this year

        # 4. Determine the season's start and end
        # The freeze_start is the beginning of the FIRST significant block
        first_significant_block_id = significant_blocks.index[0]
        freeze_start = bad_blocks[block_grouper == first_significant_block_id].index[0]

        # The thaw_end is the LAST bad day in the entire water year
        thaw_end = bad_days.index.max()

        seasons.append((freeze_start, thaw_end))

    return seasons

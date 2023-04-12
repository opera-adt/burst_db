import math
from collections import Counter, namedtuple
from functools import lru_cache
from itertools import groupby
from typing import List

from itertools import repeat
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from tqdm.auto import tqdm

MIN_FRAME = 5
MAX_FRAME = 12
TARGET_FRAME = 10

FrameSlice = namedtuple("FrameSlice", ["start_idx", "end_idx", "is_land"])


def create_frame_to_burst_mapping(
    is_in_land, target_frame: int, min_frame: int, max_frame: int
) -> pd.DataFrame:
    """Create the JOIN table between frames_number and burst_id."""
    frame_slices = create_frame_slices(is_in_land, min_frame=min_frame)

    # parallelize the solve function calls
    with ProcessPoolExecutor() as executor:
        cumulative_slice_idxs = list(
            tqdm(
                executor.map(
                    _process_slice,
                    frame_slices,
                    repeat(target_frame),
                    repeat(min_frame),
                    repeat(max_frame),
                ),
                desc="Solving frame sizes",
                total=len(frame_slices),
            )
        )

    # Flatten the list of lists, since each result from executor.map is a list
    cumulative_slice_idxs = sorted(
        [item for sublist in cumulative_slice_idxs for item in sublist]
    )

    # Create the frame IDs mapping to burst_id
    # (frame_id, OGC_FID)
    frame_ogc_fid_tuples = []
    for frame_id, (start_idx, end_idx, is_land) in enumerate(
        cumulative_slice_idxs, start=1
    ):
        for burst_id in range(start_idx + 1, end_idx + 1):
            for ogc_fid in range(1 + 3 * (burst_id - 1), 4 + 3 * (burst_id - 1)):
                frame_ogc_fid_tuples.append((frame_id, ogc_fid, is_land))

    df_frame_to_burst_id = pd.DataFrame(
        frame_ogc_fid_tuples, columns=["frame_fid", "burst_ogc_fid", "is_land"]
    )
    return df_frame_to_burst_id


def _process_slice(slice_info, target_frame, min_frame, max_frame):
    start_idx, end_idx, is_land = slice_info
    n = end_idx - start_idx
    cur_slices = solve(
        n,
        target=target_frame,
        min_frame=min_frame,
        max_frame=max_frame,
    )
    # bump up so they refer to rows, instead of being from 0
    return [(s + start_idx, e + start_idx, is_land) for (s, e) in cur_slices]


@lru_cache(maxsize=1000)
def solve(n, target=TARGET_FRAME, max_frame=MAX_FRAME, min_frame=MIN_FRAME):
    """Solve the dynamic programming problem to find the best frame sizes.

    Parameters
    ----------
    n : int
        The number of bursts to split into frames
    target : int
        The target number of bursts per frame
    max_frame : int
        The maximum number of bursts per frame
    min_frame : int
        The minimum number of bursts per frame

    Returns
    -------
    list of tuples
        The start and end indices of each frame

    Notes
    -----
    This is posed as the same problem as the text justification problem.
    "words" == bursts. "lines" == Frames. Where to "break the lines" == "group the frames"

    The differences between here and reference [1]_ are

    1. Our "badness" has both a maximum and minimum size, beyond which the
    badness is infinite.  We also have a "target" so that the majority of
    frames are exactly that size, and only occasionally to we adjust the
    size to be smaller or larger.
    2. The slices are overlapping, so we need to add 1 to the length of each frame.
    This is accounted for in the "badness" function.

    Reference
    ---------
    ..[1] https://ocw.mit.edu/courses/6-006-introduction-to-algorithms-fall-2011/resources/mit6_006f11_lec20/
    """
    # DP[i][0] is the minimum badness of the frames starting at i
    # and DP[i][1] is the index of the next frame (to backtrack and get the slices)
    DP = [None] * (n + 1)
    DP[n] = (0, None)
    for i in range(n - 1, -1, -1):
        DP[i] = min(
            (
                (
                    DP[j][0]
                    + _badness(
                        i, j, target=target, max_frame=max_frame, min_frame=min_frame
                    ),
                    j,
                )
                for j in range(i + 1, n + 1)
            ),
            key=lambda x: x[0],
        )

    # backtrack and get the slices
    slices = []
    i = 0
    while i < n:
        j = DP[i][1]
        # To make the frames overlapping, make it j + 1
        # Just need to account for the final one, which can't be
        # bigger than n.
        end = min(j + 1, n)
        slices.append((i, end))
        i = j

    return slices


def _badness(i, j, target=TARGET_FRAME, max_frame=MAX_FRAME, min_frame=MIN_FRAME):
    """Compute the badness of a frame of length j - i.

    To account for the overlap, each frame will really be 1 bigger (except last)
    so say "max_frame" is 12, even though it'll be 13
    If target is a fraction (like 9.5), than either side of the target is allowed.
    with zero badness
    """
    n = j - i
    # make n+1 cuz it's bigger
    if (n + 1) > max_frame or (n + 1) < min_frame:
        return math.inf
    else:
        return math.floor(abs((n + 1) - target)) ** 3


def create_frame_slices(is_land_indicator, min_frame=MIN_FRAME) -> List[FrameSlice]:
    """Group adjacent frames that are too small."""
    indicator = is_land_indicator.copy()
    ii = 0
    # First iter: make sure all land sequences are at least min_frame long
    for is_land, v in groupby(indicator):
        n_frames = len(list(v))
        ii += n_frames
        if is_land and n_frames < min_frame:
            indicator[ii - min_frame // 2 : ii + min_frame // 2 + 1] = True

    # Second iter: keep looping while there's any small water sequences. make them land
    keep_looping = True
    while keep_looping:
        keep_looping = False
        ii = 0
        for is_land, v in groupby(indicator):
            n_frames = len(list(v))
            ii += n_frames
            if not is_land and n_frames < min_frame:
                keep_looping = True
                indicator[ii - min_frame // 2 : ii + min_frame // 2 + 1] = True
            # loop will break when we didn't adjust any water sequences

    consecutive_land_frames = Counter()
    consecutive_water_frames = Counter()
    frame_slices: List[FrameSlice] = []
    ii, i_prev = 0, 0
    for is_land, cur_indicators in groupby(indicator):
        n_frames = len(list(cur_indicators))
        i_prev = ii
        ii += n_frames
        frame_slices.append(FrameSlice(i_prev, ii, bool(is_land)))

        if is_land:
            consecutive_land_frames[n_frames] += 1
        else:
            consecutive_water_frames[n_frames] += 1

    print("Number of occurrences with consecutive land bursts:")
    print(sorted(consecutive_land_frames.items())[:5], end=",... ")
    print(sorted(consecutive_land_frames.items())[-5:])
    print("Number of occurrences with consecutive water bursts:")
    print(sorted(consecutive_water_frames.items())[:5], end=",... ")
    print(sorted(consecutive_water_frames.items())[-5:])

    return frame_slices

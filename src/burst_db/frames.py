import math
from collections import Counter
from itertools import groupby

MIN_FRAME = 5
MAX_FRAME = 12
TARGET_FRAME = 10


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

    The differences between here and reference ..[1] are

    1. Our "badness" has both a maximum and minumum size, beyond which the
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



def _buffer_small_frames(indicator, min_frame=MIN_FRAME):
    indicator2 = indicator.copy()
    ii = 0
    for k, v in groupby(indicator):
        n_frames = len(list(v))
        ii += n_frames
        if k and n_frames < min_frame:
            indicator2[ii - min_frame // 2 : ii + min_frame // 2 + 1] = True

    consecutive_land_frames = Counter()
    land_slices = []
    ii, i_prev = 0, 0
    for k, v in groupby(indicator2):
        n_frames = len(list(v))
        i_prev = ii
        ii += n_frames
        if not k:
            continue
        land_slices.append((i_prev, ii))
        consecutive_land_frames[n_frames] += 1

    return indicator2, consecutive_land_frames, land_slices


def _make_frame_tuples(land_slices):
    frame_slices = []
    for start_idx, end_idx in land_slices:
        cur_slices = solve(end_idx - start_idx)
        # bump up so they refer to rows, instead of being from 0
        cur_slices = [(s + start_idx, e + start_idx) for (s, e) in cur_slices]
        frame_slices.extend(cur_slices)

    # Create the frame IDs mapping to burst_id
    # (frame_id, OGC_FID)
    frame_tuples = []
    for frame_id, (start_idx, end_idx) in enumerate(frame_slices, start=1):
        for burst_id in range(start_idx + 1, end_idx + 1):
            for ogc_fid in range(1 + 3 * (burst_id - 1), 4 + 3 * (burst_id - 1)):
                frame_tuples.append((frame_id, ogc_fid))
    return frame_tuples

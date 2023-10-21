-- Get the granule where each burst ID was seen for the first time in north america.
-- These are NOT unique granules
SELECT
    burst_id_jpl,
    min(sensing_time) AS first_seen_time,
    granule
FROM
    bursts
WHERE
    is_north_america = 1
GROUP BY
    1;

-- Get the unique granules covering all global burst IDS
-- which were seen after 2017
WITH first_seen AS (
    SELECT
        burst_id_jpl,
        min(sensing_time) AS first_seen_time,
        granule
    FROM
        bursts
    WHERE
        is_north_america = 1
        AND sensing_time > '2017-01-01'
    GROUP BY
        1
)
SELECT
    DISTINCT granule;

-- Make a table grouped by (frame_id, burst_id_jpl) with the
-- count of number of acquisition dates from 2016-06-01 to 2023-06-01
COPY (
    WITH frame_ids AS(
        SELECT
            DISTINCT frame_id
        FROM
            (
                SELECT
                    DISTINCT min_frame_id AS frame_id
                FROM
                    bursts_denormalized
                UNION
                SELECT
                    DISTINCT max_frame_id AS frame_id
                FROM
                    bursts_denormalized
            )
    ),
    unioned_bursts AS (
        SELECT
            min_frame_id AS frame_id,
            burst_id_jpl,
            sensing_time
        FROM
            bursts_denormalized
        UNION
        SELECT
            max_frame_id AS frame_id,
            burst_id_jpl,
            sensing_time
        FROM
            bursts_denormalized
    )
    SELECT
        frame_id,
        burst_id_jpl,
        count(*) AS num_acquisition_dates
    FROM
        unioned_bursts
        JOIN frame_ids USING (frame_id)
    WHERE
        sensing_time BETWEEN '2016-06-01'
        AND '2023-06-01'
    GROUP BY
        ALL
    ORDER BY
        ALL
) TO 'frame_counts.csv' (FORMAT CSV);

-- Now we want, for each frame, to go through the unique
-- acquisition dates and count: how many bursts were seen
-- on this date? (As in, X/27 were present for Frame ID Y on
-- date Z)
-- First, we make an "unstacked" table, where (burst_id, sensing_time)
-- pairs are repeated if they are in different frame IDs
CREATE TABLE bursts_with_frame_ids AS(
    WITH frame_ids AS(
        SELECT
            DISTINCT frame_id
        FROM
            (
                SELECT
                    DISTINCT min_frame_id AS frame_id
                FROM
                    bursts_denormalized
                UNION
                SELECT
                    DISTINCT max_frame_id AS frame_id
                FROM
                    bursts_denormalized
            )
    ),
    unioned_bursts AS (
        SELECT
            min_frame_id AS frame_id,
            burst_id_jpl,
            sensing_time
        FROM
            bursts_denormalized
        UNION
        SELECT
            max_frame_id AS frame_id,
            burst_id_jpl,
            sensing_time
        FROM
            bursts_denormalized
    )
    SELECT
        frame_id,
        burst_id_jpl,
        sensing_time
    FROM
        unioned_bursts
        JOIN frame_ids USING (frame_id)
    GROUP BY
        ALL
    ORDER BY
        ALL
): -- Now we can do the actual counting
COPY (
    SELECT
        frame_id,
        date_trunc('day', sensing_time) AS sensing_date,
        count(*) AS num_bursts
    FROM
        bursts_with_frame_ids
    GROUP BY
        ALL
    ORDER BY
        ALL
) TO 'frame_burst_count_by_date.csv' (FORMAT CSV);

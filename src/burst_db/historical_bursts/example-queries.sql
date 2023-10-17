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

set -e
set -x

in_file=${1-opera-s1-disp.gpkg}
out_file=${2-opera-s1-disp-2d.gpkg}
# Convert geometries to multipolygons
# Convert to 2D (no Z)

if [ -f "$out_file" ]; then
    echo "File exists"
else
    echo "Creating $out_file"

    query_base="SELECT AsGPB(CastToXY(GeomFromGPB(geom))) AS geom, * from"
    # Also simplify to a 0.1 degree tolerance
    # precision=${3-"0.1"}
    # query_base="SELECT AsGPB(CastToXY(st_simplify(GeomFromGPB(geom), $precision))) AS geom, * from"

    ogr2ogr -nln frames -nlt MULTIPOLYGON -dialect sqlite \
        -sql "$query_base frames" \
        "$out_file" "$in_file"
    table_name="burst_id_map"
    ogr2ogr -append -nln $table_name -nlt MULTIPOLYGON -dialect sqlite \
        -sql "$query_base $table_name" \
        "$out_file" "$in_file"
    table_name="frames_bursts"
    ogr2ogr -append -nln $table_name -dialect sqlite \
        -sql "SELECT * from $table_name" \
        "$out_file" "$in_file"
fi

script_dir=$(dirname $(readlink -f $0))

# sql_file="$script_dir/make_geojson.sql"
# duckdb $out_file <$sql_file

sql_query_burst_id=$(
    cat <<'EOF'
SELECT b.fid,
    b.burst_id_jpl,
    is_land,
    is_north_america,
    b.orbit_pass,
    st_simplify(b.geom, 0.1) AS geometry
FROM frames f
    JOIN frames_bursts fb ON fb.frame_fid = f.fid
    JOIN burst_id_map b ON fb.burst_ogc_fid = b.fid
WHERE is_land=1;
EOF
)

sql_query_frame=$(
    cat <<'EOF'
SELECT f.fid,
    is_land,
    is_north_america,
    orbit_pass,
    st_simplify(f.geom, 0.1) AS geometry
FROM frames f
WHERE is_land=1;
EOF
)

ogr2ogr -f GeoJSON -dialect sqlite -sql "$sql_query_frame" \
    "frame_geometries_simple.geojson" $out_file
zip frame_geometries_simple.geojson.zip frame_geometries_simple.geojson

ogr2ogr -f GeoJSON -dialect sqlite -sql "$sql_query_burst_id" \
    "burst_id_geometries_simple.geojson" $out_file
zip burst_id_geometries_simple.geojson.zip burst_id_geometries_simple.geojson

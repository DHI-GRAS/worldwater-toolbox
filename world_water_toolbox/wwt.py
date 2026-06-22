# -*- coding: utf-8 -*-
"""
Created on Jun 2023
@authors: Andrea Sulova, Walid Ghariani, Alkis Kokkalis - DHI GRAS

Updated on 28-05-2026 by Walid Ghariani:
- Compatibility with OpenEO UDFs: https://github.com/DHI-GRAS/worldwater-toolbox/issues/11
- Some clean up and code refactoring
- def added a world water indicators udp 
"""

from typing import Union
import os
from dateutil.relativedelta import relativedelta

import geopandas as gpd
import openeo
from openeo import DataCube
from openeo.api.process import Parameter
from openeo.extra.spectral_indices.spectral_indices import (
    append_indices,
    compute_indices,
)
from openeo.processes import if_, exp, log, eq, date_shift
from openeo.rest.udp import build_process_dict

from .openeo_filters import (
    FILTER_TYPES,
    HYSTERESIS_DEFAULTS,
    apply_hysteresis,
    apply_spatial_filter,
)


# Collections
COLLECTION_S2 = "SENTINEL2_L2A"
COLLECTION_S1 = "SENTINEL1_GRD"
COLLECTION_DEM = "COPERNICUS_30"
COLLECTION_WORLDCOVER = "ESA_WORLDCOVER_10M_2020_V1"

# DEM temporal extent
DEM_DATE_START = "2010-01-01"
DEM_DATE_END = "2030-12-31"

# OpenEO job executor resource options
JOB_EXECUTOR_MEMORY = "1g"
JOB_EXECUTOR_MEMORY_OVERHEAD = "3g"
JOB_EXECUTOR_CORES = 1

# Neighborhood processing window size and time unit
NEIGHBORHOOD_WINDOW_PX = 256
NEIGHBORHOOD_TIME_UNIT = "P1D"

# CLP cloud probability threshold 
CLP_CLOUD_THRESHOLD = 0.3

# S2 water probability threshold applied inside s2_water_processing
# not used for now, as the S2 water probability is not thresholded but 
# merged with S1 using the LOOKUPTABLE-based fusion method
#S2_WATER_THRESHOLD = 0.75

# SAR backscatter processing options
SAR_IMPLEMENTATION_VERSION = "2"
SAR_TILE_SIZE = 256
SAR_OTB_MEMORY = 1024

# UDF file pth: resolved relative to this module
UDF_FILE_PATH = os.path.join(os.path.dirname(__file__), "udf.py")


LOOKUPTABLE = {
    "Deserts": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-7.03 + (-0.44 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.133 + (-5.92 * ndvi) + (14.82 * ndwi))))
        ),
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-3.69 + (-0.25 * vv) + (0.47 * ndvi) + (15.3 * ndwi))))
        ),
    },
    "Mountain": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-3.76 + (-0.262 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.262 + (0.75 * ndvi) + (12.65 * ndwi))))
        ),
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-1.13 + (-0.11 * vv) + (3.03 * ndvi) + (13.21 * ndwi))))
        ),
    },
    "Tropical forest": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-5.8 + (-0.415 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.344 + (2.886 * ndvi) + (11.91 * ndwi))))
        )
        * 100,
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-3.25 + (-0.23 * vv) + (4.17 * ndvi) + (9.5 * ndwi))))
        ),
    },
    "Tropical savanna": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-7.0 + (-0.444 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.344 + (2.886 * ndvi) + (11.91 * ndwi))))
        )
        * 100,
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-1.06 + (-0.17 * vv) + (3.82 * ndvi) + (14.4 * ndwi))))
        ),
    },
    "Subtropical savanna": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-7.17 + (-0.48 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.845 + (2.14 * ndvi) + (13.5 * ndwi))))
        ),
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-2.64 + (-0.23 * vv) + (8.6 * ndwi))))
        ),
    },
    "Subtropical forest": {
        "S1": lambda vh, vv: (1 / (1 + exp(-(-6.67 + (-6.67 * vv))))),
        "S2": lambda ndvi, ndwi: (
            1 / (1 + exp(-(0.712 + (-1.133 * ndvi) + (7.16 * ndwi))))
        ),
        "S1_S2": lambda vv, ndvi, ndwi: (
            1 / (1 + exp(-(-2.72 + (-0.22 * vv) + (-0.49 * ndvi) + (4.55 * ndwi))))
        ),
    },
    "Temperate broadleaf": {
        "S1": lambda vh, vv: 1 / (1 + exp(-(-8.82 + (-0.58 * vv)))),
        "S2": lambda ndvi, ndwi: 1
        / (1 + exp(-(-0.013 + (5.38 * ndvi) + (13.79 * ndwi)))),
        "S1_S2": lambda vv, ndvi, ndwi: 1
        / (1 + exp(-(-2.7 + (-0.2 * vv) + (3.6 * ndvi) + (9.73 * ndwi)))),
    },
    "Temperate grassland": {
        "S1": lambda vh, vv: 1 / (1 + exp(-(-7.01 + (-0.426 * vv)))),
        "S2": lambda ndvi, ndwi: 1
        / (1 + exp(-(1.286 + (8.74 * ndvi) + (23.217 * ndwi)))),
        "S1_S2": lambda vv, ndvi, ndwi: 1
        / (1 + exp(-(-3.43 + (-0.25 * vv) + (11.74 * ndvi) + (22.035 * ndwi)))),
    },
    "Tundra": {
        "S1": lambda vh, vv: 1 / (1 + exp(-(-4.57 + (-0.13 * vh) + (-0.11 * vv)))),
        "S2": lambda ndvi, ndwi: 1 / (1 + exp(-(1.07 + (3.6 * ndwi)))),
        #             "S1_S2": lambda hh, hv, ndwi: 1 / (1 + exp(- (-3.8 + (-0.15 * hh) + (-0.03 * hv) + (2.5 * ndwi))))
        "S1_S2": lambda vv, ndvi, ndwi: 1
        / (1 + exp(-(-3.8 + (-0.15 * vv) + (-0.03 * ndvi) + (2.5 * ndwi)))),
    },
}


def hillshade_mask(
    connection: openeo.Connection,
    spatial_extent: dict,
    start_date_exclusion,
    month_end,
    cloud_cover: int,
) -> DataCube:

    # Hill-shade masking
    collection = COLLECTION_S2

    # angles are lower resolution, so we load them separately
    s2_angles = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=["B01", "sunAzimuthAngles", "sunZenithAngles"],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover},
    )

    dem_cube = connection.load_collection(
        COLLECTION_DEM,
        spatial_extent=spatial_extent,
        temporal_extent=[DEM_DATE_START, DEM_DATE_END],
    )
    dem_cube = dem_cube.max_time()
    # Resample s2 cube (Azimuth and Zenith) to 30m
    s2_cube_30 = s2_angles.resample_spatial(resolution=30, method="average")
    # DEM WGS84 to DEM UTM from S2 cube
    # dem_cube_s2 = dem_cube.resample_cube_spatial(s2_cube_30)
    # Merge 30m DEM and 30m s2 cube due Azimuth and Zenith
    merged_cube = s2_cube_30.merge_cubes(dem_cube)
    # Apply the hill-shade udf
    #     print('udf', os.path.dirname(os.path.abspath(__file__)) + '/udf.py')
    process = openeo.UDF.from_file(UDF_FILE_PATH, runtime="Python")
    hillshade = merged_cube.apply_neighborhood(
        process=process,
        size=[
            {"dimension": "t", "value": NEIGHBORHOOD_TIME_UNIT},
            {"dimension": "x", "unit": "px", "value": NEIGHBORHOOD_WINDOW_PX},
            {"dimension": "y", "unit": "px", "value": NEIGHBORHOOD_WINDOW_PX},
        ],
        overlap=[
            {"dimension": "x", "unit": "px", "value": "8"},
            {"dimension": "y", "unit": "px", "value": "8"},
        ],
    )
    # Rename bands in a hill shade cube
    hillshade = hillshade.rename_labels(
        "bands", ["hillshade_mask", "sunAzimuthAngles", "sunZenithAngles", "DEM"]
    )
    return hillshade


def masked_s2_cube(
    connection: openeo.Connection,
    spatial_extent: dict,
    start_date_exclusion,
    month_end,
    cloud_cover: int,
) -> DataCube:

    collection = COLLECTION_S2

    col_metadata = connection.describe_collection(collection)
    bands = col_metadata.get("cube:dimensions", {}).get("bands", {}).get("values", [])

    # Loading S2 collection
    s2_cube = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=["B02", "B03", "B04", "B08"],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover},
    )
    # Loading S2 collection for cloud masking processing
    if "CLP" in bands:
        cloudmask_bands = ["CLP", "SCL"]
    else:
        cloudmask_bands = ["SCL"]

    s2_cube_masking = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=cloudmask_bands,
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover},
    )
    # Applying mask from sen2cor
    scl = s2_cube_masking.band("SCL")
    mask_scl = (scl == 3) | (scl == 8) | (scl == 9) | (scl == 10) | (scl == 11) 
    if "CLP" in bands:
        clp = s2_cube_masking.band("CLP")
        mask_clp = mask_scl | (clp / 255) > CLP_CLOUD_THRESHOLD
    else:
        mask_clp = mask_scl

    if "sunAzimuthAngles" in bands and "sunZenithAngles" in bands:
        hillshade = (
            hillshade_mask(
                connection, spatial_extent, start_date_exclusion, month_end, cloud_cover
            )
            .filter_bands("hillshade_mask")
            .linear_scale_range(0, 10, 0, 10)
        )
        the_mask = mask_clp.merge_cubes(
            hillshade.drop_dimension("bands"), overlap_resolver="sum"
        )
    else:
        the_mask = mask_clp
    # it is better to create one combined mask, which can be used to reduce data loading in full resolution cube
    # Mask s2 cube with a hillshade mask
    s2_cube_masked = s2_cube.mask(the_mask)
    # Replace 0 to nan in s2 cubes
    s2_cube = s2_cube_masked.apply(lambda x: if_(x.neq(0), x))
    return s2_cube


def _water_extent_multiple_months(
    connection: openeo.Connection,
    month_start,
    month_end,
    geometry: str,
    region: Union[str, Parameter],
    cloud_cover: Union[int, Parameter],
    threshold: Union[int, float, Parameter] = 75,
    only_s1: Union[bool, Parameter] = False,
    filter_type: Union[str, Parameter] = "none",
    gaussian_sigma: float = 0.85,
    uniform_size: int = 3,
    median_size: int = 3,
    use_hysteresis: Union[bool, Parameter] = False,
    high_thr: Union[float, Parameter] = HYSTERESIS_DEFAULTS["high_thr"],
    low_thr: Union[float, Parameter] = HYSTERESIS_DEFAULTS["low_thr"],
    max_dist_m: Union[float, Parameter] = HYSTERESIS_DEFAULTS["max_dist_m"],
    connectivity: Union[int, Parameter] = HYSTERESIS_DEFAULTS["connectivity"],
) -> DataCube:
    spatial_extent = geometry if isinstance(geometry, Parameter) else _get_spatial_extent(geometry)
    start_date_exclusion = date_shift(month_start, value=-1, unit="month")

    _skip_s2 = isinstance(only_s1, bool) and only_s1
    if not _skip_s2:
        s2_cube = masked_s2_cube(connection, spatial_extent, start_date_exclusion, month_end, cloud_cover)
        s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube, region)
        s2_median_water = s2_cube_water.filter_temporal([month_start, month_end]).aggregate_temporal_period("month", "median")
        ndxi_median = ndxi_cube.filter_temporal([month_start, month_end]).aggregate_temporal_period("month", "median")
    else:
        s2_median_water = None
        ndxi_median = None

    s1_cube = sentinel1_preprocessing(connection, month_end, month_start, spatial_extent, region)
    s1_median = s1_cube.aggregate_temporal_period("month", "median")

    water_prob = _water_probability(s1_median, s2_median_water, ndxi_median, region, only_s1)

    builtup_mask_cube = None
    if COLLECTION_WORLDCOVER in connection.list_collection_ids():
        worldcover_cube = connection.load_collection(
            COLLECTION_WORLDCOVER, spatial_extent=spatial_extent, bands=["MAP"]
        )
        builtup_mask_cube = (
            (worldcover_cube.band("MAP") == 50)
            .max_time()
            .resample_cube_spatial(water_prob)
        )
        water_probability = water_prob.mask(builtup_mask_cube)
    else:
        water_probability = water_prob
    water_probability = water_probability.rename_labels("bands", ["water_prob_sum"])

    smoothed = apply_spatial_filter(
        water_probability,
        filter_type,
        gaussian_sigma=gaussian_sigma,
        uniform_size=uniform_size,
        median_size=median_size,
    )
    smoothed = smoothed.rename_labels("bands", ["water_prob_sum"])

    if builtup_mask_cube is not None:
        smoothed = smoothed.mask(builtup_mask_cube)
        smoothed = smoothed.rename_labels("bands", ["water_prob_sum"])

    water_mask_simple = smoothed.apply(
        lambda x: x["water_prob_sum"] * 100 > threshold
    ).rename_labels("bands", ["surface_water"])

    water_mask_hyst = apply_hysteresis(
        smoothed.filter_bands(["water_prob_sum"]),
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    ).rename_labels("bands", ["surface_water"])

    if isinstance(use_hysteresis, bool):
        water_mask = water_mask_hyst if use_hysteresis else water_mask_simple
    else:
        water_mask_raw = if_(use_hysteresis, water_mask_hyst, water_mask_simple)
        water_mask = DataCube(
            water_mask_raw.pgnode,
            water_mask_hyst.connection,
            metadata=water_mask_hyst.metadata,
        )

    output = smoothed.merge_cubes(water_mask)
    output = output * 1.0

    if builtup_mask_cube is not None:
        output = output.mask(builtup_mask_cube)

    return output


def _monthly_water_probability(
    connection: openeo.Connection,
    spatial_extent: dict,
    start_date,
    end_date,
    cloud_cover: Union[int, Parameter],
    region: Union[str, Parameter],
    only_s1: Union[bool, Parameter] = False,
) -> DataCube:
    """Monthly water probability for a multi-month period.

    Returns a DataCube with a time dimension (one step per month) and a single
    band holding the fused S1/S2 water probability in ~[0, 1].
    """
    start_date_exclusion = date_shift(start_date, value=-1, unit="month")

    s1_cube = sentinel1_preprocessing(connection, end_date, start_date, spatial_extent, region)
    s1_median = s1_cube.aggregate_temporal_period("month", "median")

    _skip_s2 = isinstance(only_s1, bool) and only_s1
    if not _skip_s2:
        s2_cube = masked_s2_cube(connection, spatial_extent, start_date_exclusion, end_date, cloud_cover)
        _, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube, region)
        s2_median_water = s2_cube_water.filter_temporal([start_date, end_date]).aggregate_temporal_period("month", "median")
        ndxi_median = ndxi_cube.filter_temporal([start_date, end_date]).aggregate_temporal_period("month", "median")
    else:
        s2_median_water = None
        ndxi_median = None

    return _water_probability(s1_median, s2_median_water, ndxi_median, region, only_s1)


def _compute_water_indicators(
    binary_monthly: DataCube,
) -> DataCube:
    """Derive annual water indicators from monthly binary water masks.

    Returns a 5-band DataCube (no time dimension):
        mean_swf, annual_classification, min_water_extent,
        max_water_extent, water_seasonality
    """

    # Count water months per pixel
    water_count = binary_monthly.reduce_dimension(
        reducer=lambda data: data.sum(),
        dimension="t",
    ).rename_labels("bands", ["water_count"])

    # Mean surface water frequency
    mean_swf = binary_monthly.mean_time().rename_labels(
        "bands", ["mean_swf"]
    )
    # Annual classification: 
    # 1 = Non-water (0-1 months) 
    # 2 = Seasonal water (2-9 months) 
    # 3 = Permanent water (10-12 months)
    annual_class = water_count.apply(
        lambda x: if_(
            x["water_count"] >= 10,
            3,
            if_(
                x["water_count"] >= 2,
                2,
                1,
            ),
        )
    ).rename_labels("bands", ["annual_classification"])

    # Minimum water extent
    # 1 = Non-water (0-9 months)
    # 2 = Minimum water (10-12 months)
    min_water = water_count.apply(
        lambda x: if_(
            x["water_count"] >= 10,
            2,
            1,
        )
    ).rename_labels("bands", ["min_water_extent"])

    # Maximum water extent
    # 1 = Non-water (0-1 months)
    # 2 = Maximum water (2-12 months)
    max_water = water_count.apply(
        lambda x: if_(
            x["water_count"] >= 2,
            2,
            1,
        )
    ).rename_labels("bands", ["max_water_extent"])

    # Water seasonality
    # 1 = 0-1 months
    # 2 = 2-3 months
    # 3 = 4-6 months
    # 4 = 7-9 months
    # 5 = 10-12 months
    seasonality = water_count.apply(
        lambda x: if_(
            x["water_count"] >= 10,
            5,
            if_(
                x["water_count"] >= 7,
                4,
                if_(
                    x["water_count"] >= 4,
                    3,
                    if_(
                        x["water_count"] >= 2,
                        2,
                        1,
                    ),
                ),
            ),
        )
    ).rename_labels("bands", ["water_seasonality"])


    return (
        mean_swf
        .merge_cubes(annual_class)
        .merge_cubes(min_water)
        .merge_cubes(max_water)
        .merge_cubes(seasonality)
    )


def _water_extent(
    connection: openeo.Connection,
    month_start,
    month_end,
    start,
    end,
    geometry: str,
    region: str,
    threshold: Union[int, float],
    cloud_cover: int,
    rgb_processing: bool,
    only_s1: bool = False,
    output_name_o: str = "",
    filter_type: str = "none",
    gaussian_sigma: float = 0.85,
    uniform_size: int = 3,
    median_size: int = 3,
    use_hysteresis: bool = False,
    high_thr: float = HYSTERESIS_DEFAULTS["high_thr"],
    low_thr: float  = HYSTERESIS_DEFAULTS["low_thr"],
    max_dist_m: float = HYSTERESIS_DEFAULTS["max_dist_m"],
    connectivity: int = HYSTERESIS_DEFAULTS["connectivity"],
) -> str:
    """
    Calculate water extent (binary mask) using S1 and S2 collection.

    Parameters
    ----------
    connection : string
        Connection to manage and persist settings when interacting with the OpenEO API.

    month_start : date
        Start date of current processing.

    month_end : date
        End date of current processing.

    start : date
        The start date for the entire period that will be processed.

    end : date
         The end date for the entire period that will be processed.

    geometry : Geojson
        Geojson file path to the area of interest (AOI).

    region : string
        Choose an eco-region from the options.

    threshold : integer
        Custom threshold used for the water probability layer.

    cloud_cover : integer
        Maximum cloud cover allowed for the Sentinel-2 collection.


    Returns
    -------
    output_folder : string
         Path to the storage of water masks.

    """

    spatial_extent = _get_spatial_extent(geometry)
    output, s2_cube = _water_extent_for_month(
        connection,
        spatial_extent,
        region,
        month_start,
        month_end,
        cloud_cover,
        threshold,
        only_s1,
        filter_type=filter_type,
        gaussian_sigma=gaussian_sigma,
        uniform_size=uniform_size,
        median_size=median_size,
        use_hysteresis=use_hysteresis,
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    )
    # Output folder
    region_naming = "_".join(region.split(" "))
    output_folder = (
        output_name_o
        + "/"
        + f'WWT_{region_naming}_{start.strftime("%Y_%m")}_{end.strftime("%Y_%m")}'
    )
    prefix = "water_onlyS1" if only_s1 else "water_"
    month_start_str = month_start.strftime("%Y_%m")
    filename = prefix + "_" + month_start_str + "_" + month_start_str
    # Send jobs to the backend
    print(
        "Water mask processing between: "
        + month_start_str
        + " and "
        + str(month_end.strftime("%Y_%m"))
    )
    my_job = output.create_job(
        title=prefix + month_start_str + "_" + str(month_end.strftime("%Y_%m")),
        out_format="GTiff",
        job_options={
            "executor-memory": JOB_EXECUTOR_MEMORY,
            "executor-memoryOverhead": JOB_EXECUTOR_MEMORY_OVERHEAD,
            "executor-cores": JOB_EXECUTOR_CORES,
        },
        **{"filename_prefix": filename},
    )
    results = my_job.start_and_wait().get_results()
    results.download_files(output_folder)

    if rgb_processing and s2_cube is not None:
        # Process S2 median image if rgb_processing is True
        print(
            "S2 median image processing between: "
            + month_start_str
            + " and "
            + str(month_end.strftime("%Y_%m"))
        )
        s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time()
        my_job = s2_cube_median.create_job(
            title="median_" + month_start_str + "_" + str(month_end.strftime("%Y_%m")),
            out_format="GTiff",
            job_options={"node_caching": True},
            **{
                "filename_prefix": "median_"
                + month_start_str
                + "_"
                + str(month_end.strftime("%Y_%m"))
            },
        )
        results = my_job.start_and_wait().get_results()
        results.download_files(output_folder)

    return output_folder


def _water_extent_for_month(
    connection: openeo.Connection,
    spatial_extent: dict,
    region: Union[str, Parameter],
    month_start,
    month_end,
    cloud_cover: Union[int, Parameter],
    threshold: Union[int, float, Parameter],
    only_s1: Union[bool, Parameter] = False,
    filter_type: Union[str, Parameter] = "none",
    gaussian_sigma: float = 0.85,
    uniform_size: int = 3,
    median_size: int = 3,
    use_hysteresis: Union[bool, Parameter] = False,
    high_thr: Union[float, Parameter] = HYSTERESIS_DEFAULTS["high_thr"],
    low_thr: Union[float, Parameter]  = HYSTERESIS_DEFAULTS["low_thr"],
    max_dist_m: Union[float, Parameter] = HYSTERESIS_DEFAULTS["max_dist_m"],
    connectivity: Union[int, Parameter] = HYSTERESIS_DEFAULTS["connectivity"],
):
    month_start = (
        str(month_start) if not isinstance(month_start, Parameter) else month_start
    )
    start_date_exclusion = date_shift(month_start, value=-1, unit="month")

    # Normalized radar back-scatter
    s1_cube = sentinel1_preprocessing(
        connection, month_end, month_start, spatial_extent, region
    )
    # Calculate S1 median mosaic
    s1_median = s1_cube.median_time()

    # Skip S2 when only_s1 is a Python True at call time.
    # When only_s1 is an OpenEO Parameter (UDP export) the full graph must be built.
    _skip_s2 = isinstance(only_s1, bool) and only_s1
    if not _skip_s2:
        s2_cube = masked_s2_cube(
            connection, spatial_extent, start_date_exclusion, month_end, cloud_cover
        )
        s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube, region)
        s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time()
        s2_median_water = s2_cube_water.filter_temporal(
            [month_start, month_end]
        ).median_time()
        ndxi_median = ndxi_cube.filter_temporal([month_start, month_end]).median_time()
    else:
        s2_cube_median = None
        s2_median_water = None
        ndxi_median = None
    merge_all = _water_probability(
        s1_median, s2_median_water, ndxi_median, region, only_s1
    )

    builtup_mask_cube = None
    if COLLECTION_WORLDCOVER in connection.list_collection_ids():
        # Mask built-up area using ESA world cover layer
        worldcover_cube = connection.load_collection(
            COLLECTION_WORLDCOVER, spatial_extent=spatial_extent, bands=["MAP"]
        )
        builtup_mask_cube = (
            (worldcover_cube.band("MAP") == 50)
            .max_time()
            .resample_cube_spatial(merge_all)
        )
        water_probability = merge_all.mask(builtup_mask_cube)
    else:
        water_probability = merge_all
    water_probability = water_probability.rename_labels("bands", ["water_prob_sum"])

    # optional spatial smoothing before thresholding
    smoothed = apply_spatial_filter(
        water_probability,
        filter_type,
        gaussian_sigma=gaussian_sigma,
        uniform_size=uniform_size,
        median_size=median_size,
    )
    smoothed = smoothed.rename_labels("bands", ["water_prob_sum"])

    # Re-apply built-up mask after spatial filter: filters can bleed-fill values
    # into masked pixels, so we must exclude built-up pixels again before thresholding.
    if builtup_mask_cube is not None:
        smoothed = smoothed.mask(builtup_mask_cube)
        smoothed = smoothed.rename_labels("bands", ["water_prob_sum"])

    # Simple threshold mask (water_prob_sum is [0,1]; threshold is [0,100])
    water_mask_simple = smoothed.apply(lambda x: x["water_prob_sum"] * 100 > threshold)
    water_mask_simple = water_mask_simple.rename_labels("bands", ["surface_water"])

    # Hysteresis mask (thresholds are [0,1], matching the raw probability scale)
    water_mask_hyst = apply_hysteresis(
        smoothed.filter_bands(["water_prob_sum"]),
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    )
    water_mask_hyst = water_mask_hyst.rename_labels("bands", ["surface_water"])

    # Select mask
    if isinstance(use_hysteresis, bool):
        water_mask = water_mask_hyst if use_hysteresis else water_mask_simple
    else:
        water_mask_raw = if_(use_hysteresis, water_mask_hyst, water_mask_simple)
        water_mask = DataCube(
            water_mask_raw.pgnode,
            water_mask_hyst.connection,
            metadata=water_mask_hyst.metadata,
        )

    output = smoothed.merge_cubes(water_mask)
    output = output * 1.0

    # Final built-up mask applied to both output bands (water_prob_sum + surface_water).
    if builtup_mask_cube is not None:
        output = output.mask(builtup_mask_cube)

    return output, s2_cube_median


def _water_probability(
    s1_median: DataCube,
    s2_median_water: DataCube,
    ndxi_median: DataCube,
    region: Union[str, Parameter],
    only_s1: Union[bool, Parameter],
) -> DataCube:
    def log_(x):
        return 10 * log(x, 10)

    s1_median = s1_median.apply(log_)
    is_tundra = eq(region, "Tundra")

    # S1-only water probability - Tundra path: VH=data[0], VV=data[1]
    s1_tundra_water = 0
    for key in LOOKUPTABLE.keys():
        s1_tundra_water = if_(
            eq(region, key),
            s1_median.filter_bands(["VH", "VV"]).reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S1"](vh=data[0], vv=data[1]),
                dimension="bands",
            ),
            s1_tundra_water,
        )

    # S1-only water probability - non-Tundra path: VV only, passed as both vh and vv args
    s1_nontundra_water = 0
    for key in LOOKUPTABLE.keys():
        s1_nontundra_water = if_(
            eq(region, key),
            s1_median.filter_bands(["VV"]).reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S1"](vh=data[0], vv=data[0]),
                dimension="bands",
            ),
            s1_nontundra_water,
        )

    s1_median_water_raw = if_(is_tundra, s1_tundra_water, s1_nontundra_water)
    s1_median_water = DataCube(
        s1_median_water_raw.pgnode,
        s1_median.connection,
        metadata=s1_median.metadata.reduce_dimension("bands"),
    )

    # S1-only result (selected at openEO runtime when only_s1=True)
    s1_only_result = s1_median_water.add_dimension("bands", "var", type="bands")

    # When only_s1 is a Python True (not a UDP Parameter), skip S1+S2 and S2 graph construction entirely.
    if isinstance(only_s1, bool) and only_s1:
        return s1_only_result

    # S1+S2 fusion - single path for all regions (Tundra S1_S2 lambda uses vv+ndvi+ndwi, not HH/HV)
    # After merge_cubes: data[0]=VV, data[1]=NDWI, data[2]=NDVI
    s1_s2_cube = s1_median.filter_bands(["VV"]).merge_cubes(ndxi_median)
    s1_s2_water_raw = 0
    for key in LOOKUPTABLE.keys():
        s1_s2_water_raw = if_(
            eq(region, key),
            s1_s2_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S1_S2"](
                    vv=data[0], ndvi=data[2], ndwi=data[1]
                ),
                dimension="bands",
            ),
            s1_s2_water_raw,
        )

    s1_s2_water = DataCube(
        s1_s2_water_raw.pgnode,
        s1_median.connection,
        metadata=s1_s2_cube.metadata.reduce_dimension("bands"),
    )

    merged = (
        s1_s2_water.add_dimension("bands", "s1_s2_water", type="bands")
        .merge_cubes(s2_median_water)
        .merge_cubes(s1_median_water)
    )

    def combined(bands):
        s1_s2_water = bands.array_element(0)
        s1_s2_mask = s1_s2_water >= 0

        s2_median_water = bands.array_element(1)
        s2_mask = if_(s1_s2_mask == 0, s2_median_water) >= 0

        s1_median_water = bands.array_element(2)
        s1_mask = if_(s2_mask == 0, if_(s1_s2_mask == 0, s1_median_water)) >= 0

        s1_s2_masked = if_(s1_s2_mask != 0, s1_s2_water, 0)
        s2_masked = if_(s2_mask != 0, s2_median_water, 0)
        s1_masked = if_(s1_mask != 0, s1_median_water, 0)

        return s1_s2_masked + s2_masked + s1_masked

    merge_all = merged.apply_dimension(combined, dimension="bands")

    # Select S1-only or full fusion at openEO runtime - if_() works for both Python bools and UDP Parameters
    result_raw = if_(only_s1, s1_only_result, merge_all)
    return DataCube(
        result_raw.pgnode, s1_median.connection, metadata=s1_only_result.metadata
    )


def sentinel1_preprocessing(
    connection: openeo.Connection,
    month_end,
    month_start,
    spatial_extent: dict,
    region: Union[str, Parameter],
) -> DataCube:
    # Load VH+VV: universally available in IW mode (dual-pol VV+VH).
    # Tundra S1-only path uses VH as cross-pol proxy; all other paths use VV only.
    s1_cube = connection.load_collection(
        COLLECTION_S1,
        spatial_extent=spatial_extent,
        temporal_extent=[month_start, month_end],
        bands=["VH", "VV"],
    )
    s1_cube = s1_cube.sar_backscatter(
        coefficient="sigma0-ellipsoid",
        mask=False,
        elevation_model=COLLECTION_DEM,
        options={
            "implementation_version": SAR_IMPLEMENTATION_VERSION,
            "tile_size": SAR_TILE_SIZE,
            "otb_memory": SAR_OTB_MEMORY,
            "debug": True,
        },
    )
    return s1_cube


def _get_spatial_extent(geometry: str) -> dict:
    gdf = gpd.read_file(geometry)
    gdf = gdf.to_crs(4326)
    bbox = gdf.geometry.total_bounds
    spatial_extent = {
        "west": bbox[0],
        "east": bbox[2],
        "south": bbox[1],
        "north": bbox[3],
        "crs": 4326,
    }
    return spatial_extent


def s2_water_processing(s2_cube: DataCube, region: Union[str, Parameter]) -> tuple:
    # Always compute both indices for UDP compatibility (Tundra S2 lambda accepts but ignores ndvi)
    # ndxi_cube bands: data[0]=NDWI, data[1]=NDVI
    s2_indices = ["NDWI", "NDVI"]
    s2_cube = append_indices(s2_cube, s2_indices)
    ndxi_cube = compute_indices(s2_cube, s2_indices)
    s2_cube = s2_cube.rename_labels("bands", ["B02", "B03", "B04", "B08"] + s2_indices)

    # Create water extent for the S2 collection using logistic expressions from Lookup Table
    s2_cube_water = 0
    for key in LOOKUPTABLE.keys():
        # the if/else lookup is done in the process graph, for the UDP, and can not yet be done in the callback
        s2_cube_water = if_(
            eq(region, key),
            ndxi_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S2"](ndwi=data[0], ndvi=data[1]),
                dimension="bands",
            ),
            s2_cube_water,
        )
    s2_cube_water = DataCube(
        s2_cube_water.pgnode,
        ndxi_cube.connection,
        metadata=ndxi_cube.metadata.reduce_dimension("bands"),
    )
    s2_cube_water = s2_cube_water.add_dimension("bands", "water_prob", type="bands")
    # # Generate a binary water mask by applying a threshold to the water probability band
    # s2_cube_water_threshold = s2_cube_water.apply_dimension(
    #     dimension="bands", process=lambda x: if_(x > S2_WATER_THRESHOLD, x, 0)
    # )
    # s2_cube_water_threshold = s2_cube_water_threshold.rename_labels("bands", ["w_T75"])
    return s2_cube, ndxi_cube, s2_cube_water


def generate_water_extent_udp(connection: openeo.Connection):
    """ water extent UDP for a given month"""
    start_date = Parameter.date(name="start_date", description="The start date.")
    spatial_extent = Parameter.spatial_extent()

    region = Parameter.string(
        name="region",
        description="Eco-Region on which to compute water probability",
        default="Deserts",
        values=list(LOOKUPTABLE.keys()),
    )
    only_s1 = Parameter.boolean(
        name="only_s1",
        description="Boolean variable to specify if only Sentinel-1 data will be used (True) or both Sentinel-1 and Sentinel-2",
        default=False,
    )
    rgb_processing = Parameter.boolean(
        name="rgb_processing",
        description="Boolean variable to specify if Sentinel-2 rgb image will be generated",
        default=False,
    )
    cloud_cover = Parameter.number(
        name="cloud_cover", description="Maximum cloud cover percentage", default=85
    )
    water_threshold = Parameter.number(
        name="water_threshold",
        description="Water probability threshold (0-100). Used when use_hysteresis=false.",
        default=75,
    )
    filter_type = Parameter.string(
        name="filter_type",
        description=(
            "Spatial smoothing filter applied to the water probability before thresholding. "
            f"One of {FILTER_TYPES}. Use 'none' to skip smoothing."
        ),
        default="none",
        values=FILTER_TYPES,
    )
    use_hysteresis = Parameter.boolean(
        name="use_hysteresis",
        description=(
            "Use hysteresis (double) thresholding instead of a single threshold. "
            "Seeds pixels >= high_thr and grows into neighbours >= low_thr."
        ),
        default=False,
    )
    high_thr = Parameter.number(
        name="high_thr",
        description="Hysteresis high (seed) threshold in [0, 1]. Default 0.4.",
        default=HYSTERESIS_DEFAULTS["high_thr"],
    )
    low_thr = Parameter.number(
        name="low_thr",
        description="Hysteresis low (extension) threshold in [0, 1]. Default 0.2.",
        default=HYSTERESIS_DEFAULTS["low_thr"],
    )
    max_dist_m = Parameter.number(
        name="max_dist_m",
        description="Maximum grow distance in metres for hysteresis. Default 250 m.",
        default=HYSTERESIS_DEFAULTS["max_dist_m"],
    )
    connectivity = Parameter.number(
        name="connectivity",
        description="Hysteresis binary structure connectivity: 1=cross (4-connected), 2=full 3x3 (8-connected). Default 2.",
        default=HYSTERESIS_DEFAULTS["connectivity"],
    )
    output, s2_cube = _water_extent_for_month(
        connection,
        spatial_extent,
        region,
        start_date,
        date_shift(start_date, value=1, unit="month"),
        cloud_cover,
        water_threshold,
        only_s1,
        filter_type=filter_type,
        use_hysteresis=use_hysteresis,
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    )
    # Merge S2 RGB bands into the UDP output when rgb_processing is True.
    # s2_cube is always available in the UDP path (only_s1 is a Parameter here, not a Python bool).
    s2_rgb = s2_cube.filter_bands(["B04", "B03", "B02"])
    output = if_(rgb_processing, output.merge_cubes(s2_rgb), output)
    returns = {
        "description": "A data cube with the newly computed values.\n\nAll dimensions stay the same, except for the dimensions specified in corresponding parameters. There are three cases how the dimensions can change:\n\n1. The source dimension is the target dimension:\n   - The (number of) dimensions remain unchanged as the source dimension is the target dimension.\n   - The source dimension properties name and type remain unchanged.\n   - The dimension labels, the reference system and the resolution are preserved only if the number of values in the source dimension is equal to the number of values computed by the process. Otherwise, all other dimension properties change as defined in the list below.\n2. The source dimension is not the target dimension. The target dimension exists with a single label only:\n   - The number of dimensions decreases by one as the source dimension is 'dropped' and the target dimension is filled with the processed data that originates from the source dimension.\n   - The target dimension properties name and type remain unchanged. All other dimension properties change as defined in the list below.\n3. The source dimension is not the target dimension and the latter does not exist:\n   - The number of dimensions remain unchanged, but the source dimension is replaced with the target dimension.\n   - The target dimension has the specified name and the type other. All other dimension properties are set as defined in the list below.\n\nUnless otherwise stated above, for the given (target) dimension the following applies:\n\n- the number of dimension labels is equal to the number of values computed by the process,\n- the dimension labels are incrementing integers starting from zero,\n- the resolution changes, and\n- the reference system is undefined.",
        "schema": {"type": "object", "subtype": "datacube"},
    }
    udp = build_process_dict(
        output,
        "worldwater_water_extent",
        "Computes water extent for a given month, provided by DHI.",
        description="Computes water extent for a given month, provided by DHI.",
        parameters=[
            start_date,
            spatial_extent,
            region,
            only_s1,
            rgb_processing,
            cloud_cover,
            water_threshold,
            filter_type,
            use_hysteresis,
            high_thr,
            low_thr,
            max_dist_m,
            connectivity,
        ],
        returns=returns,
    )
    return udp


def generate_water_indicators_udp(connection: openeo.Connection):
    """ water indicators UDP that processes multiple months and returns water indicators."""
    start_date = Parameter.date(name="start_date", description="Start date of the processing period.")
    end_date = Parameter.date(name="end_date", description="End date of the processing period (exclusive).")
    spatial_extent = Parameter.spatial_extent()

    region = Parameter.string(
        name="region",
        description="Eco-Region on which to compute water probability",
        default="Deserts",
        values=list(LOOKUPTABLE.keys()),
    )
    only_s1 = Parameter.boolean(
        name="only_s1",
        description="Use only Sentinel-1 (True) or both Sentinel-1 and Sentinel-2 (False).",
        default=False,
    )
    cloud_cover = Parameter.number(
        name="cloud_cover",
        description="Maximum cloud cover percentage for Sentinel-2.",
        default=85,
    )
    filter_type = Parameter.string(
        name="filter_type",
        description=(
            "Spatial smoothing filter applied to the water probability before thresholding. "
            f"One of {FILTER_TYPES}. Use 'none' to skip smoothing."
        ),
        default="none",
        values=FILTER_TYPES,
    )
    use_hysteresis = Parameter.boolean(
        name="use_hysteresis",
        description="Use hysteresis (double) thresholding instead of a single threshold.",
        default=False,
    )
    water_threshold = Parameter.number(
        name="water_threshold",
        description="Water probability threshold (0-100). Used when use_hysteresis=False.",
        default=75,
    )
    high_thr = Parameter.number(
        name="high_thr",
        description="Hysteresis high (seed) threshold in [0, 1]. Default 0.6.",
        default=HYSTERESIS_DEFAULTS["high_thr"],
    )
    low_thr = Parameter.number(
        name="low_thr",
        description="Hysteresis low (extension) threshold in [0, 1]. Default 0.4.",
        default=HYSTERESIS_DEFAULTS["low_thr"],
    )
    max_dist_m = Parameter.number(
        name="max_dist_m",
        description="Maximum grow distance in metres for hysteresis. Default 250 m.",
        default=HYSTERESIS_DEFAULTS["max_dist_m"],
    )
    connectivity = Parameter.number(
        name="connectivity",
        description="Hysteresis connectivity: 1=cross (4-connected), 2=full 3x3 (8-connected). Default 2.",
        default=HYSTERESIS_DEFAULTS["connectivity"],
    )
    output = _water_extent_multiple_months(
        connection,
        start_date,
        end_date,
        spatial_extent,
        region,
        cloud_cover,
        threshold=water_threshold,
        only_s1=only_s1,
        filter_type=filter_type,
        gaussian_sigma=0.85,
        uniform_size=3,
        median_size=3,
        use_hysteresis=use_hysteresis,
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    ) 
    binary_monthly = output.filter_bands(["surface_water"])
    indicators = _compute_water_indicators(binary_monthly)
    output = output.merge_cubes(indicators) 

    # PS: merging the data need to be fixed so because the water indicators do not share the same time dim as the monthly layers 
    # which results in duplicated water indicators across the time dim.
    returns = {
        "description": (
            "A data cube with monthly water probailities and water maks (band: water_prob_sum, surface_water"
            "time-stacked) merged with annual water indicator bands: mean_swf, "
            "annual_classification, min_water_extent, max_water_extent, water_seasonality."
        ),
        "schema": {"type": "object", "subtype": "datacube"},
    }
    udp = build_process_dict(
        output,
        "worldwater_water_indicators",
        "Computes multi-month water indicators, provided by DHI.",
        description=(
            "Computes monthly water probabilities, binary water masks and annual water indicators "
            "(Annual Water Classification, Seasonal Water Classification, Minimum Water Extent, Maximum Water Extent, Average Surface Water Frequency) for a given period."
            "provided by DHI."
        ),
        parameters=[
            start_date,
            end_date,
            spatial_extent,
            region,
            only_s1,
            cloud_cover,
            filter_type,
            use_hysteresis,
            water_threshold,
            high_thr,
            low_thr,
            max_dist_m,
            connectivity,
        ],
        returns=returns,
    )
    return udp


def main(
    backend: str,
    start,
    end,
    region: str,
    geometry: str,
    rgb_processing: bool,
    cloud_cover: int,
    threshold: Union[int, float],
    only_s1: bool = False,
    filter_type: str = "none",
    gaussian_sigma: float = 0.85,
    uniform_size: int = 3,
    median_size: int = 3,
    use_hysteresis: bool = False,
    high_thr: float = HYSTERESIS_DEFAULTS["high_thr"],
    low_thr: float  = HYSTERESIS_DEFAULTS["low_thr"],
    max_dist_m: float = HYSTERESIS_DEFAULTS["max_dist_m"],
    connectivity: int = HYSTERESIS_DEFAULTS["connectivity"],
) -> str:

    connection = openeo.connect(backend).authenticate_oidc()
    if os.path.isfile(geometry):
        print("Processing")
    else:
        print("The AOI file does not exist. Please insert the correct GeoJSON file.")
        exit()

    print("Start and end dates:", start, end)
    print("Region:", region)
    print("Threshold:", threshold)
    print("Cloud Cover:", cloud_cover)

    months_list = [start]
    current_month = start

    # Get a list of months to be processed
    while current_month < end:
        current_month += relativedelta(months=1)
        months_list.append(current_month)

    # Iterate toolbox for each month
    for month_start in months_list:
        month_end = month_start + relativedelta(months=1)
        output_folder = _water_extent(
            connection,
            month_start,
            month_end,
            start,
            end,
            geometry,
            region,
            threshold,
            cloud_cover,
            rgb_processing,
            only_s1=only_s1,
            filter_type=filter_type,
            gaussian_sigma=gaussian_sigma,
            uniform_size=uniform_size,
            median_size=median_size,
            use_hysteresis=use_hysteresis,
            high_thr=high_thr,
            low_thr=low_thr,
            max_dist_m=max_dist_m,
            connectivity=connectivity,
        )

    print("Successfully finished! The output files are located at:", output_folder)
    return output_folder

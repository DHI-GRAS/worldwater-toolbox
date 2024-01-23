# -*- coding: utf-8 -*-
"""
Created on Jun 2023
@author: Andrea Sulova
"""

import os
import glob
import numpy as np
import rioxarray
import xarray as xr
import rasterio
import geopandas as gpd
from dateutil.relativedelta import *
import openeo
from openeo import DataCube
from openeo.api.process import Parameter
from openeo.extra.spectral_indices.spectral_indices import append_indices, compute_indices
from openeo.processes import if_, exp, array_element, log, count, gte, eq, sum, date_shift


LOOKUPTABLE = {

        "Deserts": {
            "S1": lambda vh, vv: (1 / (1 + exp(- (-7.03 + (-0.44 * vv))))),
            "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.133 + (-5.92 * ndvi) + (14.82 * ndwi))))),
            "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-3.69 + (-0.25 * vv) + (0.47 * ndvi) + (15.3 * ndwi))))),
        },
        "Mountain": {
            "S1": lambda vh, vv: (1 / (1 + exp(- (-3.76 + (-0.262 * vv))))),
            "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.262 + (0.75 * ndvi) + (12.65 * ndwi))))),
            "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-1.13 + (-0.11 * vv) + (3.03 * ndvi) + (13.21 * ndwi)))))
        },
        "Tropical forest":
            {
                "S1": lambda vh, vv: (1 / (1 + exp(- (-5.8 + (-0.415 * vv))))),
                "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.344 + (2.886 * ndvi) + (11.91 * ndwi)))))*100,
                "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-3.25 + (-0.23 * vv) + (4.17 * ndvi) + (9.5 * ndwi)))))
        },
         "Tropical savanna":
            {
                "S1": lambda vh, vv: (1 / (1 + exp(- (-7.0 + (-0.444 * vv))))),
                "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.344 + (2.886 * ndvi) + (11.91 * ndwi)))))*100,
                "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-1.06 + (-0.17 * vv) + (3.82* ndvi) + (14.4* ndwi)))))
        },
        "Subtropical savanna":
            {
                "S1": lambda vh, vv: (1 / (1 + exp(- (-7.17 + (-0.48 * vv))))),
                "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.845 + (2.14 * ndvi) + (13.5 * ndwi))))),
                "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-2.64 + (-0.23 * vv) + (8.6 * ndwi))))),
        },
        "Subtropical forest":
            {
                "S1": lambda vh, vv: (1 / (1 + exp(- (-6.67 + (-6.67* vv))))),
                "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.712 + (-1.133 * ndvi) + (7.16 * ndwi))))),
                "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-2.72 + (-0.22 * vv) + (-0.49  * ndvi) + 4.55 * ndwi))))
        },
        "Temperate broadleaf":
            {
            "S1": lambda  vh, vv: 1 / (1 + exp(- (-8.82 + (-0.58 * vv)))),
            "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (-0.013 + (5.38 * ndvi) + (13.79 * ndwi)))),
            "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-2.7 + (-0.2 * vv) + (3.6 * ndvi) + (9.73 * ndwi))))
        },
        "Temperate grassland":
            {
            "S1": lambda  vh, vv: 1 / (1 + exp(- (-7.01 + (-0.426 * vv)))),
            "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (1.286 + (8.74 * ndvi) + (23.217 * ndwi)))),
            "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-3.43 + (-0.25 * vv) + (11.74 * ndvi) + (22.035 * ndwi))))
            },
        "Tundra":
            {
            "S1": lambda vh, vv: 1 / (1 + exp(- (-4.57 + (-0.13 * vh) + (-0.11 * vv)))),
            "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (1.07 + (3.6 * ndwi)))),
#             "S1_S2": lambda hh, hv, ndwi: 1 / (1 + exp(- (-3.8 + (-0.15 * hh) + (-0.03 * hv) + (2.5 * ndwi))))
            "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-3.8 + (-0.15 * vv) + (-0.03 * ndvi) + (2.5 * ndwi))))
            }
        }
        
def _water_indicators(path, prefix):
    
    """
    Calculate Calculate water frequency (0-1), minimum water extent, maximum water extent and water classification
        (1) Permanent water
        (2) Seasonal water
        (0) No water

    Parameters
    ----------
    path : STRING
        Path where water masks are stored

    Returns
    -------
    water_frequency : NUMPY ARRAY
        Water frequency as int
        
    max_extent : NUMPY ARRAY
        Minimum water extent as int
        
    min_extent : NUMPY ARRAY
        Maximum water extent as int
        
    classification : NUMPY 
        Water classification as int

    """

    raster_list = []
    watermasks = glob.glob(path + '/water*.tif')
    print(watermasks)
    for watermask in watermasks:
        raster = rioxarray.open_rasterio(watermask)
        raster_list.append(raster)

    # Concatenate the rioxarray objects along the 'band' dimension
    stacked_raster = xr.concat(raster_list, dim='band')

    # Sum the pixel values along the 'band' dimension
    sum_raster = stacked_raster.sum(dim='band')
    n_months = len(raster_list)

    # Calculate indices of occurrence, max and min extent
    water_occurrence = sum_raster / n_months
    min_extent = xr.where(sum_raster == n_months, 1, 0)
    max_extent = xr.where(water_occurrence > 0, 1, 0)

    # Classification: Permanent water (1) and seasonal water (2)
    classification = xr.zeros_like(water_occurrence).where(max_extent == 0, 2, 0).where(min_extent == 0, 1, 0)

    # Save files to folder
    water_occurrence.astype(np.float32).to_dataset(name='water_occurrence').rio.to_raster(
        path + f'/index_water_occurrence_{prefix}.tif')
    min_extent.astype(np.float32).to_dataset(name='min_extent').rio.to_raster(path + f'/water_min_extent_{prefix}.tif')
    max_extent.astype(np.float32).to_dataset(name='max_extent').rio.to_raster(path + f'/water_max_extent_{prefix}.tif')
    classification.astype(np.float32).to_dataset(name='classification').rio.to_raster(
        path + f'/water_classification_{prefix}.tif')

    return water_occurrence, min_extent, max_extent, classification
                                         

def hillshade_mask(connection, spatial_extent, start_date_exclusion, month_end, cloud_cover, use_sentinelhub=True):

    # Hill-shade masking
    collection = 'SENTINEL2_L2A'
    if use_sentinelhub:
        collection = 'SENTINEL2_L2A_SENTINELHUB'
        
    # angles are lower resolution, so we load them separately
    s2_angles = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=['B01','sunAzimuthAngles', 'sunZenithAngles'],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover}
    )

    dem_cube = connection.load_collection("COPERNICUS_30",
                                          spatial_extent=spatial_extent,
                                          temporal_extent=["2010-01-01", "2030-12-31"])
    dem_cube = dem_cube.max_time()
    # Resample s2 cube (Azimuth and Zenith) to 30m
    s2_cube_30 = s2_angles.resample_spatial(resolution=30, method='average')
    # DEM WGS84 to DEM UTM from S2 cube
    # dem_cube_s2 = dem_cube.resample_cube_spatial(s2_cube_30)
    # Merge 30m DEM and 30m s2 cube due Azimuth and Zenith
    merged_cube = s2_cube_30.merge_cubes(dem_cube)
    # Apply the hill-shade udf
#     print('udf', os.path.dirname(os.path.abspath(__file__)) + '/udf.py')
    process = openeo.UDF.from_file(os.path.dirname(os.path.abspath(__file__)) + '/udf.py', runtime="Python")
    hillshade = merged_cube.apply_neighborhood(
        process=process,
        size=[
            {"dimension": "t", "value": "P1D"},
            {"dimension": "x", "unit": "px", "value": 256},
            {"dimension": "y", "unit": "px", "value": 256}
        ],
        overlap=[
            {"dimension": "x", "unit": "px", "value": "8"},
            {"dimension": "y", "unit": "px", "value": "8"}
        ]
    )
    # Rename bands in a hill shade cube
    hillshade = hillshade.rename_labels("bands",
                                        ["hillshade_mask", "sunAzimuthAngles", "sunZenithAngles", "DEM"])
    return hillshade

def masked_s2_cube(connection: openeo.Connection, spatial_extent, start_date_exclusion, month_end, cloud_cover, use_sentinelhub=True):

    collection = 'SENTINEL2_L2A'
    if use_sentinelhub:
        collection = 'SENTINEL2_L2A_SENTINELHUB'

    col_metadata = connection.describe_collection(collection)
    bands = col_metadata.get('cube:dimensions',{}).get('bands',{}).get('values',[])


    # Loading S2 collection
    s2_cube = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=['B02', 'B03', 'B04', 'B08'],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover}
    )

    # Loading S2 collection for cloud masking processing
    if 'CLP' in bands:
        cloudmask_bands = ['CLP', 'SCL']
    else:
        cloudmask_bands = ['SCL']

    s2_cube_masking = connection.load_collection(
        collection,
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=cloudmask_bands,
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover}
    )
    # Applying mask from sen2cor
    scl = s2_cube_masking.band("SCL")
    mask_scl = (scl == 3) | (scl == 8) | (scl == 9) | (scl == 10) | (scl == 11)

    if 'CLP' in bands:
        clp = s2_cube_masking.band("CLP")
        mask_clp = mask_scl | (clp / 255) > 0.3
    else:
        mask_clp = mask_scl

    if 'sunAzimuthAngles' in bands and 'sunZenithAngles' in bands:
        hillshade = hillshade_mask(connection,spatial_extent,start_date_exclusion,month_end,cloud_cover).filter_bands(
            "hillshade_mask").linear_scale_range(0,10,0,10)
        the_mask = mask_clp.merge_cubes(hillshade.drop_dimension("bands"), overlap_resolver="sum")
    else:
        the_mask = mask_clp

    #it is better to create one combined mask, which can be used to reduce data loading in full resolution cube

    # Mask s2 cube with a hillshade mask
    s2_cube_masked = s2_cube.mask(the_mask)


    # Replace 0 to nan in s2 cubes
    s2_cube = s2_cube_masked.apply(lambda x: if_(x.neq(0),x))
    return s2_cube


def generate_water_extent_udp(connection: openeo.Connection):
    from openeo.rest.udp import build_process_dict


    DATE_SCHEMA = {
                    "type": "string",
                    "format": "date",
                    "subtype": "date"
                }
    bbox_schema = {
        "title": "Bounding Box",
        "type": "object",
        "subtype": "bounding-box",
        "required": [
            "west",
            "south",
            "east",
            "north"
        ],
        "properties": {
            "west": {
                "description": "West (lower left corner, coordinate axis 1).",
                "type": "number"
            },
            "south": {
                "description": "South (lower left corner, coordinate axis 2).",
                "type": "number"
            },
            "east": {
                "description": "East (upper right corner, coordinate axis 1).",
                "type": "number"
            },
            "north": {
                "description": "North (upper right corner, coordinate axis 2).",
                "type": "number"
            },
            "crs": {
                "description": "Coordinate reference system of the extent, specified as as [EPSG code](http://www.epsg-registry.org/) or [WKT2 CRS string](http://docs.opengeospatial.org/is/18-010r7/18-010r7.html). Defaults to `4326` (EPSG code 4326) unless the client explicitly requests a different coordinate reference system.",
                "anyOf": [
                    {
                        "title": "EPSG Code",
                        "type": "integer",
                        "subtype": "epsg-code",
                        "minimum": 1000,
                        "examples": [
                            3857
                        ]
                    },
                    {
                        "title": "WKT2",
                        "type": "string",
                        "subtype": "wkt2-definition"
                    }
                ],
                "default": 4326
            }
        }
    }
    start_date = Parameter(
        name="start_date", description="The start date.",
        schema=DATE_SCHEMA
    )
    spatial_extent = Parameter(name="bbox", schema=bbox_schema, description="The spatial extent, as a bounding box")
    region = Parameter.string("region",description="Eco-Region on which to compute water probability", default="Deserts",values=LOOKUPTABLE.keys())
    only_s1 = Parameter.boolean("only_s1",description="Boolean variable to specifiy if only Sentinel-1 data will be used (True) or both Sentinel-1 and Sentinel-2", default=False)
    rgb_processing = Parameter.boolean("rgb_processing",description="Boolean variable to specifiy if Sentinel-2 rgb image will be generated")
    output, s2_cube = _water_extent_for_month(connection, spatial_extent, region, start_date, date_shift(start_date,value=1,unit="month"), 85, 75, True, rgb_processing, only_s1)

    udp = build_process_dict(output,"worldwater_water_extent","Computes water extent for a given month.")
    return udp


def _water_extent_multiple_months(connection: openeo.Connection, month_start, month_end, geometry, region, cloud_cover, only_s1, use_sentinelhub=True):
    spatial_extent = _get_spatial_extent(geometry)
    start_date_exclusion = date_shift(month_start,value=-1,unit="month")

    s2_cube = masked_s2_cube(connection, spatial_extent, start_date_exclusion, month_end, cloud_cover,
                             use_sentinelhub=use_sentinelhub)

    s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube, region)

    # Monthly median s2 image
    s2_median_water = s2_cube_water.aggregate_temporal_period("month","median")
    ndxi_median = ndxi_cube.aggregate_temporal_period("month","median")

    s1_cube = sentinel1_preprocessing(connection, month_end, month_start, spatial_extent, use_sentinelhub, region)

    # Normalized radar back-scatter
    # Calculate S1 median mosaic
    s1_median = s1_cube.aggregate_temporal_period("month","median")

    merge_all = _water_probability(s1_median, s2_median_water, ndxi_median,region,only_s1)
    return merge_all


def _water_extent(connection: openeo.Connection, month_start, month_end, start, end, geometry, region, threshold, 
                  cloud_cover, rgb_processing, use_sentinelhub=True, only_s1=False, output_name_o=''):
        
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

    # Read geometry from GeoJSON file
    spatial_extent = _get_spatial_extent(geometry)
    output, s2_cube = _water_extent_for_month(connection, spatial_extent, region, month_start, month_end, cloud_cover, threshold,
                                              use_sentinelhub, rgb_processing, only_s1)
    
    # Output folder
    region_naming = '_'.join(region.split(" "))
    output_folder = output_name_o + '/' + f'WWT_{region_naming}_{start.strftime("%Y_%m")}_{end.strftime("%Y_%m")}'

    prefix = 'water_onlyS1' if only_s1 else 'water_'
    filename = prefix + '_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_start.strftime("%Y_%m"))
    # Send jobs to the backend
    print('Water mask processing between: ' + str(month_start.strftime("%Y_%m")) + ' and ' + str(month_end.strftime("%Y_%m")))
    my_job = output.create_job(
        title = prefix + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m")),
        out_format="GTiff",
        job_options={
            "executor-memory": "1g",
            "executor-memoryOverhead": "3g",
            "executor-cores": 1
        },
        **{'filename_prefix': filename}
    )
    results = my_job.start_and_wait().get_results()
    results.download_files(output_folder)

    if rgb_processing:
        # Process S2 median image if rgb_processing is True
        print('S2 median image processing between: ' + str(month_start.strftime("%Y_%m")) + ' and ' + str(month_end.strftime("%Y_%m")))
        s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time()
        my_job = s2_cube_median.create_job(
            title='median_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m")),
            out_format="GTiff",
            job_options={"node_caching": True},
            **{'filename_prefix': 'median_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m"))}
        )
        results = my_job.start_and_wait().get_results()
        results.download_files(output_folder)

    return output_folder


def _water_extent_for_month(connection, spatial_extent, region, month_start, month_end, cloud_cover, threshold, use_sentinelhub, rgb_processing, only_s1=False):
    month_start = str(month_start) if not isinstance(month_start,Parameter) else month_start
    start_date_exclusion = date_shift(month_start, value=-1, unit="month")

    # Normalized radar back-scatter
    s1_cube = sentinel1_preprocessing(connection, month_end, month_start, spatial_extent, use_sentinelhub, region)
    # Calculate S1 median mosaic
    s1_median = s1_cube.median_time()

    if only_s1:
        merge_all = _water_probability(s1_median, None, None, region, True)    
        if rgb_processing:
            s2_cube = masked_s2_cube(connection, spatial_extent, start_date_exclusion, month_end,cloud_cover, use_sentinelhub = use_sentinelhub)
            s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube,region)
            s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time()
        else:
            s2_cube_median = None
    else:
        s2_cube = masked_s2_cube(connection, spatial_extent, start_date_exclusion, month_end, cloud_cover,
                                use_sentinelhub=use_sentinelhub)
        s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(s2_cube, region)
        
        s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time() if rgb_processing else None
        s2_median_water = s2_cube_water.filter_temporal([month_start, month_end]).median_time()
        ndxi_median = ndxi_cube.filter_temporal([month_start, month_end]).median_time()
        merge_all = _water_probability(s1_median, s2_median_water, ndxi_median, region, False)
    
    if 'ESA_WORLDCOVER_10M_2020_V1' in connection.list_collection_ids():
        # Mask built-up area using ESA world cover layer
        worldcover_cube = connection.load_collection(
            "ESA_WORLDCOVER_10M_2020_V1",
            temporal_extent=['2020-12-30', '2021-01-01'],
            spatial_extent=spatial_extent,
            bands=["MAP"]
        )

        builtup_mask = worldcover_cube.band("MAP") == 50
        water_probability = merge_all.mask(builtup_mask.max_time().resample_cube_spatial(merge_all))
    else:
        water_probability = merge_all
    water_probability = water_probability.rename_labels("bands", ["water_prob_sum"])
    
    # Apply custom threshold to water probability layer
    water_mask = water_probability.apply(lambda x: x['water_prob_sum']*100 > threshold)
    water_mask = water_mask.rename_labels("bands", ["surface_water"])
    output = water_probability.merge_cubes(water_mask)
    output = output * 1.0
    
    return output, s2_cube_median


def _water_probability(s1_median, s2_median_water, ndxi_median, region, only_s1):
    def log_(x):
        return 10 * log(x, 10)

    s1_median = s1_median.apply(log_)

    # Create a water extent for the S1 collection using logistic expression
    s1_median_water = 0
    for key in LOOKUPTABLE.keys():
        #the if/else lookup is done in the process graph, for the UDP, and can not yet be done in the callback
        if region=='Tundra':
            s1_median_water = if_(eq(region, key), s1_median.reduce_dimension(reducer=lambda data: LOOKUPTABLE[key]["S1"](vh=data[0], vv=data[1]), dimension="bands"), s1_median_water)
        else:
            s1_median_water = if_(eq(region, key), s1_median.reduce_dimension(reducer=lambda data: LOOKUPTABLE[key]["S1"](vh=data[0], vv=data[0]), dimension="bands"), s1_median_water)
    
    if only_s1:
#         merge_all = s1_median_water.add_dimension("bands", "var", type="bands")
        s1_water = DataCube(s1_median_water.pgnode,s1_median.connection, metadata=s1_median.metadata.reduce_dimension("bands"))
        return s1_water.add_dimension("bands", "var", type="bands")
    
    # exclusion_mask = (s1_median_water > 0.5) & (s2_cube_swf < 0.33)
    # TODO: below line was not used
    # s1_median_water_mask = s1_median_water.mask(exclusion_mask.resample_cube_spatial(s1_median_water))
    # Calculate water extent for the S1 & S2 collection using logistic expression from the lookup table
    if region=='Tundra':
        s1_s2_cube = (
            s1_median.filter_bands(["HH", "HV"])
            .merge_cubes(ndxi_median)
        )
        s1_s2_water = 0
        for key in LOOKUPTABLE.keys():
            # the if/else lookup is done in the process graph, for the UDP, and can not yet be done in the callback
            s1_s2_water = if_(eq(region, key), s1_s2_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S1_S2"](vv=data[0], ndvi=data[1], ndwi=data[2]), dimension="bands"), s1_s2_water)

    else:
        s1_s2_cube = (
            s1_median.filter_bands(["VV"])
            .merge_cubes(ndxi_median)
        )
        s1_s2_water = 0
        for key in LOOKUPTABLE.keys():
            # the if/else lookup is done in the process graph, for the UDP, and can not yet be done in the callback
            s1_s2_water = if_(eq(region, key), s1_s2_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S1_S2"](vv=data[0], ndvi=data[2], ndwi=data[1]), dimension="bands"), s1_s2_water)
    
    s1_s2_water = DataCube(s1_s2_water.pgnode,s1_s2_cube.connection, metadata=s1_s2_cube.metadata.reduce_dimension("bands"))
    merged = s1_s2_water.add_dimension("bands", "s1_s2_water", type="bands").merge_cubes(s2_median_water).merge_cubes(s1_median_water)

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

    merge_all = merged.apply_dimension(combined, dimension='bands')
    # Composites water probabilities according to data availability and accuracy.
    """
        s1_s2_mask = s1_s2_water >= 0
        s2_mask = s2_median_water.mask(s1_s2_mask) >= 0
        s1_mask = s1_median_water.mask(s1_s2_mask).mask(s2_mask) >= 0
        s1_s2_masked = s1_s2_water.mask(s1_s2_mask.apply(lambda x: x.eq(0)), replacement=0)
        s2_masked = s2_median_water.mask(s2_mask.apply(lambda x: x.eq(0)), replacement=0)
        s1_masked = s1_median_water.mask(s1_mask.apply(lambda x: x.eq(0)), replacement=0)
        merge_all = (
            s1_s2_masked.merge_cubes(s2_masked, overlap_resolver="sum")
            .merge_cubes(s1_masked, overlap_resolver="sum")
        )
        """
    
    
    return merge_all


def sentinel1_preprocessing(connection, month_end, month_start, spatial_extent, use_sentinelhub, region):
    # Loading S1 collection
    pol = 'DH' if region == 'Tundra' else 'DV' 
    s1_bands = ['HH', 'HV'] if region == 'Tundra' else ['VV'] 
    s1_cube = connection.load_collection(
        'SENTINEL1_GRD',
        spatial_extent=spatial_extent,
        temporal_extent=[month_start, month_end],
        bands=s1_bands,
        properties={"polarization": lambda p: p == pol}
    )
    
    coefficient = 'gamma0-terrain' if use_sentinelhub else 'sigma0-ellipsoid'
    # Apply terrain correction for back-scatter values
    s1_cube = s1_cube.sar_backscatter(coefficient=coefficient, mask=use_sentinelhub, elevation_model="COPERNICUS_30",
                                      options=
                                      {"implementation_version": "2", "tile_size": 256, "otb_memory": 1024,
                                       "debug": True})
    s1_cube = s1_cube.rename_labels("bands", s1_bands + ["mask"])
    if use_sentinelhub:
        def apply_mask(bands):
            return if_(bands.array_element(1) != 2, bands)

        s1_cube = s1_cube.apply_dimension(apply_mask, dimension="bands")
    else:
        pass
    return s1_cube


def _get_spatial_extent(geometry):
    gdf = gpd.read_file(geometry)
    gdf = gdf.to_crs(4326)
    bbox = gdf.geometry.total_bounds
    spatial_extent = {'west': bbox[0], 'east': bbox[2], 'south': bbox[1], 'north': bbox[3], 'crs': 4326}
    return spatial_extent


def s2_water_processing(s2_cube,region):
    # Adding NDVI and NDWI to datacube
    s2_indices = ["NDWI"] if region == 'Tundra' else ["NDWI", "NDVI"]
    s2_cube = append_indices(s2_cube, s2_indices)
    ndxi_cube = compute_indices(s2_cube, s2_indices)
    s2_cube = s2_cube.rename_labels("bands", ["B02", "B03", "B04", "B08"]+s2_indices)

    # Create water extent for the S2 collection using logistic expressions from Lookup Table

    s2_cube_water = 0
    for key in LOOKUPTABLE.keys():
        # the if/else lookup is done in the process graph, for the UDP, and can not yet be done in the callback
        if region == 'Tundra':
            s2_cube_water = if_(eq(region, key), ndxi_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S2"](ndwi=data[0], ndvi=data[0]), dimension="bands"),s2_cube_water)
        else:
            s2_cube_water = if_(eq(region, key), ndxi_cube.reduce_dimension(
                reducer=lambda data: LOOKUPTABLE[key]["S2"](ndwi=data[0], ndvi=data[1]), dimension="bands"),s2_cube_water)

    s2_cube_water = DataCube(s2_cube_water.pgnode, ndxi_cube.connection,
                           metadata=ndxi_cube.metadata.reduce_dimension("bands"))

    s2_cube_water = s2_cube_water.add_dimension("bands", "water_prob", type="bands")
    # Generate a binary water mask by applying a threshold to the water probability band
    s2_cube_water_threshold = s2_cube_water.apply_dimension(dimension="bands", process=lambda x: if_(x > 0.75, x, 0))
    s2_cube_water_threshold = s2_cube_water_threshold.rename_labels("bands", ["w_T75"])
    # Calculate the water frequency
    s2_cube_water_sum = s2_cube_water_threshold.reduce_dimension(reducer="sum", dimension="t")
    s2_cube_water_sum = s2_cube_water_sum.rename_labels("bands", ["sum"])
    s2_count = s2_cube.band("B08")
    s2_count = s2_count.reduce_dimension(reducer=lambda data: data.count(), dimension="t")
    # TODO: sum divided by count equals mean??
    # s2_cube_swf = s2_cube_water_threshold.reduce_dimension(reducer="mean", dimension="t")
    # TODO s2_cube_swf not used in the end?
    # s2_cube_swf = s2_cube_water_sum / s2_count
    # s2_cube_swf = s2_cube_swf.rename_labels("bands", ["swf"])
    return s2_cube, ndxi_cube , s2_cube_water


def main(backend, start, end, region, geometry, rgb_processing, cloud_cover, threshold, use_sentinelhub=True, only_s1=False):
    # Connect to openEO backends

    connection = openeo.connect(backend).authenticate_oidc()

    # Check if the AOI file exists
    if os.path.isfile(geometry):
        print('Processing')
    else:
        print('The AOI file does not exist. Please insert the correct GeoJSON file.')
        exit()

    # Summary of the query
    print('Start and end dates:', start, end)
    print('Region:', region)
    print('Threshold:', threshold)
    print('Cloud Cover:', cloud_cover)

    months_list = [start]
    current_month = start

    # Get a list of months to be processed
    while current_month < end:
        current_month += relativedelta(months=1)
        months_list.append(current_month)

    # Iterate toolbox for each month
    for month_start in months_list:
        month_end = month_start + relativedelta(months=1)
        output_folder = _water_extent(connection, month_start, month_end, start, end, geometry, region, threshold, cloud_cover, rgb_processing, use_sentinelhub=use_sentinelhub, only_s1=only_s1)

    # Calculate water frequency (0-1), minimum water extent, maximum water extent, and water classification
    _water_indicators(output_folder, only_s1)

    print('Successfully finished! The output files are located at:', output_folder)

 

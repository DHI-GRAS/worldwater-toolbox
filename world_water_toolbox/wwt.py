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
from openeo.extra.spectral_indices.spectral_indices import append_indices
from openeo.processes import if_, exp, array_element, log, count, gte, eq, sum


LOOKUPTABLE = {
            "Deserts": {
                "S1": lambda vh, vv: 1 / (1 + exp(- (-7.03 + (-0.44 * vv)))),
                "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (0.133 + (-5.92 * ndvi) + (14.82 * ndwi)))),
                "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-3.69 + (-0.25 * vv) + (0.47 * ndvi) + (15.3 * ndwi)))),
            },

            "Mountain": {
                "S1": lambda vh, vv: 1 / (1 + exp(- (-3.76 + (-0.262 * vv)))),
                "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (0.262 + (0.75 * ndvi) + (12.65 * ndwi)))),
                "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-1.13 + (-0.11 * vv) + (3.03 * ndvi) + (13.21 * ndwi)))),
            },

            "Tropical_forest":
                {
                    "S1": lambda vh, vv: (1 / (1 + exp(- (-5.8 + (-0.415 * vv))))),
                    "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.344 + (2.886 * ndvi) + (11.91 * ndwi)))))*100,
                    "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-3.25 + (-0.23 * vv) + (4.17 * ndvi) + (9.5 * ndwi))))),
            },

             "Tropical_savanna":
                {
                    "S1": lambda vh, vv: (1 / (1 + exp(- (-7.0 + (-0.444 * vv))))),
                    "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.344 + (2.886 * ndvi) + (11.91 * ndwi)))))*100,
                    "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-1.06 + (-0.17 * vv) + (3.82* ndvi) + (14.4* ndwi))))),
            },

            "Subtropical_savanna":
                {
                    "S1": lambda vh, vv: (1 / (1 + exp(- (-7.17 + (-0.48 * vv))))),
                    "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.845 + (2.14 * ndvi) + (13.5 * ndwi))))),
                    "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-2.64 + (-0.23 * vv) + (8.6 * ndwi))))),
            },

            "Subtropical_forest":
                {
                    "S1": lambda vh, vv: (1 / (1 + exp(- (-6.67 + (-6.67* vv))))),
                    "S2": lambda ndvi, ndwi: (1 / (1 + exp(- (0.712 + (-1.133 * ndvi) + (7.16 * ndwi))))),
                    "S1_S2": lambda vv, ndvi, ndwi: (1 / (1 + exp(- (-2.72 + (-0.22 * vv) + (-0.49  * ndvi) + 4.55 * ndwi)))),
            },

            "Temperate_broadleaf":
                {
                "S1": lambda  vh, vv: 1 / (1 + exp(- (-8.82 + (-0.58 * vv)))),
                "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (-0.013 + (5.38 * ndvi)) + (13.79 * ndwi))),
                "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-2.7 + (-0.2 * vv)) + (3.6 * ndvi)) + (9.73 * ndwi))
            },

            "Temperate_grassland":
                {
                "S1": lambda  vh, vv: 1 / (1 + exp(- (-7.01 + (-0.426 * vv)))),
                "S2": lambda ndvi, ndwi: 1 / (1 + exp(- (1.286 + (8.74 * ndvi)) + (23.217 * ndwi))),
                "S1_S2": lambda vv, ndvi, ndwi: 1 / (1 + exp(- (-3.43 + (-0.25 * vv)) + (11.74 * ndvi)) + (22.035 * ndwi))
                }
            }


def _water_indicators(path):
    
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
        path + '/index_water_occurrence.tif')
    min_extent.astype(np.float32).to_dataset(name='min_extent').rio.to_raster(path + '/index_min_extent.tif')
    max_extent.astype(np.float32).to_dataset(name='max_extent').rio.to_raster(path + '/index_max_extent.tif')
    classification.astype(np.float32).to_dataset(name='classification').rio.to_raster(
        path + '/index_classification.tif')

    return water_occurrence, min_extent, max_extent, classification
                                         

def _water_extent(connection, month_start, month_end, start, end, geometry, region, threshold, cloud_cover, rgb_processing): 
        
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
    gdf = gpd.read_file(geometry)
    bbox = gdf.geometry.total_bounds    
    spatial_extent = {'west': bbox[0], 'east': bbox[2], 'south': bbox[1], 'north': bbox[3], 'crs': 4326} 
    start_date_exclusion = (month_start + relativedelta(months=-1)) 

    # Loading S2 collection
    s2_cube = connection.load_collection(
        'SENTINEL2_L2A_SENTINELHUB',
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=['B02', 'B03', 'B04', 'B08', 'sunAzimuthAngles', 'sunZenithAngles'],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover}
    )
 
      # Loading S2 collection for cloud masking processing
    s2_cube_masking = connection.load_collection(
        'SENTINEL2_L2A_SENTINELHUB',
        spatial_extent=spatial_extent,
        temporal_extent=[start_date_exclusion, month_end],
        bands=['CLP', 'SCL'],
        properties={"eo:cloud_cover": lambda v: v <= cloud_cover}
    )

    # Applying mask from sen2cor
    scl = s2_cube_masking.band("SCL")
    mask_scl = (scl == 3) | (scl == 8) | (scl == 9) | (scl == 10) | (scl == 11)
    clp = s2_cube_masking.band("CLP")
    mask_clp = mask_scl | (clp / 255) > 0.3

    # Hill-shade masking
    dem_cube = connection.load_collection("COPERNICUS_30",
                                          spatial_extent=spatial_extent,
                                          temporal_extent=["2010-01-01", "2030-12-31"])

    dem_cube = dem_cube.max_time()

    # Resample s2 cube (Azimuth and Zenith) to 30m
    s2_cube_30 = s2_cube.resample_spatial(resolution = 30, method = 'average')
    
    # DEM WGS84 to DEM UTM from S2 cube
    dem_cube_s2 = dem_cube.resample_cube_spatial(s2_cube_30)
    
    # Merge 30m DEM and 30m s2 cube due Azimuth and Zenith
    merged_cube = s2_cube_30.merge_cubes(dem_cube_s2)

    # Apply the hill-shade udf  
    print('udf', os.path.dirname(os.path.abspath(__file__)) + '/udf.py')
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
    hillshade = hillshade.rename_labels("bands", ["hillshade_mask", "B03", "B04", "B08", "sunAzimuthAngles", "sunZenithAngles", "DEM"])
                   
    # Select a hill-shade band from hill-shade cube and resample it to 10m using s2 cube 
    hillshade_mask = hillshade.band("hillshade_mask").resample_cube_spatial(s2_cube)
    
    # Mask s2 cube with a hillshade mask
    s2_cube_hillshade = s2_cube.mask(hillshade_mask) 

    # Mask s2 cube with a cloud mask
    s2_cube = s2_cube_hillshade.mask(mask_clp.resample_cube_spatial(s2_cube_hillshade))
    
    # Replace 0 to nan in s2 cubes
    s2_cube = s2_cube.mask(s2_cube.apply(lambda x: x.eq(0)), replacement = None)
    
    # Adding NDVI and NDWI to datacube 
    s2_cube = append_indices(s2_cube, ["NDWI", "NDVI"]) 
    s2_cube = s2_cube.rename_labels("bands", ["B02", "B03", "B04", "B08", "sunAzimuthAngles", "sunZenithAngles", "NDWI", "NDVI"]) 

    # Create water extent for the S2 collection using logistic expressions from Lookup Table
    def water_function(data):
        return LOOKUPTABLE[region]["S2"](ndwi=data[6], ndvi=data[7])

    s2_cube_water = s2_cube.reduce_dimension(reducer=water_function, dimension="bands")
    s2_cube_water = s2_cube_water.add_dimension("bands", "water_prob", type="bands")

    # Generate a binary water mask by applying a threshold to the water probability band
    s2_cube_water_threshold = s2_cube_water.apply_dimension(dimension="bands", process=lambda x: if_(x > 0.75, x, 0))
    s2_cube_water_threshold = s2_cube_water_threshold.rename_labels("bands", ["w_T75"])

    # Calculate the water frequency
    s2_cube_water_sum = s2_cube_water_threshold.reduce_dimension(reducer="sum", dimension="t")
    s2_cube_water_sum = s2_cube_water_sum.rename_labels("bands", ["sum"])

    s2_count = s2_cube.band("B08")
    s2_count = s2_count.reduce_dimension(reducer=lambda data: data.count(), dimension="t")
    s2_cube_swf = s2_cube_water_sum.resample_cube_spatial(s2_count) / s2_count
    s2_cube_swf = s2_cube_swf.rename_labels("bands", ["swf"])
    
    # Monthly median s2 image
    s2_median_water = s2_cube_water.filter_temporal([month_start, month_end]).median_time()
    s2_cube_median = s2_cube.filter_temporal([month_start, month_end]).median_time()
    
    # Loading S1 collection
    s1_cube = connection.load_collection(
        'SENTINEL1_GRD',
        spatial_extent=spatial_extent,
        temporal_extent=[month_start, month_end],
        bands=['VH', 'VV'],
        properties={"polarization": lambda p: p == "DV"}
    )

    # Apply terrain correction for back-scatter values
    s1_cube = s1_cube.sar_backscatter(coefficient="gamma0-terrain", mask=True, elevation_model="COPERNICUS_30")
    s1_cube = s1_cube.rename_labels("bands", ["VH", "VV", "mask", "incidence_angle"])
    s1_cube_mask = s1_cube.band("mask")

    def apply_mask(bands):
        return if_(bands.array_element(2) != 2, bands)

    s1_cube = s1_cube.apply_dimension(apply_mask, dimension="bands")

    # Normalized radar back-scatter
    def log_(x):
        return 10 * log(x, 10)

    # Calculate S1 median mosaic
    s1_median = s1_cube.median_time().apply(log_)

    # Create a water extent for the S1 collection using logistic expression
    def s1_water_function(data):
        return LOOKUPTABLE[region]["S1"](vh=data[0], vv=data[1])

    s1_median_water = s1_median.reduce_dimension(reducer=s1_water_function, dimension="bands")
    exclusion_mask = (s1_median_water.resample_cube_spatial(s2_cube_swf) > 0.5) & (s2_cube_swf < 0.33)
    s1_median_water_mask = s1_median_water.mask(exclusion_mask.resample_cube_spatial(s1_median_water))


    # Calculate water extent for the S1 & S2 collection using logistic expression from the lookup table
    def s1_s2_water_function(data):
        return LOOKUPTABLE[region]["S1_S2"](vv=data[0], ndvi=data[1], ndwi=data[2])

    s1_s2_cube = (
        s1_median.filter_bands(["VV"])
        .resample_cube_spatial(s2_cube_median)
        .merge_cubes(s2_cube_median.filter_bands(["NDVI", "NDWI"]))
    )
    s1_s2_water = (
        s1_s2_cube.reduce_dimension(reducer=s1_s2_water_function, dimension="bands")
        .add_dimension("bands", "var", type="bands")
    )

    # Composites water probabilities according to data availability and accuracy.
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

    # Mask built-up area using ESA world cover layer
    worldcover_cube = connection.load_collection(
        "ESA_WORLDCOVER_10M_2020_V1",
        temporal_extent=['2020-12-30', '2021-01-01'],
        spatial_extent=spatial_extent,
        bands=["MAP"]
    )

    builtup_mask = worldcover_cube.band("MAP") == 50
    water_probability = merge_all.mask(builtup_mask.max_time().resample_cube_spatial(merge_all))
    water_probability = water_probability.rename_labels("bands", ["water_prob_sum"])

    # Apply custom threshold to water probability layer
    output = water_probability > (threshold / 100)
    output = output.rename_labels("bands", ["surface_water"])
    output = output * 1.0

    # Output folder
    region_naming = '_'.join(region.split(" "))
    output_folder = os.path.dirname(geometry) + '/' + f'WWT_{region_naming}_{start.strftime("%Y_%m")}_{end.strftime("%Y_%m")}'

    # Send jobs to the backend
    print('Water mask processing between: ' + str(month_start.strftime("%Y_%m")) + ' and ' + str(month_end.strftime("%Y_%m")))
    my_job = output.create_job(
        title='water_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m")),
        out_format="GTiff",
        job_options={"node_caching": True},
        **{'filename_prefix': 'water_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m"))}
    )
    results = my_job.start_and_wait().get_results()
    results.download_files(output_folder)

    if rgb_processing:
        # Process S2 median image if rgb_processing is True
        print('S2 median image processing between: ' + str(month_start.strftime("%Y_%m")) + ' and ' + str(month_end.strftime("%Y_%m")))
        my_job = s2_cube_median.create_job(
            title='median_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m")),
            out_format="GTiff",
            job_options={"node_caching": True},
            **{'filename_prefix': 'median_' + str(month_start.strftime("%Y_%m")) + '_' + str(month_end.strftime("%Y_%m"))}
    )
    results = my_job.start_and_wait().get_results()
    results.download_files(output_folder)

    return output_folder


def main(backend, start, end, region, geometry, rgb_processing, cloud_cover, threshold):
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
        output_folder = _water_extent(connection, month_start, month_end, start, end, geometry, region, threshold, cloud_cover, rgb_processing)

    # Calculate water frequency (0-1), minimum water extent, maximum water extent, and water classification
    _water_indicators(output_folder)

    print('Successfully finished! The output files are located at:', output_folder)

 
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:21:38 2022

@author: rmgu
"""
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))

import rioxarray as rxr
from worldwater_toolbox.prepare_terrain_shadow_mask import terrain_shadow_mask

dem_path = "/home/ubuntu/senet_python_packages/worldwater-toolbox/test_data/Copernicus_DSM_10_N47_00_E007_00_DEM.dt2"
mask_path = "/home/ubuntu/senet_python_packages/worldwater-toolbox/test_data/Copernicus_DSM_10_N47_00_E007_00_MASK.tif"
sun_zenith = 40.0
sun_azimuth = 120.0


dem = rxr.open_rasterio(dem_path, chunks=100)

# Make sure we work in metric projection
if dem.rio.crs != dem.rio.estimate_utm_crs():
    original_crs = dem.rio.crs
    dem = dem.rio.reproject(dem.rio.estimate_utm_crs())
else:
    original_crs = None

transform = dem.rio.transform()
resolution = dem.rio.resolution()
shadow_mask = terrain_shadow_mask(dem, sun_zenith, sun_azimuth, transform, resolution)

# Reproject back to original projection
if original_crs is not None:
    shadow_mask = shadow_mask.where(dem != dem.rio.nodata, 255)
    shadow_mask = shadow_mask.rio.set_nodata(255)
    shadow_mask = shadow_mask.rio.reproject(original_crs)

shadow_mask.rio.to_raster(mask_path)

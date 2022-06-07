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

dem_path = ""
# Data cube must contain sunAzimuthAngles and sunZenithAngles datasets
s2_cube_path = ""
mask_path = ""


dem = rxr.open_rasterio(dem_path, chunks=100)
s2_cube = rxr.open_rasterio(s2_cube_path, chunks=100)
sun_azimuth = s2_cube["sunAzimuthAngles"]
sun_zenith = s2_cube["sunZenithAngles"]

# Make sure we work in S2 image projection and subset
if dem.rio.crs != s2_cube.rio.crs:
    original_crs = dem.rio.crs
    dem_data = dem.rio.reproject(s2_cube.rio.crs,
                                 shape=s2_cube.rio.shape,
                                 transform=s2_cube.rio.transform())
    new_dem_path = dem_path[:-4]+"_reprojected.tif"
    dem_data.rio.to_raster(new_dem_path)
    dem_data.close()
    dem.close()
    dem = rxr.open_rasterio(new_dem_path, chunks=100)
else:
    original_crs = None

transform = dem.rio.transform()
resolution = dem.rio.resolution()
shadow_mask = terrain_shadow_mask(dem, sun_zenith, sun_azimuth, transform, resolution)

shadow_mask.rio.to_raster(mask_path)

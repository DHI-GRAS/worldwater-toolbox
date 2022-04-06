# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:30:39 2022

@author: rmgu
"""

import xarray as xr
import numpy as np

from hillshade.hillshade import hillshade
from hillshade.scripts.cli import rasterize


def _run_shader(elevation_model, resolution_x, resolution_y, zenith, ray_xdir, ray_ydir):

    # Assume chunking is already done by Dask
    ystart = 0
    yend = elevation_model.shape[1]
    # Make sure inputs have the right data type
    resolution = (float(resolution_x), float(resolution_y))
    zenith = float(zenith)
    ray = (float(ray_xdir), float(ray_ydir))

    shadow = hillshade(elevation_model[0, ].astype(np.float32),
                       resolution,
                       zenith,
                       ray,
                       ystart,
                       yend)
    return shadow.astype(np.int8).reshape(elevation_model.shape)


def terrain_shadow_mask(dem, sun_zenith, sun_azimuth, dem_transform, dem_resolution):
    """Calculate shaded regions based on the elevation model and the incident angles of the sun.
    Params:
        dem (xr.DataArray):
            Two-dimensional array specifying the elevation (in m) at each point of the grid.
        sun_zenith (float):
            Sun zenith angle in degrees
        sun_azimuth (float):
            Sun azimuth in degrees (north is 0, increasing clockwise)
        transform (rasterio.affine.Affine):
            Affine transform of the dem
        resolution (tuple):
            Resolution in meters of the dem
    Returns:
        shadow (xr.DataArray):
            An array of ones where there is shadow and zeros otherwise
    """

    # Convert the azimuth into its components on the XY-plane
    ray_xdir, ray_ydir = rasterize(sun_azimuth, dem_transform)

    shadow = xr.apply_ufunc(_run_shader,
                            dem,
                            dem_resolution[0],
                            dem_resolution[1],
                            sun_zenith,
                            ray_xdir,
                            ray_ydir,
                            dask="parallelized",
                            output_dtypes=[int])

    return shadow

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:30:39 2022

@author: rmgu
"""

import dask.array as da
import numpy as np

from hillshade.hillshade import hillshade
from hillshade.scripts.cli import rasterize


def _run_shader(sun_zenith, sun_azimuth, elevation_model, resolution_x, resolution_y, dem_transform):

    shadow_stack = np.zeros(sun_azimuth.shape)

    for step in range(sun_azimuth.shape[0]):
        azimuth = np.nanmean(sun_azimuth[step, ...])
        zenith = np.nanmean(sun_zenith[step, ...])
        if np.all(np.isnan(zenith)):
            continue
        # Convert the azimuth into its components on the XY-plane
        ray_xdir, ray_ydir = rasterize(azimuth, dem_transform)

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
        shadow_stack[step, ...] = shadow.astype(np.int16).reshape(elevation_model.shape)

    shadow_stack[np.isnan(sun_azimuth)] = 255
    return shadow_stack


def terrain_shadow_mask(dem, sun_zenith, sun_azimuth, dem_transform, dem_resolution):
    """Calculate shaded regions based on the elevation model and the incident angles of the sun.
    Params:
        dem (xr.DataArray):
            Two-dimensional array specifying the elevation (in m) at each point of the grid.
        sun_zenith (xr.DataArray):
            Three-dimensional array (time, x, y) specifying the sun zenith angle (in degrees)
        sun_azimuth (float):
            Three-dimensional array (time, x, y) specifying the sun zenith angle
            (in degrees, north is 0, increasing clockwise)
        transform (rasterio.affine.Affine):
            Affine transform of the dem
        resolution (tuple):
            Resolution in meters of the dem
    Returns:
        shadow (xr.DataArray):
            An array of ones where there is shadow and zeros otherwise
    """

    dem_overlapped = da.overlap.overlap(dem.data, depth={0: 0, 1: 10, 2: 10}, boundary=0)
    sun_azimuth_overlapped = da.overlap.overlap(sun_azimuth.data,
                                                depth={0: 0, 1: 10, 2: 10},
                                                boundary=0)
    sun_zenith_overlapped = da.overlap.overlap(sun_zenith.data,
                                               depth={0: 0, 1: 10, 2: 10},
                                               boundary=0)
    shadow_overlapped = da.map_blocks(_run_shader,
                                      sun_zenith_overlapped,
                                      sun_azimuth_overlapped,
                                      dem_overlapped,
                                      dem_resolution[0],
                                      dem_resolution[1],
                                      dem_transform)
    shadow_da = da.overlap.trim_internal(shadow_overlapped, {0: 0, 1: 10, 2: 10})
    shadow = sun_zenith.copy(data=shadow_da).astype(np.int16)

    shadow = shadow.rio.set_nodata(255)
    return shadow

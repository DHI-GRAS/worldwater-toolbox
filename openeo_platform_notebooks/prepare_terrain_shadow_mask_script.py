# -*- coding: utf-8 -*-
"""
Created on Mon Apr4 15:30:39 2022

@author: rmgu
"""

import dask.array as da
import numpy as np
from numba import jit

#from hillshade.hillshade import hillshade
#from hillshade.scripts.cli import rasterize

def rasterize(azimuth, transform=None):
    """Convert the azimuth into its components on the XY-plane. Depending on the value of the
    azimuth either the x or the y component of the resulting vector is scaled to 1, so that
    it can be used conveniently to walk a grid.
    """
    azimuth = np.deg2rad(azimuth)
    xdir, ydir = np.sin(azimuth), np.cos(azimuth)

    if transform is not None:
        xdir, ydir = transform * (xdir, ydir)
        xdir -= transform.xoff
        ydir -= transform.yoff

    slope = ydir / xdir
    if slope < 1. and slope > -1.:
        xdir = 1.
        ydir = slope
    else:
        xdir = 1. / slope
        ydir = 1.
    return xdir, ydir

@jit('f8(f8, f8, f8)', nopython=True, nogil=True)
def height(xbase, ybase, angle):
    """Calculate the height of a right triangle whose base is defined by a vector in the xy-plane.
    The offset is added to the triangle height.
    """
    rho = (xbase**2 + ybase**2)**.5
    return rho * np.tan(angle)


@jit('boolean(f8, f8, Tuple((i8, i8)))', nopython=True, nogil=True)
def within_bounds(pixel_x, pixel_y, bounds):
    """Check whether x and y pixel coordinates are within bounds"""
    return (pixel_y < bounds[0]) and (pixel_x < bounds[1]) and (pixel_y > 0) and (pixel_x > 0)


@jit('i8[:,:](f4[:,:], Tuple((f8,f8)), f8, Tuple((f8, f8)), i8, i8)', nopython=True, nogil=True)
def hillshade(elevation_model, resolution, zenith, ray, ystart, yend):
    """Calculate a shaded region for elevation_model[ystart:yend] by looping over every
    pixel of the elevation model and tracing the path towards the sun until an obstacle is
    hit or the maximum elevation of the model is reached. The path is defined by a rasterized
    direction of the sun ray in the xy-plane (ray) and the zenith.
    Params:
        elevation_model (np.ndarray):
            Two-dimensional array specifying the elevation at each point of the grid.
        resolution (tuple):
            resolution in meters of the elevation_model
        zenith (float):
            zenith in degrees
        ray (tuple):
            rasterized XY-direction of azimuth
        ystart (int):
            y-chunk starting index
        yend (int):
            y-chunk ending index
    Returns:
        shadow (np.ndarray):
            an array of ones where there is shade and zeros otherwise
    """
    if max(ray) != 1.:
        raise ValueError("xy-direction is not rasterized.")
    shadow = np.zeros((yend - ystart, elevation_model.shape[1]), dtype=np.int64)
    zenith = np.deg2rad(90 - zenith)
    dx, dy = ray
    xres, yres = resolution
    z_max = elevation_model.max()
    bounds = elevation_model.shape

    for pixel_y in range(ystart, yend):
        for pixel_x in range(elevation_model.shape[1]):

            pixel_z = elevation_model[pixel_y, pixel_x]
            ray_x = float(pixel_x)
            ray_y = float(pixel_y)
            intersection = None

            while within_bounds(ray_x, ray_y, bounds):
                xbase = (ray_x - pixel_x) * xres
                ybase = (ray_y - pixel_y) * yres
                ray_z = height(xbase, ybase, zenith) + pixel_z
                if ray_z > z_max:
                    break
                if ray_z < elevation_model[int(ray_y), int(ray_x)]:
                    intersection = (ray_y, ray_x)
                    break
                ray_x += dx
                ray_y += dy

            if intersection is not None:
                shadow[pixel_y - ystart, pixel_x] = 1
    return shadow


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
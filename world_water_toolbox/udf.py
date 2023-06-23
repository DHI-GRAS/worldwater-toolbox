import numpy as np
from openeo.udf import XarrayDataCube
from openeo.udf.debug import inspect
from hillshade.hillshade import hillshade


def rasterize(azimuth, resolution=None):
    azimuth = np.deg2rad(azimuth)
    xdir, ydir = np.sin(azimuth), np.cos(azimuth)

    if resolution is not None:
        xdir = xdir * resolution[0]
        ydir = ydir * resolution[1]
        signx = np.sign(xdir)
        signy = np.sign(ydir)
    slope = abs(ydir / xdir)

    if slope < 1. and slope > -1.:
        xdir = 1.
        ydir = slope
    else:
        xdir = 1. / slope
        ydir = 1.

    return xdir * signx, ydir * signx


def _run_shader(sun_zenith, sun_azimuth, elevation_model, resolution_x, resolution_y):
    azimuth = np.nanmean(sun_azimuth.astype(np.float32))
    zenith = np.nanmean(sun_zenith.astype(np.float32))
    if np.isnan(azimuth):
        shadow = np.zeros(elevation_model.shape) + 255
    else:
        resolution = (float(resolution_x), float(resolution_y))
        ray_xdir, ray_ydir = rasterize(azimuth, resolution)

        # Assume chunking is already done by Dask
        ystart = 0
        yend = elevation_model.shape[0]

        # Make sure inputs have the right data type
        zenith = float(zenith)
        ray = (float(ray_xdir), float(ray_ydir))
        shadow = hillshade(
            elevation_model.astype(np.float32),
            resolution,
            zenith,
            ray,
            ystart,
            yend
        )
        shadow = shadow.reshape(elevation_model.shape)
        shadow[np.isnan(sun_azimuth)] = 255
    return shadow


def apply_datacube(cube: XarrayDataCube, context: dict) -> XarrayDataCube:
    in_xarray = cube.get_array()
    sun_zenith = in_xarray.sel({"bands": "sunZenithAngles"}).values.astype(np.float32)
    sun_azimuth = in_xarray.sel({"bands": "sunAzimuthAngles"}).values.astype(np.float32)
    elevation_model = in_xarray.sel({"bands": "DEM"}).values.astype(np.float32)
    res_y = in_xarray.coords["y"][int(len(in_xarray.coords["y"]) / 2) + 1] - in_xarray.coords["y"][
        int(len(in_xarray.coords["y"]) / 2)]
    res_x = in_xarray.coords["x"][int(len(in_xarray.coords["x"]) / 2) + 1] - in_xarray.coords["x"][
        int(len(in_xarray.coords["x"]) / 2)]

    sun_zenith = sun_zenith

    shadow = _run_shader(sun_zenith, sun_azimuth, elevation_model, res_x, res_x)
    cube.get_array().values[0] = shadow

    return cube

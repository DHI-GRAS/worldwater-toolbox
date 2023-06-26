import numpy as np
from openeo.udf import XarrayDataCube
from openeo.udf.debug import inspect
from hillshade.hillshade import hillshade


def rasterize(azimuth, resolution=None):
    """
    Rasterize the azimuth angles from Sentinel-2 metadata.
    
    Parameters
    ----------
    azimuth : numpy array
        Azimuth angles from Sentinel-2 metadata.

    resolution : tuple, optional
        Resolution of the azimuth file, e.g., [10, 10].

    Returns
    -------
    water_frequency : numpy array
        Water frequency as an integer.
    """
    
    # Convert azimuth angles to radians
    azimuth = np.deg2rad(azimuth)
    
    # Calculate x and y direction components using sine and cosine
    xdir, ydir = np.sin(azimuth), np.cos(azimuth)

    # Adjust x and y direction components based on the resolution
    if resolution is not None:
        xdir = xdir * resolution[0]
        ydir = ydir * resolution[1]
        signx = np.sign(xdir)
        signy = np.sign(ydir)
    
    # Calculate the slope of the azimuth angles
    slope = abs(ydir / xdir)

    # Adjust x and y direction components based on the slope
    if slope < 1. and slope > -1.:
        xdir = 1.
        ydir = slope
    else:
        xdir = 1. / slope
        ydir = 1.

    return xdir * signx, ydir * signx



def _run_shader(sun_zenith, sun_azimuth, elevation_model, resolution_x, resolution_y):
    """
    Calculate shaded regions based on the elevation model and the incident angles of the sun.

    Parameters:
    -----------
    sun_zenith : xr.DataArray
        Three-dimensional array (time, x, y) specifying the sun zenith angle (in degrees).

    sun_azimuth : float
        Three-dimensional array (time, x, y) specifying the sun zenith angle (in degrees, north is 0, increasing clockwise)

    elevation_model : xr.DataArray
        Two-dimensional array specifying the elevation (in m) at each point of the grid.

    resolution_x : float
        Resolution in meters of the dem in the x-direction.

    resolution_y : float
        Resolution in meters of the dem in the y-direction.

    Returns
    -------
    shadow : xr.DataArray
        An array of ones where there is shadow and zeros otherwise.
    """
    
    # Calculate mean azimuth and zenith
    azimuth = np.nanmean(sun_azimuth.astype(np.float32))
    zenith = np.nanmean(sun_zenith.astype(np.float32))

    if np.isnan(azimuth):
        # If azimuth is NaN, set shadow to 255
        shadow = np.zeros(elevation_model.shape) + 255
    else:
        # Calculate resolution
        resolution = (float(resolution_x), float(resolution_y))
        
        # Rasterize azimuth to obtain ray directions
        ray_xdir, ray_ydir = rasterize(azimuth, resolution)

        # Assume chunking is already done by Dask
        ystart = 0
        yend = elevation_model.shape[0]

        # Make sure inputs have the right data type
        zenith = float(zenith)
        ray = (float(ray_xdir), float(ray_ydir))
        
        # Compute hill shade using elevation model and incident angles
        shadow = hillshade(
            elevation_model.astype(np.float32),
            resolution,
            zenith,
            ray,
            ystart,
            yend
        )
        
        # Reshape shadow to match elevation model shape
        shadow = shadow.reshape(elevation_model.shape)
        
        # Set shadow to 255 where sun azimuth is NaN
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


    shadow = _run_shader(sun_zenith, sun_azimuth, elevation_model, res_x, res_x)
    cube.get_array().values[0] = shadow

    return cube

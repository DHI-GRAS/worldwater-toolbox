import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter


def xr_filter(da: xr.DataArray, filter_func, **kwargs) -> xr.DataArray:
    """
    Apply a scipy spatial filter in a NaN-safe manner.

    Pixels that are NaN in the input remain NaN in the output.

    Parameters
    ----------
    da : xr.DataArray
        2-D input slice.
    filter_func : callable
        A scipy.ndimage filter (e.g. ``gaussian_filter``, ``uniform_filter``).
    **kwargs
        Passed directly to *filter_func*.

    Returns
    -------
    xr.DataArray
    """
    a = da.values
    mask = np.isfinite(a)
    a_filled = np.where(mask, a, 0.0)

    filtered_values = filter_func(a_filled, **kwargs)
    filtered_mask = filter_func(mask.astype(float), **kwargs)

    with np.errstate(invalid="ignore", divide="ignore"):
        result = filtered_values / filtered_mask

    result[~mask] = np.nan
    return xr.DataArray(result, coords=da.coords, dims=da.dims)


def median_smooth(swf: xr.DataArray, size: int = 3) -> xr.DataArray:
    """
    Apply a median filter; NaN positions are preserved.

    Parameters
    ----------
    swf : xr.DataArray
        2-D input slice.
    size : int
        Kernel size (default 3 -> 3x3).

    Returns
    -------
    xr.DataArray
    """
    a = swf.values
    smoothed = median_filter(a, size=size)
    smoothed[~np.isfinite(a)] = np.nan
    return xr.DataArray(smoothed, coords=swf.coords, dims=swf.dims)


def gaussian_smooth(swf: xr.DataArray, sigma: float = 0.8) -> xr.DataArray:
    """
    Apply a Gaussian filter; NaN positions are preserved.

    Parameters
    ----------
    swf : xr.DataArray
        2-D input slice.
    sigma : float
        Standard deviation for the Gaussian kernel (default 0.8).

    Returns
    -------
    xr.DataArray
    """
    return xr_filter(swf, gaussian_filter, sigma=sigma)


def uniform_smooth(swf: xr.DataArray, size: int = 3) -> xr.DataArray:
    """
    Apply a uniform (box) filter; NaN positions are preserved.

    Parameters
    ----------
    swf : xr.DataArray
        2-D input slice.
    size : int
        Kernel size (default 3 -> 3x3).

    Returns
    -------
    xr.DataArray
    """
    return xr_filter(swf, uniform_filter, size=size, mode="nearest")


def smooth_timeseries(
    water_ts: xr.DataArray,
    method: str = "gaussian",
    sigma: float = 0.8,
    filter_size: int = 3,
) -> xr.DataArray:
    """
    Apply spatial smoothing to every time slice of *water_ts*.

    Parameters
    ----------
    water_ts : xr.DataArray
        Dims: (time, y, x).
    method : {"gaussian", "median", "uniform"}
    sigma : float
        Gaussian sigma (used only when method="gaussian").
    filter_size : int
        Kernel size for median/uniform filters.

    Returns
    -------
    xr.DataArray
    """
    if method not in ("gaussian", "median", "uniform"):
        raise ValueError(f"Unknown smooth method: {method!r}. Choose gaussian/median/uniform.")

    slices = []
    for t in water_ts.time:
        sl = water_ts.sel(time=t)
        if method == "gaussian":
            slices.append(gaussian_smooth(sl, sigma=sigma))
        elif method == "median":
            slices.append(median_smooth(sl, size=filter_size))
        elif method == "uniform":
            slices.append(uniform_smooth(sl, size=filter_size))

    return xr.concat(slices, dim="time")

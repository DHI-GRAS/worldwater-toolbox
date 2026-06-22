"""
openeo_filters.py
Author: Walid Ghariani
Date: 19-05-2026
Description: OpenEO spatial filters and hysteresis thresholding for water probability cubes.
"""

from typing import Union

import numpy as np
from scipy.ndimage import gaussian_filter

import openeo
from openeo import DataCube
from openeo.api.process import Parameter
from openeo.processes import if_, eq


FILTER_TYPES = ["none", "gaussian", "uniform", "median"]

HYSTERESIS_DEFAULTS = {
    "high_thr":    0.6,
    "low_thr":     0.4,
    "max_dist_m":  250.0,
    "connectivity": 2,
}

# Tile / overlap for apply_neighborhood calls.
# overlap_px must be >= ceil(max_dist_m / pixel_size_m):
# 250 m / 10 m = 25 px -> 32 for a small safety margin. 
# PS: we might need to check/test this further
_MEDIAN_TILE_PX    = 128
_MEDIAN_OVERLAP_PX = 8
_HYST_TILE_PX      = 1024
_HYST_OVERLAP_PX   = 32


def _gaussian_kernel(sigma: float = 0.85, truncate: float = 4.0) -> list:
    """Build a 2-D Gaussian kernel as a nested list for apply_kernel."""
    radius = int(truncate * sigma + 0.5)
    size = 2 * radius + 1
    impulse = np.zeros((size, size))
    impulse[radius, radius] = 1.0
    k = gaussian_filter(impulse, sigma=sigma)
    return (k / k.sum()).tolist()


def _uniform_kernel(size: int = 3) -> list:
    """Build a uniform (box) kernel as a nested list for apply_kernel."""
    return (np.ones((size, size)) / size ** 2).tolist()


_MEDIAN_UDF = """
from openeo.udf import XarrayDataCube
import numpy as np
from scipy.ndimage import median_filter

def apply_datacube(cube: XarrayDataCube, context: dict) -> XarrayDataCube:
    array = cube.get_array()
    size = context.get("size", 3)
    values = array.values.copy()
    smoothed = median_filter(values, size=size)
    smoothed[~np.isfinite(values)] = np.nan
    return XarrayDataCube(array.copy(data=smoothed))
"""

_HYSTERESIS_UDF = """
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure
from openeo.udf import XarrayDataCube


def _pixel_size_m(da):
    dx, dy = 10.0, 10.0  # fallback: 10 m (Sentinel-2 native)
    for xdim in ("x", "lon"):
        if xdim in da.coords and da.sizes.get(xdim, 0) > 1:
            dx = abs(float(da.coords[xdim].values[1] - da.coords[xdim].values[0]))
            break
    for ydim in ("y", "lat"):
        if ydim in da.coords and da.sizes.get(ydim, 0) > 1:
            dy = abs(float(da.coords[ydim].values[1] - da.coords[ydim].values[0]))
            break
    return dx, dy


def _hyst_2d(arr2d, max_steps, struct, high_thr, low_thr):
    valid   = np.isfinite(arr2d)
    seeds   = (arr2d >= high_thr) & valid
    allowed = (arr2d >= low_thr)  & valid
    if seeds.sum() == 0:
        return np.zeros_like(arr2d, dtype=np.float32)
    current = seeds.copy()
    for _ in range(max_steps):
        grown = binary_dilation(current, structure=struct)
        new   = grown & allowed & ~current
        if not np.any(new):
            break
        current |= new
    return current.astype(np.float32)


def apply_datacube(cube: XarrayDataCube, context: dict) -> XarrayDataCube:
    high_thr     = float(context.get("high_thr",    0.6))
    low_thr      = float(context.get("low_thr",     0.4))
    max_dist_m   = float(context.get("max_dist_m",  250.0))
    connectivity =  int(context.get("connectivity", 2))

    da = cube.get_array()
    has_bands = "bands" in da.dims
    if has_bands:
        da = da.squeeze("bands", drop=True)

    dx, dy = _pixel_size_m(da)
    max_steps = int(np.ceil(max_dist_m / min(dx, dy)))
    struct = generate_binary_structure(2, connectivity)

    tdim = next((d for d in ("t", "time") if d in da.dims), None)
    if tdim is None:
        result = _hyst_2d(da.values, max_steps, struct, high_thr, low_thr)
        result_da = da.copy(data=result)
    else:
        slices = [
            _hyst_2d(da.isel({tdim: i}).values, max_steps, struct, high_thr, low_thr)
            for i in range(da.sizes[tdim])
        ]
        result_da = da.copy(data=np.stack(slices, axis=0))

    if has_bands:
        result_da = result_da.expand_dims("bands")
    return XarrayDataCube(result_da)
"""


def apply_spatial_filter(
    cube: DataCube,
    filter_type: Union[str, Parameter] = "gaussian",
    gaussian_sigma: float = 0.85,
    uniform_size: int = 3,
    median_size: int = 3,
    median_tile_px: int = _MEDIAN_TILE_PX,
    median_overlap_px: int = _MEDIAN_OVERLAP_PX,
) -> DataCube:
    """
    Apply a spatial smoothing filter to *cube*.

    Parameters
    ----------
    cube : DataCube
        Input water probability cube (bands preserved as-is).
    filter_type : str or Parameter
        One of "none", "gaussian", "uniform", "median".
        Pass an openeo.Parameter for UDP compatibility - all three filter
        branches are wired into the process graph and the backend selects at
        runtime via if_() nodes.
    gaussian_sigma : float
        Gaussian kernel std-dev (default 0.85).
    uniform_size : int
        Box kernel side length in pixels (default 3 -> 3x3).
    median_size : int
        Median window size (default 3 -> 3x3).
    median_tile_px : int
        apply_neighborhood tile size for the median path.
    median_overlap_px : int
        apply_neighborhood overlap for the median path.

    Returns
    -------
    DataCube
        Smoothed cube with the same band structure as *cube*.
        Returns *cube* unchanged when filter_type is "none".
    """
    gaussian_cube = cube.apply_kernel(kernel=_gaussian_kernel(sigma=gaussian_sigma))
    uniform_cube  = cube.apply_kernel(kernel=_uniform_kernel(size=uniform_size))
    median_cube   = cube.apply_neighborhood(
        process=openeo.UDF(_MEDIAN_UDF, runtime="Python", context={"size": median_size}),
        size=[
            {"dimension": "x", "value": median_tile_px,    "unit": "px"},
            {"dimension": "y", "value": median_tile_px,    "unit": "px"},
        ],
        overlap=[
            {"dimension": "x", "value": median_overlap_px, "unit": "px"},
            {"dimension": "y", "value": median_overlap_px, "unit": "px"},
        ],
    )

    if isinstance(filter_type, str):
        choices = {
            "none":     cube,
            "gaussian": gaussian_cube,
            "uniform":  uniform_cube,
            "median":   median_cube,
        }
        if filter_type not in choices:
            raise ValueError(
                f"filter_type must be one of {FILTER_TYPES}, got {filter_type!r}"
            )
        return choices[filter_type]

    result_raw = if_(
        eq(filter_type, "none"), cube,
        if_(
            eq(filter_type, "gaussian"), gaussian_cube,
            if_(eq(filter_type, "uniform"), uniform_cube, median_cube),
        ),
    )
    return DataCube(result_raw.pgnode, cube.connection, metadata=cube.metadata)


def apply_hysteresis(
    cube: DataCube,
    high_thr: Union[float, Parameter]  = HYSTERESIS_DEFAULTS["high_thr"],
    low_thr: Union[float, Parameter]   = HYSTERESIS_DEFAULTS["low_thr"],
    max_dist_m: Union[float, Parameter] = HYSTERESIS_DEFAULTS["max_dist_m"],
    connectivity: Union[int, Parameter] = HYSTERESIS_DEFAULTS["connectivity"],
    tile_px: int    = _HYST_TILE_PX,
    overlap_px: int = _HYST_OVERLAP_PX,
) -> DataCube:
    """
    Apply hysteresis (double) thresholding to a smoothed probability cube.

    Seeds are pixels >= high_thr.  Connected neighbours within max_dist_m that
    are >= low_thr are retained.  This fills in weak-signal water adjacent to
    confident water without pulling in isolated noise.

    Input values must be in [0, 1] — matching the raw logistic-regression
    water probability output from _water_probability().

    Parameters
    ----------
    cube : DataCube
        Smoothed water probability in [0, 1].  A single-band cube is expected;
        the UDF squeezes the "bands" dimension if present.
    high_thr : float or Parameter
        Seed threshold (confident water).  Default 0.6.
    low_thr : float or Parameter
        Extension threshold (candidate water adjacent to seeds).  Default 0.4.
    max_dist_m : float or Parameter
        Maximum grow distance in metres.  Default 250 m.
        overlap_px must be >= ceil(max_dist_m / pixel_size_m).
    connectivity : int or Parameter
        Binary structure connectivity: 1=cross (4-connected), 2=full 3x3
        (8-connected, default).
    tile_px : int
        apply_neighborhood tile size.  Increase for large AOIs to avoid
        water bodies being cut at tile boundaries.
    overlap_px : int
        apply_neighborhood overlap.  Must be >= ceil(max_dist_m / pixel_size_m);
        at 10 m / 250 m that is 25 px - default 32 gives a small margin.

    Returns
    -------
    DataCube
        Binary water mask (0.0 / 1.0 float32).
    """
    context = {
        "high_thr":    high_thr,
        "low_thr":     low_thr,
        "max_dist_m":  max_dist_m,
        "connectivity": connectivity,
    }
    return cube.apply_neighborhood(
        process=openeo.UDF(_HYSTERESIS_UDF, runtime="Python", context=context),
        size=[
            {"dimension": "x", "value": tile_px,    "unit": "px"},
            {"dimension": "y", "value": tile_px,    "unit": "px"},
        ],
        overlap=[
            {"dimension": "x", "value": overlap_px, "unit": "px"},
            {"dimension": "y", "value": overlap_px, "unit": "px"},
        ],
    )

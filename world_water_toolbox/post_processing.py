"""
post_processing.py
Author: Walid Ghariani
Date: 20-04-2026
Description: water Indicators Post-Processing Pipeline
"""

import argparse
import re
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import rioxarray  
import xarray as xr
from scipy.ndimage import binary_dilation, generate_binary_structure

from .filters import smooth_timeseries

HIGH_THR: float = 0.6
LOW_THR: float = 0.4
MAX_DIST_M: float = 250.0
CONNECTIVITY: int = 2
NODATA: int = 255
SEASONALITY_NODATA: int = 255

COMPRESS: str = "ZSTD"
BLOCKSIZE: int = 256
PREDICTOR: int = 2

_CREATION_OPTIONS = {
    "TILED": "YES",
    "COMPRESS": COMPRESS,
    "BLOCKXSIZE": str(BLOCKSIZE),
    "BLOCKYSIZE": str(BLOCKSIZE),
    "PREDICTOR": str(PREDICTOR),
}

def start_time_from_path(path) -> pd.Timestamp:
    """
    Parse a timestamp from a water mask filename.

    Supports two naming conventions:
      - Input files:  ``water_YYYY_MM_YYYY_MM.tif``
      - Output files: ``<prefix>_YYYYMMDD.tif``

    Parameters
    ----------
    path : str or Path
        File path to parse.

    Returns
    -------
    pd.Timestamp
    """
    p = str(path).replace("\\", "/")
    # Output pattern: _YYYYMMDD.tif
    m = re.search(r"_(\d{8})\.tiff?$", p, flags=re.IGNORECASE)
    if m:
        return pd.to_datetime(m.group(1), format="%Y%m%d")
    # Input pattern: water_YYYY_MM_YYYY_MM.tif
    m = re.search(r"water_(\d{4})_(\d{2})_\d{4}_\d{2}\.tif$", p)
    if m:
        year, month = int(m.group(1)), int(m.group(2))
        return pd.Timestamp(year=year, month=month, day=1)

    raise ValueError(f"Could not parse timestamp from: {path}")


def pixel_size_m(da: xr.DataArray) -> tuple:
    """
    Return the pixel size (x, y) in metres from a DataArray's rio metadata.

    Parameters
    ----------
    da : xr.DataArray

    Returns
    -------
    tuple of float
        (x_resolution, y_resolution) both positive.
    """
    xres, yres = da.rio.resolution()
    return float(abs(xres)), float(abs(yres))


def set_geo_metadata(
    target: Union[xr.DataArray, xr.Dataset],
    source: Union[xr.DataArray, xr.Dataset],
) -> Union[xr.DataArray, xr.Dataset]:
    """
    Copy CRS and affine transform from *source* to *target*.

    Parameters
    ----------
    target : xr.DataArray or xr.Dataset
    source : xr.DataArray or xr.Dataset

    Returns
    -------
    xr.DataArray or xr.Dataset
        *target* with spatial metadata written.
    """
    if hasattr(source, "rio") and hasattr(target, "rio"):
        if source.rio.crs is not None:
            target = target.rio.write_crs(source.rio.crs)
        try:
            if source.rio.transform() is not None:
                target = target.rio.write_transform(source.rio.transform())
        except Exception:
            pass
    return target


def load_monthly_water_maps(
    data_dir,
    pattern: str = "water_*.tif",
) -> xr.DataArray:
    """
    Scan *data_dir* for monthly water probability GeoTIFFs, load them, and
    return a time-stacked DataArray.

    Parameters
    ----------
    data_dir : str or Path
        Directory containing the monthly TIFF files.
    pattern : str
        Glob pattern used to find files (default ``"water_*.tif"``).

    Returns
    -------
    xr.DataArray
        Dims: (time, y, x), sorted by time.
    """
    data_dir = Path(data_dir)
    tif_files = sorted(data_dir.glob(pattern))

    if not tif_files:
        tif_files = sorted(data_dir.rglob(pattern))

    if not tif_files:
        raise FileNotFoundError(
            f"No files matching {pattern!r} found under {data_dir}"
        )

    rasters = []
    for f in tif_files:
        da = rioxarray.open_rasterio(f).sel(band=1)
        t = start_time_from_path(f)
        da = da.expand_dims(time=[t])
        rasters.append(da)

    water_ts = xr.concat(rasters, dim="time").sortby("time")
    water_ts.name = "water_probability"
    print(f"Loaded {len(rasters)} time slices from {data_dir}")
    return water_ts


def hysteresis_timeseries(
    probas: xr.DataArray,
    high_thr: float = HIGH_THR,
    low_thr: float = LOW_THR,
    max_dist_m: float = MAX_DIST_M,
    connectivity: int = CONNECTIVITY,
    extra_mask: Union[xr.DataArray, None] = None,
    out_nodata: int = NODATA,
) -> xr.DataArray:
    """
    Apply double-threshold hysteresis thresholding per time slice.

    Seeds are pixels above *high_thr*; they are grown into adjacent pixels
    above *low_thr* within a maximum Euclidean distance of *max_dist_m*
    metres.  Each month is processed independently.

    Parameters
    ----------
    probas : xr.DataArray
        Dims: (time, y, x).  Values typically in [0, 1] or [0, 100].
    high_thr : float
        High-confidence water threshold (seeds).
    low_thr : float
        Minimum threshold for region growth.
    max_dist_m : float
        Maximum growth radius in metres.
    connectivity : {1, 2}
        Pixel connectivity: 1 = 4-connected, 2 = 8-connected.
    extra_mask : xr.DataArray or None
        Optional spatial mask (True = allowed pixels).
    out_nodata : int
        NoData value written into the output uint8 array.

    Returns
    -------
    xr.DataArray
        uint8 array with values 0 (non-water), 1 (water), *out_nodata* (nodata).
    """
    dx, dy = pixel_size_m(probas.isel(time=0))
    step_m = min(dx, dy)
    max_steps = int(np.ceil(max_dist_m / step_m))
    struct = generate_binary_structure(2, connectivity)

    in_nodata = probas.rio.nodata
    extra_mask_np = extra_mask.values.astype(bool) if extra_mask is not None else None

    def _hyst_np(prob2d):
        valid = np.ones_like(prob2d, dtype=bool)
        if in_nodata is not None:
            valid &= prob2d != in_nodata
        if np.issubdtype(prob2d.dtype, np.floating):
            valid &= np.isfinite(prob2d)

        seeds = (prob2d >= high_thr) & valid
        allowed = (prob2d >= low_thr) & valid

        if extra_mask_np is not None:
            seeds &= extra_mask_np
            allowed &= extra_mask_np

        if seeds.sum() == 0:
            return np.zeros_like(seeds, dtype=np.uint8)

        current = seeds.copy()
        for _ in range(max_steps):
            grown = binary_dilation(current, structure=struct)
            new = grown & allowed & ~current
            if not np.any(new):
                break
            current |= new

        return current.astype(np.uint8)

    spatial_dims = [d for d in probas.dims if d != "time"]
    if len(spatial_dims) != 2:
        raise ValueError(f"Expected dims (time, y, x). Got {probas.dims!r}")

    out_u8 = xr.apply_ufunc(
        _hyst_np,
        probas,
        input_core_dims=[spatial_dims],
        output_core_dims=[spatial_dims],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.uint8],
        keep_attrs=True,
    ).rename("water_hysteresis")

    valid_da = xr.ones_like(probas.isel(time=0), dtype=bool)
    if in_nodata is not None:
        valid_da = valid_da & (probas.isel(time=0) != in_nodata)
    if np.issubdtype(probas.dtype, np.floating):
        valid_da = valid_da & xr.apply_ufunc(np.isfinite, probas.isel(time=0))

    out_u8 = out_u8.where(valid_da, other=np.uint8(out_nodata))
    out_u8 = out_u8.rio.write_nodata(np.uint8(out_nodata), encoded=True)
    out_u8.attrs.update(
        {
            "method": "double_threshold_hysteresis_per_month",
            "high_thr": float(high_thr),
            "low_thr": float(low_thr),
            "max_dist_m": float(max_dist_m),
            "connectivity": int(connectivity),
            "input_nodata": None if in_nodata is None else int(in_nodata),
            "output_nodata": int(out_nodata),
        }
    )
    return out_u8


def compute_annual_classification(
    water_bin: xr.DataArray,
    time_dim: str = "time",
    seasonal_min_months: int = 2,
    seasonal_max_months: int = 9,
    permanent_min_months: int = 10,
) -> tuple:
    """
    Classify pixels as non-water, seasonal, or permanent from monthly
    binary water masks.

    Parameters
    ----------
    water_bin : xr.DataArray
        Monthly binary masks (0/1, NaN for missing).  Dims: (time, y, x).
    time_dim : str
        Name of the time dimension.
    seasonal_min_months : int
        Lower bound (inclusive) for seasonal class.
    seasonal_max_months : int
        Upper bound (inclusive) for seasonal class.
    permanent_min_months : int
        Pixels with water count ≥ this value are classified as permanent.

    Returns
    -------
    water_sum : xr.DataArray (uint16)
        Count of water months per pixel.
    valid_n : xr.DataArray (uint16)
        Count of valid (non-NaN) months per pixel.
    occurrence_pct : xr.DataArray (float32)
        Water occurrence as a percentage of valid months.
    classification : xr.DataArray (uint8)
        0 = non-water, 1 = seasonal, 2 = permanent, 255 = nodata.
    """
    wb = water_bin.where(np.isfinite(water_bin))

    valid = np.isfinite(wb)
    valid_n = valid.sum(time_dim).astype("uint16")

    wb01 = xr.where(valid, (wb > 0).astype("uint8"), 0)
    water_sum = wb01.sum(time_dim).astype("uint16")
    water_sum.name = "water_month_count"

    occurrence_pct = xr.where(valid_n > 0, (water_sum / valid_n) * 100.0, np.nan).astype("float32")
    occurrence_pct.name = "water_occurrence_pct"

    classification = xr.full_like(water_sum, fill_value=0, dtype=np.uint8).rename("water_class")
    seasonal = (water_sum >= seasonal_min_months) & (water_sum <= seasonal_max_months)
    permanent = water_sum >= permanent_min_months

    classification = classification.where(~seasonal, 1)
    classification = classification.where(~permanent, 2)
    classification = classification.where(valid_n > 0, 255).astype(np.uint8)
    classification.attrs.update(
        {
            "classes": {0: "non-water", 1: "seasonal", 2: "permanent", 255: "nodata"},
            "seasonal_range_months": f"{seasonal_min_months}-{seasonal_max_months}",
            "permanent_min_months": permanent_min_months,
        }
    )
    return water_sum, valid_n, occurrence_pct, classification


def compute_minmax_water_extent(
    water_sum: xr.DataArray,
    seasonal_min_months: int = 2,
    permanent_min_months: int = 10,
) -> tuple:
    """
    Derive minimum (permanent) and maximum (seasonal + permanent) water extent
    binary masks from monthly water counts.

    Parameters
    ----------
    water_sum : xr.DataArray (uint16)
        Count of wet months per pixel, as returned by
        ``compute_annual_classification``.
    seasonal_min_months : int
        Pixels wet ≥ this many months are included in the maximum extent.
    permanent_min_months : int
        Pixels wet ≥ this many months are included in the minimum extent.

    Returns
    -------
    min_water : xr.DataArray (uint8)
        1 where water is permanent, 0 elsewhere.
    max_water : xr.DataArray (uint8)
        1 where water is seasonal or permanent, 0 elsewhere.
    """
    max_water = xr.where(water_sum >= seasonal_min_months, 1, 0).astype("uint8")
    max_water.name = "max_water"
    max_water.attrs.update(
        {
            "description": "Maximum water extent (seasonal + permanent water)",
            "threshold_months": f">= {seasonal_min_months}",
        }
    )

    min_water = xr.where(water_sum >= permanent_min_months, 1, 0).astype("uint8")
    min_water.name = "min_water"
    min_water.attrs.update(
        {
            "description": "Minimum water extent (permanent water only)",
            "threshold_months": f">= {permanent_min_months}",
        }
    )

    return min_water, max_water


def compute_water_seasonality(
    water_sum: xr.DataArray,
    valid_n: xr.DataArray,
) -> xr.DataArray:
    """
    Bin monthly water counts into five seasonality recurrence classes.

    Classes
    -------
    1 : 0–1 months 
    2 : 2–3 months
    3 : 4–6 months
    4 : 7–9 months
    5 : 10–12 months 
    255 : nodata (no valid observations)

    Parameters
    ----------
    water_sum : xr.DataArray (uint16)
        Count of wet months per pixel.
    valid_n : xr.DataArray (uint16)
        Count of valid (non-NaN) months per pixel.

    Returns
    -------
    xr.DataArray (uint8)
        Seasonality class raster.
    """
    ws = water_sum
    water_binned = xr.full_like(ws, fill_value=SEASONALITY_NODATA, dtype="uint8")
    water_binned = water_binned.where(~((ws >= 0) & (ws <= 1)), 1)
    water_binned = water_binned.where(~((ws >= 2) & (ws <= 3)), 2)
    water_binned = water_binned.where(~((ws >= 4) & (ws <= 6)), 3)
    water_binned = water_binned.where(~((ws >= 7) & (ws <= 9)), 4)
    water_binned = water_binned.where(~(ws >= 10), 5)
    water_binned = water_binned.where(valid_n > 0, SEASONALITY_NODATA).astype("uint8")
    water_binned.attrs.update(
        {
            "classes": {
                1: "0-1 months",
                2: "2-3 months",
                3: "4-6 months",
                4: "7-9 months",
                5: "10-12 months",
                SEASONALITY_NODATA: "nodata",
            },
            "description": "Water recurrence classes based on monthly water counts",
            "nodata": SEASONALITY_NODATA,
        }
    )
    return water_binned


def export_monthly_geotiffs(
    da: xr.DataArray,
    out_dir,
    prefix: str = "water_mask",
    cog: bool = True,
    dtype: str = "uint8",
    compress: str = COMPRESS,
    blocksize: int = BLOCKSIZE,
    predictor: int = PREDICTOR,
    overwrite: bool = False,
) -> None:
    """
    Export each time slice of *da* to a GeoTIFF (or COG) file.

    Parameters
    ----------
    da : xr.DataArray
        Dims: (time, y, x) with rioxarray spatial metadata.
    out_dir : str or Path
        Output directory (created if it does not exist).
    prefix : str
        Filename prefix; output name is ``{prefix}_{YYYYMMDD}.tif``.
    cog : bool
        Write Cloud-Optimised GeoTIFF if True, plain GeoTIFF otherwise.
    dtype : str
        Output data type (e.g. ``"uint8"``).
    compress : str
        GDAL compression (default ``"ZSTD"``).
    blocksize : int
        Tile size in pixels.
    predictor : int
        GDAL predictor (2 for integer, 3 for float).
    overwrite : bool
        Skip existing files when False.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if da.rio.crs is None:
        raise ValueError("DataArray has no CRS.  Set it first via da.rio.write_crs(...).")

    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)

    creation_options = {
        "TILED": "YES",
        "COMPRESS": compress,
        "BLOCKXSIZE": str(blocksize),
        "BLOCKYSIZE": str(blocksize),
        "PREDICTOR": str(predictor),
    }

    for t in da["time"].values:
        t_str = str(t)[:10].replace("-", "")
        out_path = out_dir / f"{prefix}_{t_str}.tif"

        if out_path.exists() and not overwrite:
            continue

        slice_da = da.sel(time=t).astype(dtype)

        try:
            driver = "COG" if cog else "GTiff"
            slice_da.rio.to_raster(out_path, driver=driver, **creation_options)
        except Exception as exc:
            slice_da.rio.to_raster(out_path, driver="GTiff", **creation_options)
            print(f"  COG write failed for {out_path.name}, wrote GTiff instead. ({exc})")

    print(f"Monthly GeoTIFFs written to: {out_dir}")


def _write_raster(da: xr.DataArray, path, dtype=None, **extra):
    """Write a single DataArray to a GeoTIFF with standard creation options."""
    opts = {**_CREATION_OPTIONS, **extra}
    kw = {"driver": "GTiff", **opts}
    if dtype:
        da = da.astype(dtype)
    da.rio.to_raster(path, **kw)


def export_annual_products(
    water_bin: xr.DataArray,
    out_dir,
    crs_source: xr.DataArray = None,
    nodata: int = NODATA,
    seasonal_min_months: int = 2,
    seasonal_max_months: int = 9,
    permanent_min_months: int = 10,
) -> None:
    """
    Derive and export all annual water products from monthly binary masks.

    Outputs under *out_dir*:
    annual_swf/annual_avrg_swf.tif               Mean water frequency [0-1]
    annual_swf/annual_median_swf.tif             Median water frequency [0-1]
    annual_classification/annual_classification  0=non-water,1=seasonal,2=permanent
    annual_water_extent/min_water_extent.tif     Permanent pixels only
    annual_water_extent/max_water_extent.tif     Seasonal + permanent pixels
    annual_seasonality/water_seasonality_class   5 recurrence bins

    Parameters
    ----------
    water_bin : xr.DataArray
        Monthly binary masks (0/1, 255 encoded as nodata), dims (time,y,x).
    out_dir : str or Path
        Root output directory.
    crs_source : xr.DataArray or None
        Source for geospatial metadata; defaults to *water_bin* itself.
    nodata : int
        NoData value encoded in *water_bin* (default 255).
    seasonal_min_months : int
        Classification threshold.
    seasonal_max_months : int
        Classification threshold.
    permanent_min_months : int
        Classification threshold.
    """
    out_dir = Path(out_dir)
    crs_source = crs_source if crs_source is not None else water_bin

    # Convert nodata -> NaN
    water_nan = water_bin.where(water_bin != nodata).astype("float32")

    # Annual surface water frequency
    swf_dir = out_dir / "annual_swf"
    swf_dir.mkdir(parents=True, exist_ok=True)

    mean_swf = water_nan.mean(dim="time", skipna=True)
    mean_swf.name = "swf"
    median_swf = water_nan.median(dim="time", skipna=True)
    median_swf.name = "swf"

    _write_raster(mean_swf, swf_dir / "annual_avrg_swf.tif")
    _write_raster(median_swf, swf_dir / "annual_median_swf.tif")
    print(f"SWF rasters written to: {swf_dir}")

    # Annual classification
    water_sum, valid_n, _occ, water_class = compute_annual_classification(
        water_nan,
        time_dim="time",
        seasonal_min_months=seasonal_min_months,
        seasonal_max_months=seasonal_max_months,
        permanent_min_months=permanent_min_months,
    )
    water_class = set_geo_metadata(water_class, crs_source)

    cls_dir = out_dir / "annual_classification"
    cls_dir.mkdir(parents=True, exist_ok=True)
    _write_raster(water_class, cls_dir / "annual_classification.tif", dtype="uint8")
    print(f"Classification raster written to: {cls_dir}")

    # Min/max water extent
    extent_dir = out_dir / "annual_water_extent"
    extent_dir.mkdir(parents=True, exist_ok=True)

    min_water, max_water = compute_minmax_water_extent(
        water_sum,
        seasonal_min_months=seasonal_min_months,
        permanent_min_months=permanent_min_months,
    )
    min_water = set_geo_metadata(min_water, crs_source)
    max_water = set_geo_metadata(max_water, crs_source)

    _write_raster(min_water, extent_dir / "min_water_extent.tif", dtype="uint8")
    _write_raster(max_water, extent_dir / "max_water_extent.tif", dtype="uint8")
    print(f"Water extent rasters written to: {extent_dir}")

    # Seasonal water classification
    seasonality_dir = out_dir / "annual_seasonality"
    seasonality_dir.mkdir(parents=True, exist_ok=True)

    water_binned = compute_water_seasonality(water_sum, valid_n)
    water_binned = set_geo_metadata(water_binned, crs_source)

    _write_raster(water_binned, seasonality_dir / "water_seasonality_class.tif", dtype="uint8")
    print(f"Seasonality raster written to: {seasonality_dir}")


def run_water_indicators(
    input_dir,
    output_dir,
    pattern: str = "water_*.tif",
    smooth_method: str = "gaussian",
    sigma: float = 0.8,
    filter_size: int = 3,
    high_thr: float = HIGH_THR,
    low_thr: float = LOW_THR,
    max_dist_m: float = MAX_DIST_M,
    connectivity: int = CONNECTIVITY,
    cog: bool = True,
    overwrite: bool = False,
    export_smoothed: bool = False,
    seasonal_min_months: int = 2,
    seasonal_max_months: int = 9,
    permanent_min_months: int = 10,
) -> None:
    """
    Run the full water indicators pipeline end-to-end.

    Steps
    -----
    1. Load monthly water probability GeoTIFFs from *input_dir*.
    2. Apply spatial smoothing (*smooth_method*).
    3. Apply hysteresis thresholding → monthly binary water masks.
    4. Export monthly binary masks to ``{output_dir}/monthly/``.
    5. Derive and export annual products (SWF, classification, extent,
       seasonality) to subdirectories under *output_dir*.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing monthly water probability TIFFs.
    output_dir : str or Path
        Root directory for all outputs.
    pattern : str
        Glob used to find input TIFFs (default ``"water_*.tif"``).
    smooth_method : {"gaussian", "median", "uniform"}
        Spatial smoothing to apply before thresholding.
    sigma : float
        Gaussian sigma (used only when *smooth_method="gaussian"*).
    filter_size : int
        Kernel size for median/uniform filters.
    high_thr : float
        Hysteresis high threshold (water seeds).
    low_thr : float
        Hysteresis low threshold (growth limit).
    max_dist_m : float
        Maximum hysteresis growth radius in metres.
    connectivity : {1, 2}
        Pixel connectivity for hysteresis growth.
    cog : bool
        Write Cloud-Optimised GeoTIFFs when True.
    overwrite : bool
        Re-create existing files when True.
    export_smoothed : bool
        When True, export the smoothed probability maps (float32) to
        ``{output_dir}/smoothed/`` before thresholding.
    seasonal_min_months : int
        Lower bound for seasonal class (class 1).
    seasonal_max_months : int
        Upper bound for seasonal class (class 1).
    permanent_min_months : int
        Pixels wet ≥ this count → permanent (class 2).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Loading monthly water maps from: {input_dir}")
    water_ts = load_monthly_water_maps(input_dir, pattern=pattern)

    print(f"[2/5] Smoothing with method={smooth_method!r}")
    smoothed = smooth_timeseries(
        water_ts,
        method=smooth_method,
        sigma=sigma,
        filter_size=filter_size,
    )

    if export_smoothed:
        print(f"[2b] Exporting smoothed probability maps")
        smoothed_export = smoothed.rio.write_crs(water_ts.rio.crs) if water_ts.rio.crs else smoothed
        export_monthly_geotiffs(
            smoothed_export,
            out_dir=output_dir / "smoothed",
            prefix="water_prob_smoothed",
            cog=cog,
            dtype="float32",
            predictor=3,
            overwrite=overwrite,
        )

    print(
        f"[3/5] Hysteresis thresholding (high={high_thr}, low={low_thr}, "
        f"max_dist={max_dist_m} m)"
    )
    water_hyst = hysteresis_timeseries(
        smoothed,
        high_thr=high_thr,
        low_thr=low_thr,
        max_dist_m=max_dist_m,
        connectivity=connectivity,
    )

    print(f"[4/5] Exporting monthly binary masks")
    export_monthly_geotiffs(
        water_hyst,
        out_dir=output_dir / "monthly",
        prefix="water_mask",
        cog=cog,
        dtype="uint8",
        overwrite=overwrite,
    )

    print(f"[5/5] Computing and exporting annual products")
    export_annual_products(
        water_hyst,
        out_dir=output_dir,
        crs_source=water_hyst,
        seasonal_min_months=seasonal_min_months,
        seasonal_max_months=seasonal_max_months,
        permanent_min_months=permanent_min_months,
    )

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Water Detection Post-Processing Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input-dir",
        required=True,
        metavar="PATH",
        help="Directory containing monthly water probability GeoTIFFs.",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        metavar="PATH",
        help="Root directory for all outputs.",
    )
    p.add_argument(
        "--pattern",
        default="water_*.tif",
        metavar="GLOB",
        help="Glob pattern for input TIFF files.",
    )
    p.add_argument(
        "--smooth",
        default="gaussian",
        choices=["gaussian", "median", "uniform"],
        dest="smooth_method",
        help="Spatial smoothing method applied before thresholding.",
    )
    p.add_argument(
        "--sigma",
        type=float,
        default=0.8,
        help="Gaussian sigma (ignored unless --smooth=gaussian).",
    )
    p.add_argument(
        "--filter-size",
        type=int,
        default=3,
        metavar="N",
        help="Kernel size for median/uniform filters.",
    )
    p.add_argument(
        "--high-thr",
        type=float,
        default=HIGH_THR,
        metavar="FLOAT",
        help="Hysteresis high threshold (water seeds).",
    )
    p.add_argument(
        "--low-thr",
        type=float,
        default=LOW_THR,
        metavar="FLOAT",
        help="Hysteresis low threshold (growth region).",
    )
    p.add_argument(
        "--max-dist",
        type=float,
        default=MAX_DIST_M,
        metavar="METRES",
        dest="max_dist_m",
        help="Maximum hysteresis growth radius in metres.",
    )
    p.add_argument(
        "--connectivity",
        type=int,
        default=CONNECTIVITY,
        choices=[1, 2],
        help="Pixel connectivity: 1=4-connected, 2=8-connected.",
    )
    p.add_argument(
        "--no-cog",
        action="store_false",
        dest="cog",
        help="Write plain GeoTIFF instead of Cloud-Optimised GeoTIFF.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-create output files even if they already exist.",
    )
    p.add_argument(
        "--export-smoothed",
        action="store_true",
        dest="export_smoothed",
        help="Export smoothed probability maps (float32) to {output-dir}/smoothed/.",
    )
    p.add_argument(
        "--nonwater-max-months",
        type=int,
        default=1,
        metavar="N",
        help="Pixels wet ≤ N months → non-water (class 0).",
    )
    p.add_argument(
        "--seasonal-min-months",
        type=int,
        default=2,
        metavar="N",
        help="Lower bound (inclusive) for seasonal class.",
    )
    p.add_argument(
        "--seasonal-max-months",
        type=int,
        default=9,
        metavar="N",
        help="Upper bound (inclusive) for seasonal class.",
    )
    p.add_argument(
        "--permanent-min-months",
        type=int,
        default=10,
        metavar="N",
        help="Pixels wet ≥ N months → permanent (class 2).",
    )
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    run_water_indicators(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        pattern=args.pattern,
        smooth_method=args.smooth_method,
        sigma=args.sigma,
        filter_size=args.filter_size,
        high_thr=args.high_thr,
        low_thr=args.low_thr,
        max_dist_m=args.max_dist_m,
        connectivity=args.connectivity,
        cog=args.cog,
        overwrite=args.overwrite,
        export_smoothed=args.export_smoothed,
        seasonal_min_months=args.seasonal_min_months,
        seasonal_max_months=args.seasonal_max_months,
        permanent_min_months=args.permanent_min_months,
    )
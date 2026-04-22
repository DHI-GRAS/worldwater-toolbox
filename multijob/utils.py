import pandas as pd
from openeo.extra.job_management import create_job_db, get_job_db
import numpy as np
import geopandas as gpd
from shapely.geometry import box
import plotly.express as px
import copy
from shapely import wkt
from shapely.geometry import box
import geopandas as gpd
import numpy as np
import os
from dateutil.relativedelta import relativedelta

import sys
sys.path.append(r"C:\Users\akou\OneDrive - DHI\Documents\repos\worldwater-toolbox\worldwater-toolbox\world_water_toolbox")
from wwt import  _water_extent_for_month 


def split_extent(spatial_extent, grid_size=10000, clip=True):
    # Build a single-geometry AOI from lon/lat bounds.
    bbox = gpd.GeoDataFrame(
        geometry=[box(spatial_extent["west"], spatial_extent["south"], spatial_extent["east"], spatial_extent["north"])],
        crs=spatial_extent.get("crs", 4326),
    )

    # Derive the local UTM zone to split in projected meter units.
    lon, lat = (spatial_extent["west"] + spatial_extent["east"]) / 2, (spatial_extent["south"] + spatial_extent["north"]) / 2
    zone = int((lon + 180) // 6) + 1
    utm = 32600 + zone if lat >= 0 else 32700 + zone

    # Convert AOI to UTM and create a regular square grid.
    b = bbox.to_crs(utm)
    xmin, ymin, xmax, ymax = b.total_bounds

    xs = np.arange(xmin, xmax, grid_size)
    ys = np.arange(ymin, ymax, grid_size)

    g = gpd.GeoDataFrame(
        geometry=[box(x, y, x + grid_size, y + grid_size) for x in xs for y in ys],
        crs=utm,
    )

    if clip:
        g = gpd.overlay(g, b, how="intersection")

    g = g.to_crs(4326)
    extents = [ {"west": geom.bounds[0], "east": geom.bounds[2], "south": geom.bounds[1], "north": geom.bounds[3],"crs": 4326} for geom in g.geometry ]
    geoms = g.geometry.tolist()
    return extents, geoms



def start_job(row, connection, **kwargs):

    row = row.to_dict()
    
    water_cube, _ = _water_extent_for_month(
        connection=connection,
        spatial_extent=row["spatial_extent"],
        region=row["region"],
        month_start=row["month_start"],
        month_end=row["month_end"],
        cloud_cover=row["cloud_cover"],
        threshold=row["threshold"],
        use_sentinelhub=row["use_sentinelhub"],
        rgb_processing=row["rgb_processing"],
        only_s1=row["only_s1"],
    )
    
    # Standardize naming so jobs and outputs can be grouped by month.
    month_tag = f"{row['month_start'][:7].replace('-', '_')}_{row['month_end'][:7].replace('-', '_')}"
    region_tag = row['region'].replace(" ", "_")
    prefix = "water_onlyS1" if row['only_s1'] else "water"
    filename_prefix = f"{prefix}_{month_tag}"

    print(f"Processing month: {month_tag}")
    # Create the backend job with GeoTIFF output and caching enabled.
    job = water_cube.create_job(
        title=f"{prefix}_{month_tag}",
        out_format="GTiff",
        job_options={"node_caching": True},
        filename_prefix=filename_prefix,
    )
    
    return job

def create_job_df(extent, region, month_start, month_end, cloud_cover, 
                   threshold, only_s1, use_sentinelhub, rgb_processing, 
                   single_month=True,output_dir='./',grid_size=5000):
    # Split AOI into manageable processing tiles.
    extents, _ = split_extent(extent, grid_size=grid_size)
    dfs = []
    
    df = pd.DataFrame({'spatial_extent':extents})   
    df['region'] = region
    
    df['cloud_cover'] = cloud_cover
    df['threshold'] = threshold
    df['only_s1'] = only_s1
    df['use_sentinelhub'] = use_sentinelhub
    df['rgb_processing'] = rgb_processing
    df['output_dir'] = output_dir

    if single_month:
        # Derive month_end from month_start for one monthly batch.
        month_end = (pd.to_datetime(month_start) + relativedelta(months=1)).strftime("%Y-%m-%d")
        df['month_start'] = month_start
        df['month_end'] = month_end
    else:
        # Build one dataframe copy per monthly interval in the requested range.
        current = pd.to_datetime(month_start)
        end = pd.to_datetime(month_end)
        while current < end:
            month_start = current.strftime("%Y-%m-%d")
            month_end = (current + relativedelta(months=1)).strftime("%Y-%m-%d")
            current += relativedelta(months=1)
            df['month_start'] = month_start
            df['month_end'] = month_end
            dfs.append(df.copy())
        df = pd.concat(dfs, ignore_index=True)

    return df

def plot_job_status(status_df, color_dict):
    import ast

    status_plot = copy.deepcopy(status_df)

    status_plot["geometry"] = status_plot["spatial_extent"].apply(ast.literal_eval).apply(
                    lambda d: box(d["west"], d["south"], d["east"], d["north"])
                )

    status_plot = gpd.GeoDataFrame(status_plot, geometry='geometry', crs='EPSG:4326')
    status_plot['color'] = status_plot['status'].map(color_dict).fillna(color_dict[None])

    minx, miny, maxx, maxy = status_plot.total_bounds
    center_lat = (miny + maxy) / 2
    center_lon = (minx + maxx) / 2

    fig = px.choropleth_map(
        status_plot,
        geojson=status_plot.geometry.__geo_interface__,
        locations=status_plot.index,
        color='status',
        color_discrete_map=color_dict,
        map_style="carto-positron",
        center={"lat": center_lat, "lon": center_lon},
        zoom=8,
        title="Job Status Overview"
    )
    fig.update_geos(fitbounds="locations")
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})

    return fig
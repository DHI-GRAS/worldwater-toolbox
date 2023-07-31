import pytest
from datetime import datetime as dt
from datetime import date as date
import openeo
from world_water_toolbox.wwt import s2_water_processing

spatial_extent = {'west': -100.08970034199719, 'east': -99.88035465800282, 'south': 19.39343565800281, 'north': 19.53214134199719,
     'crs': 4326}

connection = openeo.connect("openeo.vito.be").authenticate_oidc()

def test_hillshade():
    from world_water_toolbox import wwt

    start_date_exclusion = date(2020,12,1)
    month_end = date(2021,1,16)
    cube = wwt.hillshade_mask(connection,spatial_extent,start_date_exclusion,month_end,75)
    cube.download("hillshade.nc")

def test_masked_s2_cube():
    from world_water_toolbox import wwt

    start_date_exclusion = date(2020,12,1)
    month_end = date(2021,2,1)
    cube = wwt.masked_s2_cube(connection,spatial_extent,start_date_exclusion,month_end,75)
    cube.download("masked_clouds_hills3.nc")

def test_full_example():
    from world_water_toolbox import wwt

    start = '20210101'
    start = dt.strptime(start, "%Y%m%d").date()
    end = '20210401'
    end = dt.strptime(end, "%Y%m%d").date()
    region = 'Deserts'
    geometry = 'example/input_large/aoi.geojson'
    rgb_processing=False

    wwt.main("openeo.vito.be", start, end, region, geometry, rgb_processing, 85, 75, use_sentinelhub=True)

def test_full_example_cdse():
    from world_water_toolbox import wwt

    start = '20220101'
    start = dt.strptime(start, "%Y%m%d").date()
    end = '20220131'
    end = dt.strptime(end, "%Y%m%d").date()
    region = 'Deserts'
    geometry = 'example/input_large/aoi.geojson'
    rgb_processing=False

    wwt.main("openeo-staging.dataspace.copernicus.eu", start, end, region, geometry, rgb_processing, 85, 75, use_sentinelhub=False)

def test_multiple_months():
    from world_water_toolbox import wwt

    start = '20220101'
    start = dt.strptime(start, "%Y%m%d").date()
    end = '20220331'
    end = dt.strptime(end, "%Y%m%d").date()
    region = 'Deserts'
    geometry = 'example/input_large/aoi.geojson'
    rgb_processing=False

    connection = openeo.connect("openeo-staging.dataspace.copernicus.eu").authenticate_oidc()
    water_prob = wwt._water_extent_multiple_months(connection, start, end, geometry, region, 75, use_sentinelhub=False)
    water_prob.execute_batch("WWT probability", format="GTiff",job_options={
            "executor-memory": "2g",
            "executor-memoryOverhead": "5g",
            "executor-cores": 1
        }, filename_prefix="Worldwater_probability")





def test_local_processing():
    from openeo.local import LocalConnection

    local_conn = LocalConnection("./")
    collection = local_conn.load_collection("/home/driesj/python/worldwater-toolbox/tests/masked_clouds_hills3.nc")
    collection.metadata._orig_metadata["id"] = "SENTINEL2_L2A"
    s2_cube, ndxi_cube, s2_cube_water = s2_water_processing(collection,"Deserts")
    s2_cube_water.median_time().execute()
import pytest
from datetime import datetime as dt
from datetime import date as date
import openeo

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

    wwt.main("openeo.vito.be", start, end, region, geometry, rgb_processing, 85, 75)
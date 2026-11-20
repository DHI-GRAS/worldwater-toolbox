import os
from datetime import date
from dateutil.relativedelta import relativedelta
import openeo

from world_water_toolbox.wwt import _get_spatial_extent, _water_extent_for_month


def process_single_month(
    connection: openeo.Connection,
    month_start: str,
    month_end: str,
    region: str,
    threshold: int,
    cloud_cover: int,
    rgb_processing: bool,
    use_sentinelhub: bool,
    only_s1: bool,
    output_dir: str,
    geometry: str = None,
    spatial_extent: dict = None,
):
    if spatial_extent:
        spatial_extent = spatial_extent
    else:
        spatial_extent = _get_spatial_extent(geometry)

    water_cube, s2_cube = _water_extent_for_month(
        connection=connection,
        spatial_extent=spatial_extent,
        region=region,
        month_start=month_start,
        month_end=month_end,
        cloud_cover=cloud_cover,
        threshold=threshold,
        use_sentinelhub=use_sentinelhub,
        rgb_processing=rgb_processing,
        only_s1=only_s1,
    )

    month_tag = f"{month_start[:7].replace('-', '_')}_{month_end[:7].replace('-', '_')}"
    region_tag = region.replace(" ", "_")
    out_folder = os.path.join(output_dir, f"WWT_{region_tag}_{month_tag}")
    os.makedirs(out_folder, exist_ok=True)
    prefix = "water_onlyS1" if only_s1 else "water"
    filename_prefix = f"{prefix}_{month_tag}"

    print(f"Processing month: {month_tag}")
    job = water_cube.create_job(
        title=f"{prefix}_{month_tag}",
        out_format="GTiff",
        job_options={"node_caching": True},
        filename_prefix=filename_prefix,
    )
    results = job.start_and_wait().get_results()
    results.download_files(out_folder)

    if rgb_processing:
        s2_month = s2_cube.filter_temporal([month_start, month_end]).median_time()
        job = s2_month.create_job(
            title=f"median_{month_tag}",
            out_format="GTiff",
            job_options={"node_caching": True},
            filename_prefix=f"median_{month_tag}",
        )
        results = job.start_and_wait().get_results()
        results.download_files(out_folder)

    return out_folder


def process_annual_water(
    connection: openeo.Connection,
    year: int,
    region: str,
    threshold: int,
    cloud_cover: int,
    rgb_processing: bool,
    use_sentinelhub: bool,
    only_s1: bool,
    output_dir: str,
    geometry: str = None,
    spatial_extent: dict = None,
):
    current = date(year, 1, 1)
    end = date(year, 12, 31)

    while current <= end:
        month_start = current.strftime("%Y-%m-%d")
        month_end = (current + relativedelta(months=1)).strftime("%Y-%m-%d")

        folder = process_single_month(
            connection=connection,
            month_start=month_start,
            month_end=month_end,
            region=region,
            threshold=threshold,
            cloud_cover=cloud_cover,
            rgb_processing=rgb_processing,
            use_sentinelhub=use_sentinelhub,
            only_s1=only_s1,
            output_dir=output_dir,
            geometry=geometry,
            spatial_extent=spatial_extent,
        )

        print(f"Saved monthly result: {folder}")

        current = current + relativedelta(months=1)


def test_udp(connection: openeo.Connection):
    spatial_extent = {
        "west": -9.328873376713767,
        "east": -9.205921388808221,
        "south": 53.623510495414074,
        "north": 53.71533891769306,
        "crs": 4326,
    }
    cube = connection.datacube_from_process(
        process_id="worldwater_water_extent",
        start_date="2024-01-01",
        only_s1=False,
        use_sentinelhub=True,
        rgb_processing=False,
        cloud_cover=80,
        water_threshold=75,
        region="Temperate grassland",
        spatial_extent=spatial_extent,
        namespace="https://raw.githubusercontent.com/DHI-GRAS/worldwater-toolbox/refs/heads/new-udp/world_water_toolbox/single_month_extent_udp.json",
    )
    cube.execute_batch("WWT_Temperate_grassland_2024_01_2024_02.tif")

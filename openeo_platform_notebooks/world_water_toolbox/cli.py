import click
from datetime import datetime as dt


@click.group()
def cli():
    """Command on cli1"""
    pass


class DateType(click.ParamType):
    name = 'date'

    def convert(self, date, _, __):
        try:
            return dt.strptime(date, "%Y%m%d").date()
        except ValueError as exception:
            msg = "Not a valid date: '{0}'.".format(date)
            raise ValueError(msg) from exception
        
        
DATEIN = DateType()

@click.command()
@click.option('-s', 
              '--start',
              type = DATEIN,
              required =True,
              help = 'Start Date of Processing (yyyymmdd), e.g.: 20200101')


@click.option('-e',
              '--end', 
              type = DATEIN, 
              required = True,
              help = 'End Date of Processing (yyyymmdd, e.g.: 20201231')



@click.option(
    "-g",
    "--geometry",
    type = click.Path(dir_okay=False, file_okay=True),
    required = True,
    help = "Path to the aoi (geojson), e.g.: C:/wwt/aoi.geojson")


@click.option('-c',
              '--cloud_cover',
              type = click.FloatRange(0, 100),
              default = 95,
              show_default = True,
              required = True,
              help = 'Cloud Cover Percentage')


@click.option('-t',
              '--threshold',
              type = click.FloatRange(0, 100),
              default = 75,
              show_default = True,
              required = True,
              help = 'Threshold')

@click.option('-rgb',
              '--rgb_processing',
              is_flag = False,
              required = True,
              help = 'If you wish to process an optical(RGB) monthly median image, then it is True and "False" if you do not want to process an RGB image')


@click.option(
    "-r",
    "--region",
    type = click.Choice(['Deserts', 'Mountain','Tropical forest','Tropical savanna','Subtropical forest',
                         'Subtropical savanna','Temperate broadleaf','Temperate grassland']),
    required = True,
    help = "Define ecoregion from options")



def wwt(start, end, region, geometry, rgb_processing, cloud_cover, threshold, **kwargs):

    """Calculate water extent for each moth between start and end period for the AOI
       and generate rgb image for each month (- rgb True).
    
    Example use:

    world_water_toolbox -s 20200101 -e 20201231 -g A:/ANSU/6_Tasks/2204_WorldWater_TBX/python_package/test/aoi.geojson -rgb True -r Deserts 

    """
    
    from world_water_toolbox import wwt
    wwt.main(start, end, region, geometry, rgb_processing, cloud_cover, threshold)


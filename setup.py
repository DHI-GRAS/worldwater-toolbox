from setuptools import setup, find_packages

setup(
    name='world_water_toolbox',
    version='0.1.0',
    description='Calculate water extent for each moth between start and end period for the AOI',
    author='Andrea Sulova',
    author_email='ansu@dhigroup.com',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'openeo',
        'Click',
        'Numpy',
        'rioxarray',
        'geopandas',
        'rasterio',
        'xarray',
        
    ],
    entry_points='''
        [console_scripts]
        world_water_toolbox = world_water_toolbox.cli:wwt
    '''
)
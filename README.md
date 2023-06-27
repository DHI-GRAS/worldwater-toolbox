##  World Water Toolbox [openEO]
### An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics 


We present a simple yet robust method for fusing optical and Synthetic Aperture Radar (SAR) data for mapping surface water dynamics over.

To get started, we have included all the necessary information here, including sign in to the openEO platfom,  installation instructions and a list of fundamental commands. To provide you with a better understanding of the package's capabilities, we are in the process of preparing a set of examples that demonstrate various use cases. 

### Installing python package:

```python
pip install git+https://github.com/sulova/world_water_toolbox.git

```

After installation, the **world water python package** will contain the `world_water_toolbox` command-line interface.
```
world_water_toolbox --help
```
### Individual Commands

The individual command available within world_water_toolbox is e.g.:
```
world_water_toolbox -b openeo-dev.vito.be -s 20210101 -e 20210201 -r Deserts -g A:\ANSU\6_Tasks\2204_WorldWater_TBX\world_water_toolbox\Test\aoi.geojson -rgb True
```

### Jyputer Notebook 
Using Notebooks on openeo platform!!!


If you find this repository useful in your research, please consider citing the following papers. Thank you.

- Druce, D.; Tong, X.; Lei, X.; Guo, T.; Kittel, C.M.M.; Grogan, K.; Tottrup, C. An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics over Mainland China. Remote Sens. 2021, 13, 1663. https://doi.org/10.3390/rs13091663

Should you have any questions, require assistance, or wish to provide feedback, please don't hesitate to reach out to us.
 
Copyright (C) 2023 DHI A/S


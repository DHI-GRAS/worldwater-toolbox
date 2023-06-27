##  World Water Toolbox [openEO]
### An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics 


We present a simple yet robust method for fusing optical and Synthetic Aperture Radar (SAR) data for mapping surface water dynamics over. The full mapping workflow is implemented in openEO cloud-computing platform that provides the necessary combination of large-scale computing resources with instant access to planetary scale archives of satellite imagery.

To get started, we have included all the necessary information here, including sign in to the openEO platfom,  installation instructions and a list of fundamental commands. To provide you with a better understanding of the package's capabilities, we are in the process of preparing a set of examples that demonstrate various use cases. 


### 1) Register for openEO Platform 
Register for openEO Platform and immediately start a 30 day free trial. Free trial users receive 1000 free credits upon registration. openEO Platform is a federation of services which include:
- EGI for Authentification & Authorization via EGI Check-in(opens new window)
- Terrascope EOPlaza for Account management
- EODC for provisioning of JupyterLab
- Terrascope including Sentinel Hub connection via the EuroDataCube, EODC and Sentinel Hub backends

![image](https://github.com/DHI-GRAS/worldwater-toolbox/assets/44543964/1d83ede8-886e-4592-b98a-b4bbb9ccbab5)

### 2) Installing python package

```python
pip install git+https://github.com/sulova/world_water_toolbox.git

```

After installation, the **world water python package** will contain the `world_water_toolbox` command-line interface.
```
world_water_toolbox --help
```
### 3) Individual Commands

The individual command available:
```
world_water_toolbox -b openeo-dev.vito.be -s 20210101 -e 20210201 -r Deserts -g A:\ANSU\6_Tasks\2204_WorldWater_TBX\world_water_toolbox\Test\aoi.geojson -rgb True
```

### Jyputer Notebook 
We provide a user friendly Jyputer Notebook interface to collect the user input and see the impact the changes have on the data/results, without having to interact with the code;
![image](https://github.com/DHI-GRAS/worldwater-toolbox/assets/44543964/cbd1c03e-ca0f-4218-add8-d219da1e3ca0)


Should you have any questions, require assistance, or wish to provide feedback, please don't hesitate to reach out to us. 
If you find this repository useful in your research, please consider citing the following papers. Thank you.

- Druce, D.; Tong, X.; Lei, X.; Guo, T.; Kittel, C.M.M.; Grogan, K.; Tottrup, C. An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics over Mainland China. Remote Sens. 2021, 13, 1663. https://doi.org/10.3390/rs13091663

 
Copyright (C) 2023 DHI A/S


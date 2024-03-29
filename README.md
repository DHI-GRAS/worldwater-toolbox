##  World Water Toolbox [openEO]

![Picture2](https://github.com/DHI-GRAS/worldwater-toolbox/assets/44543964/a0d7e432-d034-4f74-aee1-6d1b603cda1a)


### An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics 

We present a simple yet robust method for fusing optical and Synthetic Aperture Radar (SAR) data for mapping surface water dynamics. The full mapping workflow is implemented in the openEO cloud-computing platform that provides the necessary combination of large-scale computing resources with instant access to planetary scale archives of satellite imagery.

To get started, you need to follow these steps: 
1)  Sign in to the openEO platform. If you do not have an account, you may need to create one by following the registration process https://docs.openeo.cloud/join/free_trial.html.

2) Choose between two options to run the World Water Toolbox:
- Option 1: Running Jupiter Notebook on the openEO platform https://github.com/DHI-GRAS/worldwater-toolbox/tree/main/Jupiter%20Notebook
- Option 2: Installing the Python package locally and leveraging an openEO backend for processing. If you prefer to work with the World Water Toolbox on your local machine, you can install the Python package and utilize one of the available openEO backends for processing. Follow the steps below to get started:

3) You can check results in your directory or on the Open Editor:  https://editor.openeo.org 
   
### 1) Register for openEO Platform: 
![logo](https://github.com/DHI-GRAS/worldwater-toolbox/assets/44543964/42768ef0-a2f8-46c5-aaf9-6099d4783f60)

Register for openEO Platform and immediately start a 30 day free trial. Register: https://docs.openeo.cloud/join/free_trial.html. Free trial users receive 1000 free credits upon registration. 
The openEO Platform is a federation of services which include:
- EGI for Authentification & Authorization via EGI Check-in
- Terrascope EOPlaza for Account management
- EODC for provisioning of JupyterLab
- Terrascope including Sentinel Hub connection via the EuroDataCube, EODC and Sentinel Hub backends



###  2) Option 1: Jupiter Notebook [eoLab]
You can run the Jupiter notebook on the openEO platform: https://eolab.eodc.eu/hub/login?next=%2Fhub%2F
We provide a user friendly Jupiter  Notebook interface to collect the user input without having to interact with the code.
![image](https://github.com/DHI-GRAS/worldwater-toolbox/assets/44543964/cbd1c03e-ca0f-4218-add8-d219da1e3ca0)

### 2) Option 2: Installing python package

```python
pip install git+https://github.com/sulova/world_water_toolbox.git

```

After installation, the **world water python package** will contain the `world_water_toolbox` command-line interface.
```
world_water_toolbox --help
```
The individual command available:
```
world_water_toolbox -b openeo-dev.vito.be -s 20210101 -e 20210201 -r Deserts -g A:\ANSU\6_Tasks\2204_WorldWater_TBX\world_water_toolbox\Test\aoi.geojson -rgb True
```
### 3) Check results 
You can check results in your directory or on the Open Editor: https://editor.openeo.org
Should you have any questions, require assistance, or wish to provide feedback, please don't hesitate to reach out to us. 

More information will be soon found on the openEO Cookbook {insert link}

If you find this repository useful in your research, please consider citing the following papers. Thank you.
- Druce, D.; Tong, X.; Lei, X.; Guo, T.; Kittel, C.M.M.; Grogan, K.; Tottrup, C. An Optical and SAR Based Fusion Approach for Mapping Surface Water Dynamics over Mainland China. Remote Sens. 2021, 13, 1663. https://doi.org/10.3390/rs13091663

 
Copyright (C) 2023 DHI A/S


# -*- coding: utf-8 -*-
"""
Created on Tue May 10 20:37:13 2022

@author: dadr
"""

import water_mask
import rioxarray as rxr
import xarray as xr
import numpy as np
# import matplotlib

## select an ecoregion from which logistic regression parameters will be extracted.
## 'tropical' or 'subtropical'
ecoregion = 'tropical'

params = water_mask.LOOKUPTABLE[ecoregion] ##shouldn't need to be called here

### Test data Sentinel-1 and Sentinel-2.
filepaths = {
    'S1': r'C:\Users\dadr\OneDrive - DHI\Desktop\projects\world-water\s1-mean.tif',
    'S2': r'C:\Users\dadr\OneDrive - DHI\Desktop\projects\world-water\s2-median.tif'
    }

img_dict = {}

for i in ['S2', 'S1']:
    
    # Get the image filepath
    f = filepaths[i]
    
    # Get the required bands for that image.
    bands = params[i]['bands']
    
    # Open the image with rioxarray
    with rxr.open_rasterio(f, chunks=2048, masked=True) as src:
        
        ### Check the required bands exist.
        if all(b in src.long_name for b in bands): 

            # create a dictionary of band:array for required inputs.
            d = {b: src[src.long_name.index(b)].rename(b) for b in bands}  
        
            # Add the image mask. Required to composite the data in the correct order.
            d[f'{i}_mask'] = np.isnan(d[list(d.keys())[0]]).rename(f'{i}_mask')
        
            # Update the dictionary.
            img_dict.update(d)

# create an empty xarray dataset.
ds = xr.Dataset()

# Add the band arrays to the dataset
for key in img_dict.keys():
    ds[key] = img_dict[key]

# Pass the xarray dataset to the water_mask function, and specify the ecoregion.
test = water_mask.water_mask(ds, 'tropical')


# optionally plot the data

# # Probability is a uint8 array (0-100) which can be threshold to define a water mask.
# test['probability'].plot.imshow()

# # sensor_rank tells us which sensor that probability came from. 
# # S2_S1: 5, S2: 4, S1: 3
# test['sensor_rank'].plot.imshow()



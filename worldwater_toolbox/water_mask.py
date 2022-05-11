# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 09:24:52 2020

@author: dadr
"""

import numpy as np
import xarray as xr


## Note this will be updated. Currently demonstration only. 
## The lookup table defines bands and logistic regression expressions 
## for each ecoregion.
LOOKUPTABLE = {
    'tropical': 
    {
        'S1': {'bands': ['VV'], 
               'expression': lambda VV: 1 / (1 + np.exp(- (-7.17 + (-0.48 * VV))))},
        'S2': {'bands': ['NDVI', 'NDWI'], 
               'expression': lambda NDVI, NDWI: 1 / (1 + np.exp(- (0.845 + (2.14 * NDVI) + (-13.5 * NDWI))))},
        'S2_S1': {'bands': ['VV', 'NDWI'], 
                  'expression': lambda VV, NDWI: 1 / (1 + np.exp(- (-2.64 + (-0.23 * VV) + (-8.6 * NDWI))))},
     },
    
    'subtropical':
    {
        'S1': {'bands': ['VV', 'VH'], 
               'expression': lambda VV, VH: 1 / (1 + np.exp(- (-8.1 + (-0.13 * VV) + (-0.27 * VH))))},
        'S2': {'bands': ['NDVI', 'NDWI'], 
               'expression': lambda NDVI, NDWI: 1 / (1 + np.exp(- (1.267 + (0.316 * NDVI) + (-11 * NDWI))))},
        'S2_S1': {'bands': ['VV', 'NDWI'], 
                  'expression': lambda VV, NDWI: 1 / (1 + np.exp(- (-2.64 + (-0.23 * VV) + (-8.6 * NDWI))))},
    }
}
    
    
def _create_img(expr_dict, dataset, sensor_value):
    """
    
    Generates water probability from a logistic regression expression.
    
    Parameters
    ----------
    expr_dict : DICTIONARY
    
        Contains keys 'expression' and 'bands'.
        - Bands value is a list of inputs to the logistic regresion expression.
        - Expression value is a lambda function that takes bands as input. 
        
    dataset : XARRAY DATASET
    
        An xarray dataset containing the required bands in expr_dict.    
    
    sensor_value : UINT8
    
        An integer value that represents the sensor inputs. Used to produce
        a useful QA output. 

    Returns
    -------
    XARRAY DATASET
    
        (1) probability - surface water probability (0-100)
        (2) sensor_rank - A QA value that represents the input imagery.

    """
    
    # lambda expression to create surface water probability
    exp = expr_dict['expression']
    
    # bands required as input to the lambda expression
    bands = expr_dict['bands']
    
    # filter the dictionary to required bands identified above. 
    filtered_dict = dict((k, dataset[k]) for k in bands)
    
    # Calculate the probability using the input dataset and lambda expression.
    probability = (exp(**filtered_dict)*100).astype(np.uint8)
    
    # Create the QA band with sensor value.
    sensor_rank = xr.full_like(probability, sensor_value)
    
    # create an empty dataset and add the probability and sensor_rank outputs.
    ds = xr.Dataset()
    ds['probability'] = probability
    ds['sensor_rank'] = sensor_rank

    return ds    

def _water_mask(dataset, region_parameters):
    """
    
    Composites water probabilities according to data availability and accuracy.
    
    Parameters
    ----------
    dataset : XARARY DATASET
    
        An xarray dataset containing the required bands in region_parameters.    
        
    region_parameters : DICTIONARY
        
        A nested dictionary with first keys as the sensor (S1, S2...).
        
        Followed by key value pairs of 'expression' and 'bands'.
        - 'bands' : Bands value is a list of inputs to the logistic regresion expression.
        - 'expression' : value is a lambda function that takes bands as input. 


    Returns
    -------
    composite : XARRAY DATASET
    
        (1) probability - surface water probability (0-100)
        (2) sensor_rank - A QA value that represents the input imagery.
            S2_S1 : 5, S2: 4, S1: 3...

    """
    
    composite = xr.where(
        
        # Where valid data for both S2 and S1.
        xr.ufuncs.logical_and(~dataset['S2_mask'], ~dataset['S1_mask']), 
        
        # Create water mask with sensor_value 5
        _create_img(region_parameters['S2_S1'], dataset, 5),
        
        # Else
        xr.where(
            
            # Where Sentinel-2 data is valid.
            ~dataset['S2_mask'], 
            
            # Create water mask from sentinel-2 with sensor value 4.
            _create_img(region_parameters['S2'], dataset, 4),
            
            # If no Sentinel-2 data, use Sentinel-1 on it's own. Sensor value 3.
            _create_img(region_parameters['S1'], dataset, 3), 
            )    
        )    
    
    return composite

def water_mask(dataset, ecoregion):
    """

    Parameters
    ----------
    dataset : XARRAY DATASET
    
        An xarray dataset containing the required bands in ecoregion.    
        
    ecoregion : STRING ['tropical', 'subtropical']
        
        Your region of interest.
        
        Used to access a LOOKUPTABLE that has logistic regression expressions
        and required bands encoded for each sensor (S2, S1, S2_S1) that is 
        parameterized for that ecoregion.

    Returns
    -------
    composite : XARRAY DATASET
    
        (1) probability - surface water probability (0-100)
        (2) sensor_rank - A QA value that represents the input imagery.
            S2_S1 : 5, S2: 4, S1: 3...
    """
    
    # Access the parameters required for that ecoregion.
    region_parameters = LOOKUPTABLE[ecoregion]
    
    # Create the water mask.
    composite = _water_mask(dataset, region_parameters)

    return composite
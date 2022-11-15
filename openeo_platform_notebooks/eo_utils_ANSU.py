# Import necessary libraries
from time import time
import rasterio
import numpy as np
# import openeo
# from openeo.rest.datacube import PGNode, THIS
# from openeo.processes import *
import math
import xarray as xr
import rioxarray
from ipyleaflet import (
    Map,
    Marker,
    TileLayer, ImageOverlay,
    Polygon, Rectangle,
    GeoJSON,
    DrawControl,
    LayersControl,
    WidgetControl,
    basemaps,
    basemap_to_tiles,
    WMSLayer,
    FullScreenControl
)
import pip  as widgets
from ipywidgets import Output, FloatSlider
from traitlets import link
import shapely.geometry
import os
import xarray_leaflet
from xarray_leaflet.transform import passthrough
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from bqplot import Lines, Figure, LinearScale, DateScale, Axis, Scatter
import bqplot.pyplot as bqplt
from datetime import datetime
import json
import hvplot.xarray
import PIL
from base64 import b64encode
from io import StringIO, BytesIO
import warnings
warnings.filterwarnings("ignore")

class openMap:
    def __init__(self,center,zoom):
        self.map = Map(center=center, zoom=zoom, scroll_wheel_zoom=True, interpolation='nearest')
        self.bbox = []
        self.point_coords = []
        self.figure = None
        self.polygon = []
        self.figure_widget = None
        feature_collection = {
            'type': 'FeatureCollection',
            'features': []
        }

        draw = DrawControl(
            circlemarker={}, polyline={}, 
            polygon={"shapeOptions": {
                       "original": {},
                       "editing": {},
            }},           
            marker= {"shapeOptions": {
                       "original": {},
                       "editing": {},
            }},
            rectangle = {"shapeOptions": {
                       "original": {},
                       "editing": {},
            }})

        self.map.add_control(draw)
        def handle_draw(target, action, geo_json):
            feature_collection['features'] = []
            feature_collection['features'].append(geo_json)
            if feature_collection['features'][0]['geometry']['type'] == 'Point':
                self.point_coords = feature_collection['features'][0]['geometry']['coordinates']
            if feature_collection['features'][0]['geometry']['type'] == 'Polygon':
               self.polygon = feature_collection['features'][0]['geometry']['coordinates']
            else:
                coords = feature_collection['features'][0]['geometry']['coordinates'][0]
                polygon = shapely.geometry.Polygon(coords)
                self.bbox = polygon.bounds

        
        layers_control = LayersControl(position='topright')
        self.map.add_control(layers_control)
        self.map.add_control(FullScreenControl())
        self.map.add_layer(basemap_to_tiles(basemaps.Esri.WorldImagery));
        draw.on_draw(handle_draw)
    
    def getBbox(self):
        if(len(self.bbox) == 0):
            mapBox = self.map.bounds     
            return [ mapBox[0][1],mapBox[0][0],mapBox[1][1],mapBox[1][0]]
        else:
            return self.bbox
            
    def getPoint(self):
        if(len(self.point_coords) != 0):
            points = self.point_coords
            return points
    
    def getPolygons(self):
        if(len(self.polygon) != 0):
            polygon = self.polygon
            return polygon

      

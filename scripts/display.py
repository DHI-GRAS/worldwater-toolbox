import folium
from folium import plugins, GeoJson
import geopandas as gpd
import os
from folium.features import GeoJsonTooltip
import ipywidgets as widgets
from datetime import date
import traitlets
import sys
from dateutil.relativedelta import *
sys.path.append("..")
from world_water_toolbox import wwt
import warnings
import glob
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import matplotlib.pyplot as plt
from datetime import date, datetime
import rioxarray
import xarray as xr
import numpy as np

# Function to create and return a Folium map with the given GeoDataFrame
class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


def define_widgets():
    # Define all widgets
    zone_w = widgets.RadioButtons(
        options=['Deserts', 'Mountain','Tropical forest','Tropical savanna','Subtropical forest',
                 'Subtropical savanna','Temperate broadleaf','Temperate grassland', 'Tundra'],
        layout={'width': 'max-content'},
        description='Ecoregions',
        disabled=False)

    start_date_w = widgets.DatePicker(
        description='Start Date',
        value=date(2023, 5, 1),
        disabled=False)

    end_date_w = widgets.DatePicker(
        description='End Date',
        value=date(2023, 5, 31),
        disabled=False)

    threshold = widgets.IntSlider(value=75, description='Threshold',)
    threshold_cloud_cover = widgets.IntSlider(value=50, description='Cloud Cover',)

    RGB_processing = widgets.RadioButtons(
        value='No', 
        options=['Yes', 'No'], 
        description='Export RGB image',
        style={'description_width': 'initial', 'width': '400px'})

    s1_box = widgets.Checkbox(value=True, description='Only Sentinel-1', indent=False)
    textbox = widgets.Text(
        value='',
        placeholder='Type a directory',
        description='Save folder:',
        disabled=False)

    button = widgets.Button(description="Run")
    output = widgets.Output()

    get_data_button = LoadedButton(description='Run',
                                   disabled=False,
                                   icon='check',
                                   value='')

    return {
        "zone_w": zone_w,
        "start_date_w": start_date_w,
        "end_date_w": end_date_w,
        "threshold": threshold,
        "threshold_cloud_cover": threshold_cloud_cover,
        "RGB_processing": RGB_processing,
        "s1_box": s1_box,
        "textbox": textbox,
        "button": button,
        "output": output,
        "get_data_button": get_data_button
    }

def initialize_map():

    
    try:
        warnings.filterwarnings("ignore")
        gdf = gpd.read_file('aoi.geojson')
        center = gdf.centroid
        latitude = center.y[0]
        longitude = center.x[0]
        zoom_start= 10
        
    except:
        latitude = 0
        longitude = 0
        zoom_start= 2 


    # Initialize a Folium map centered on the calculated mean latitude and longitude
    m = folium.Map(location= [latitude,longitude], tiles= None, zoom_start=zoom_start).add_to(folium.Figure(height = 800))
    tile_layer = folium.TileLayer( tiles = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", 
                           attr = "Tiles © Esri — Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community",
                           name = 'Satellite').add_to(m)
    
    if os.path.exists('aoi.geojson'):
        geojson_group = folium.GeoJson(gpd.read_file('aoi.geojson'), name="AoI").add_to(m)

    draw = plugins.Draw(export=True,  filename='aoi.geojson', position='topleft').add_to(m)

    return m


def handle_button_click(widgets, connection, m):
    print("Button clicked")  # Debug line
    output = widgets['output']
    
    try:
        output = widgets['output']
        spatial_extent = wwt._get_spatial_extent('aoi.geojson')

        with output:

            months_list = []
            final_end_date = widgets['end_date_w'].value
            current_month = widgets['start_date_w'].value

            while current_month < final_end_date:
                months_list.append(current_month)
                current_month = current_month + relativedelta(months = 1)

            start_date_w = widgets['start_date_w'].value
            end_date_w = widgets['end_date_w'].value
            zone_w = widgets['zone_w'].value
            threshold = widgets['threshold'].value
            threshold_cloud_cover = widgets['threshold_cloud_cover'].value
            RGB_processing = True if widgets['RGB_processing'].value == 'Yes' else False
            s1_box = widgets['s1_box'].value
            output_dir = widgets['textbox'].value

            print('Start and end dates:', start_date_w, end_date_w)
            print('Zone:', zone_w)
            print('Threshold:', threshold)
            print('Cloud Cover', threshold_cloud_cover)
            
            for month_start in months_list:
                month_end = month_start + relativedelta(months = 1) 
                output, s2_rgb =  wwt._water_extent_for_month(connection, spatial_extent, zone_w, month_start, 
                                                         month_end, threshold, threshold_cloud_cover, True, RGB_processing, s1_box)

                inp = 'onslyS1' if s1_box else 'S1_S2'
                output_dir += f'/{str(start_date_w.strftime("%Y_%m")) }_{str(end_date_w.strftime("%Y_%m"))}'
                outfile = output_dir + f'/water_{str(month_start.strftime("%Y_%m"))}_{inp}.tif'

                output.execute_batch(outputfile=outfile)

                if RGB_processing:
                    print('S2 median image processing between: ' + str(start_date_w.strftime("%Y_%m")) + ' and ' + str(end_date_w.strftime("%Y_%m")))
                    my_job = s2_rgb.create_job(title= 'median_' + str(start_date_w.strftime("%Y_%m")) + '_' + str(end_date_w.strftime("%Y_%m")), 
                                                       out_format="GTiff",job_options= {"node_caching":True},
                                                       **{'filename_prefix': 'median_' + str(start_date_w.strftime("%Y_%m")) + '_' + str(end_date_w.strftime("%Y_%m"))})
                    results = my_job.start_and_wait().get_results()
                    results.download_files(output_dir)
            # Convert tif files to wgs coord system
            print('Output folder: ', output_dir)
            print('You can check results on the Open Editor: https://editor.openeo.org/?server=https%3A%2F%2Fopeneo-dev.vito.be')

            for file_path in glob.glob(output_dir + '/*.tif'):
                if not 'wgs' in file_path:
                    output_wgs = os.path.dirname(file_path)+ '/' + file_path.split('/')[-1][:-4]+'_wgs.tif'
                    with rasterio.open(file_path) as src:
                        transform, width, height = calculate_default_transform(src.crs, 'EPSG:4326', src.width, src.height, *src.bounds)
                        kwargs = src.meta.copy()
                        kwargs.update({'crs': 'EPSG:4326','transform': transform,'width': width,'height': height})

                        with rasterio.open(output_wgs, 'w', **kwargs) as dst:
                            for i in range(1, src.count + 1):
                                reproject(source=rasterio.band(src, i),
                                            destination=rasterio.band(dst, i),
                                            src_transform=src.transform,
                                            src_crs=src.crs,
                                            dst_transform=transform,
                                            dst_crs='EPSG:4326',
                                            resampling=Resampling.nearest)


            # Adding S2 median wgs outputs to the map 
            if RGB_processing:
                for file_path in sorted(glob.glob(output_dir + '/median*wgs.tif')):
                    if os.path.isfile(file_path):
                        rgb_layer  = rioxarray.open_rasterio(file_path)
                        mlat = rgb_layer.y.values.min() 
                        mlon = rgb_layer.x.values.min()
                        xlat = rgb_layer.y.values.max()
                        xlon = rgb_layer.x.values.max()
                        gamma = 20

                        # Select the bands for RGB visualization
                        red_band = rgb_layer[2] * gamma # Select the red band
                        green_band = rgb_layer[1] * gamma  # Select the green band
                        blue_band = rgb_layer[0] * gamma # Select the blue band
                        rgb_image = np.stack((red_band, green_band, blue_band), axis=-1)
                        rgb_image = np.uint16(rgb_image)
                        m.add_child(folium.raster_layers.ImageOverlay(rgb_image, bounds=[[mlat, mlon], [xlat, xlon]], show = False,  name = 'S2 Image '+ file_path.split('/')[-1][7:14].replace("_", " ")))


            # Adding water wgs outputs to the map
            for file_path in sorted(glob.glob(output_dir + '/water*wgs.tif')):
                filename = file_path.split('/')[-1].split('.')[0]
                if os.path.isfile(file_path):
                    binary_layer  = rioxarray.open_rasterio(file_path).drop('band')[0]
                    mlat = binary_layer.y.values.min() 
                    mlon = binary_layer.x.values.min()
                    xlat = binary_layer.y.values.max()
                    xlon = binary_layer.x.values.max()

                    def colorize(array, cmap = 'Blues'):
                        normed_data = (array - array.min()) / (array.max() - array.min())  
                        cm = plt.cm.get_cmap(cmap)    
                        return cm(normed_data)  
                    if 'water' in file_path:
                        print(file_path)
                        m.add_child(folium.raster_layers.ImageOverlay(colorize(binary_layer), [[mlat, mlon], [xlat, xlon]], opacity=0.8, show = False,  name = filename))#.set_active(False)

            # wwt._water_indicators(output_dir, inp)

            # for file_path in sorted(glob.glob(output_dir + '/index*.tif')): 
            #     print(file_path)

            #     index_layer  = rioxarray.open_rasterio(file_path).drop('band')[0]
            #     mlat = binary_layer.y.values.min() 
            #     mlon = binary_layer.x.values.min()
            #     xlat = binary_layer.y.values.max()
            #     xlon = binary_layer.x.values.max()

            folium.LayerControl(collapsed=True).add_to(m)

            display(m)
    except Exception as e:
        print("Error occurred:", e)
        raise  # Optionally re-raise the exception if you want to stop execution

    display(output)
    
def generate_map(connection):
    widgets = define_widgets()
    m = initialize_map()

    display(widgets["start_date_w"])
    display(widgets["end_date_w"])
    display(widgets["zone_w"])
    display(widgets["threshold"])
    display(widgets["threshold_cloud_cover"])
    display(widgets["RGB_processing"])
    display(widgets["s1_box"])
    display(widgets["textbox"])
    display(m)
    display(widgets["button"], widgets['output'])
    
    widgets["button"].on_click(lambda b: handle_button_click(widgets, connection, m))
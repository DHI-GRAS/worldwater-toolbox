# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:02:11 2022

@author: rmgu
"""

from openeo_processes.utils import process
from worldwater_toolbox.prepare_terrain_shadow_mask import terrain_shadow_mask


@process
def terrain_shadow_mask_process():
    return TerrainShadowMaskProcess


class TerrainShadowMaskProcess:

    @staticmethod
    def exec_np(dem, sun_zenith, sun_azimuth):
        pass

    @staticmethod
    def exec_xar(dem, sun_zenith, sun_azimuth):
        # Get dem transform and resolution from datacube metadata
        dem_transform = None
        dem_resolution = None
        terrain_shadow_mask(dem, sun_zenith, sun_azimuth, dem_transform, dem_resolution)

    @staticmethod
    def exec_da(dem, sun_zenith, sun_azimuth):
        pass

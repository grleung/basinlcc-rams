# script to Download ERA-5 Data for Use in RAMS simulations
# Author: Peter J Marinescu

# Import necessary libraries
import cdsapi
import numpy as np
import os
############################################################################################
# This is for 3D DATA

# Specific save directory and filenames for saving
savedir = '/nobackup/grleung1/basinlcc/b2_test/ERA5/'
savename_p = 'ERA5_20190911_Pressure'
savename_l = 'ERA5_20190911_Land'
savename_f = 'ERA5_20190911'
filetype = 'grib'
dir_finf = savedir + savename_f + '.' + filetype

# Sepcific lat/lon bounding box
las = -6; lan = 9 # Latitude bounds (South ; North)
low = 106; loe = 120 # Longitude bounds (West ; East)

# Specify date/dates 
YY = '2019'
MM = '09'
#DD = ['13','14',]
DD = ['11','12','13']

############################################################################################
############################################################################################
# Pull 3D, hourly atmospheric conditions from ERA5 reanalysis

c = cdsapi.Client()
dir_fin1 = savedir + savename_p + '.' + filetype
c.retrieve(
'reanalysis-era5-pressure-levels',
{
    'product_type': 'reanalysis',
    'format': filetype,
    'variable': ['temperature','u_component_of_wind','v_component_of_wind',
                 'specific_humidity','geopotential',
    ],
    'pressure_level': 'all',
        'year': YY,
        'month': MM,
        'day': DD,
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],          
        'area': [
            lan, low, las,
            loe,
        ],
    },
    dir_fin1)

############################################################################################
############################################################################################
# Pull hourly surface conditions from ERA5 reanalysis

c = cdsapi.Client()
dir_fin2 = savedir + savename_l + '.' + filetype
c.retrieve(
'reanalysis-era5-single-levels',
{
    'product_type': 'reanalysis',
    'format': filetype,
    'variable': ['soil_temperature_level_1','soil_temperature_level_2',
                 'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2',
                 'snow_depth','snow_depth_water_equivalent',
    ],
        'year': YY,
        'month': MM,
        'day': DD,
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],          
        'area': [
            lan, low, las,
            loe,
        ],
    },
    dir_fin2)

os.system("cat "+dir_fin1+" "+dir_fin2+" > "+dir_finf)

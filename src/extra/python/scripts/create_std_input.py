# -*- coding: utf-8 -*-s
import numpy as np
from calendar_calc import day_number_to_date
from netCDF4 import Dataset, date2num
import pdb
import create_timeseries as cts




sigma_full = 
sigma_half = 

lats = 
lons = 

latbs = 
lonbs = 


tstd_data = np.zeros(40, 64, 128)


#Find grid and time numbers

nlon=len(lons)
nlat=len(lats)

nlonb=len(lonbs)
nlatb=len(latbs)

npfull=len(sigma_full)
nphalf=len(sigma_half)



#Output it to a netcdf file. 
file_name='ml_std_input.nc'
variable_name='tstd'

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['npfull']=npfull
number_dict['nphalf']=nphalf
number_dict['ntime']=0

time_arr = None
time_units='days since 0000-01-01 00:00:00.0'

cts.output_to_file(tstd_data,lats,lons,latbs,lonbs,sigma_full,sigma_half,time_arr,time_units,file_name,variable_name,number_dict)


# -*- coding: utf-8 -*-s
import numpy as np
from netCDF4 import Dataset, date2num
import pdb
import create_timeseries as cts
import xarray as xar

basic_dataset = xar.open_dataset('/home/links/sit204/isca_data_intel/ml_test_no_ml_1/run0001/atmos_monthly.nc', decode_times=False)

pfull = basic_dataset['pfull']
phalf = basic_dataset['phalf']

bottom_pressure = phalf.max()

sigma_full = pfull/bottom_pressure
sigma_half = phalf/bottom_pressure

lats = basic_dataset['lat'].values
lons = basic_dataset['lon'].values

latbs = basic_dataset['latb'].values
lonbs = basic_dataset['lonb'].values

variable_name_list = ['tstd', 'qstd']

tstd_data = np.zeros((40, 64, 128))+10.
qstd_data = np.zeros((40, 64, 128))+2.

output_data_dict = {'tstd': tstd_data, 'qstd':qstd_data}

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



cts.output_to_file(output_data_dict,lats,lons,latbs,lonbs,sigma_full,sigma_half,time_arr,time_units,file_name,variable_name_list,number_dict)


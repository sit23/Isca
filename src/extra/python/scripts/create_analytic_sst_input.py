import create_timeseries as cts
import gauss_grid as gg
import numpy as np
import pdb

def create_grid(t_res=42):

    t_res_dict = {42:[128,64],
                  85:[256,128],
                 }

    n_lon = t_res_dict[t_res][0]
    n_lat = t_res_dict[t_res][1]

    delta_lon = (360./n_lon)
    longitudes_in  = np.arange(0., 360., delta_lon)

    lonbs = np.arange(-delta_lon/2., 360., delta_lon)

    latitudes_in, latbs_orig  = gg.gaussian_latitudes(int(n_lat/2))
    latbs = [latbs_orig[i][0] for i in range(len(latitudes_in))]
    latbs.append(latbs_orig[-1][1])

    latbs = np.asarray(latbs)

    return longitudes_in, latitudes_in, lonbs, latbs

def prescribe_constant_sst_everywhere(sst_array, sst_value_to_use):

    sst_array = np.zeros_like(sst_array) + sst_value_to_use

    return sst_array

if __name__=="__main__":

    #Output it to a netcdf file. 
    variable_name='analytic_t_surf_1'
    file_name=variable_name+'_v1.nc'
    sst_value = 273.

    lons, lats, lonbs, latbs = create_grid()

    nlat = len(lats)
    nlon = len(lons)
    nlonb = len(lonbs)
    nlatb = len(latbs)    

    sst_in = np.zeros((12, nlat, nlon))

    sst_out = prescribe_constant_sst_everywhere(sst_in, sst_value)

    p_full = None
    p_half = None

    time_arr = np.arange(0,360,30)

    number_dict={}
    number_dict['nlat']=nlat
    number_dict['nlon']=nlon
    number_dict['nlatb']=nlatb
    number_dict['nlonb']=nlonb    
    number_dict['ntime']=len(time_arr)

    time_units='days since 0000-01-01 00:00:00.0'

    cts.output_to_file(sst_out,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict)
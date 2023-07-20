import xarray as xar
import numpy as np

def convert_column_to_input(file_name_in, nlon_out, file_name_out, p_surf):

    dataset = xar.open_dataset(file_name_in, decode_times=False)

    lats_out  = dataset['lat'].values
    latbs_out = dataset['latb'].values

    delta_lon =  (360./nlon_out)

    lons_out  = np.arange(0., 360., delta_lon)
    lonbs_out = np.arange(-delta_lon/2., 360., delta_lon)

    dataset_out = xar.Dataset(coords=dataset.coords)
    dataset_out = dataset_out.drop(('time', 'nv'))

    dataset_out['lat'] = lats_out
    dataset_out['latb'] = latbs_out
    dataset_out['lon'] = lons_out
    dataset_out['lonb'] = lonbs_out

    for var_name in ['lat', 'latb', 'lon', 'lonb']:
        dataset_out[var_name].attrs = dataset[var_name].attrs

    nlat   = dataset_out['lat'].shape[0]
    nlon   = dataset_out['lon'].shape[0]
    npfull = dataset_out['pfull'].shape[0]

    ps    = np.zeros((nlat, nlon)) +p_surf
    sphum = np.zeros((npfull, nlat, nlon))
    t     = np.zeros((npfull, nlat, nlon))
    u     = np.zeros((npfull, nlat, nlon))
    v     = np.zeros((npfull, nlat, nlon))

    av_t = dataset['temp'].mean(('time', 'lon')).values
    for lon_idx in range(nlon):
        t[:,:,lon_idx] = av_t

    dataset_out['ps']    = (('lat', 'lon'), ps)
    dataset_out['sphum'] = (('pfull', 'lat', 'lon'), sphum)
    dataset_out['t']     = (('pfull', 'lat', 'lon'), t)
    dataset_out['u']     = (('pfull', 'lat', 'lon'), u)
    dataset_out['v']     = (('pfull', 'lat', 'lon'), v)

    dataset_out.to_netcdf(file_name_out)    

if __name__=="__main__":

    month_num = 120

    for del_int_flux in [-8.]:
        base_exp = f'small_giant_planet_column_40_levels_3_bar_4.5_sh_lw_int_flux_del_int_flux_{del_int_flux}'

        file_name_in = f'/home/links/sit204/isca_data_intel/{base_exp}/run{month_num:04d}/atmos_3monthly.nc'

        p_surf = 3e5

        nlon_out = 256
        file_name_out = f'col_{month_num}_{del_int_flux}.nc'

        convert_column_to_input(file_name_in, nlon_out, file_name_out, p_surf)

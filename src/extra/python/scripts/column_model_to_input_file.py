import xarray as xar
import numpy as np

def convert_column_to_input(file_name_in, nlon_out, file_name_out, p_surf, noise_amp = 0.001):

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
        t[:,:,lon_idx] = av_t + np.random.rand(npfull, nlat)*noise_amp

    dataset_out['ps']    = (('lat', 'lon'), ps)
    dataset_out['sphum'] = (('pfull', 'lat', 'lon'), sphum)
    dataset_out['t']     = (('pfull', 'lat', 'lon'), t)
    dataset_out['u']     = (('pfull', 'lat', 'lon'), u)
    dataset_out['v']     = (('pfull', 'lat', 'lon'), v)

    dataset_out.to_netcdf(file_name_out)   

    return dataset_out 

def calculate_gs_balanced_winds(dataset_for_analysing, radius , omega, R):

    ug = xar.zeros_like(dataset_for_analysing['t'])
    int_dt_dy = xar.zeros_like(dataset_for_analysing['t'])    
    npfull = dataset_for_analysing['pfull'].shape[0]
    pfull = dataset_for_analysing['pfull'].values
    f = 2.*omega * np.sin(np.deg2rad(dataset_for_analysing['lat']))

    dt_dy = (1./radius)*dataset_for_analysing['t'].differentiate('lat')*180./np.pi

    dataset_for_analysing['dtdy'] = (dt_dy.dims, dt_dy.values)
    dataset_for_analysing['f'] = (f.dims, f.values)

    for p_idx in range(1, npfull):
        int_dt_dy[p_idx,...] = int_dt_dy[p_idx-1,...] + (R/f)*0.5*(dt_dy[p_idx,...]/pfull[p_idx]+dt_dy[p_idx-1,...]/pfull[p_idx-1])*(pfull[p_idx]-pfull[p_idx-1])

    for p_idx in range(1, npfull):
        ug[p_idx,...] = ug[0,...] + int_dt_dy[p_idx,...]

    dataset_for_analysing['int_dt_dy'] = (int_dt_dy.dims, int_dt_dy.values)
    dataset_for_analysing['ug'] = (ug.dims, ug.values)

if __name__=="__main__":

    month_num = 240

    for del_int_flux in [-8., 0., 8.]:
        base_exp = f'small_giant_planet_column_40_levels_3_bar_4.5_sh_lw_int_flux_del_int_flux_{del_int_flux}'

        file_name_in = f'/home/links/sit204/isca_data_intel/{base_exp}/run{month_num:04d}/atmos_3monthly.nc'

        p_surf = 3e5

        nlon_out = 6
        file_name_out = f'col_{month_num}_{del_int_flux}_zm.nc'

        radius = 25559.0e3
        omega  = 1.0124e-4
        R      = 3750.0

        dataset_for_analysing = convert_column_to_input(file_name_in, nlon_out, file_name_out, p_surf)

        calculate_gs_balanced_winds(dataset_for_analysing, radius , omega, R)

        dataset_for_analysing.to_netcdf(f'ug_{file_name_out}')
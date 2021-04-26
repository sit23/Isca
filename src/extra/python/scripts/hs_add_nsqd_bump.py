import numpy as np
import matplotlib.pyplot as plt
import xarray as xar
from tqdm import tqdm
import gauss_grid as gg
import create_timeseries as cts
from scipy.interpolate import CubicSpline
import pdb

def set_model_params(kappa = 2./7., g=9.81, rdgas=287.04, delh=60., t_zero=315., t_strat=200., delv=10.):

    model_params= {'g':g, 'kappa':kappa, 'rdgas':rdgas, 'pref':1000., 
    'delh':delh, 'delv':delv, 't_zero':t_zero , 't_strat':t_strat }

    return model_params

def offline_hs(dataset, model_params_internal, P00 = 1.e3, p_trop  = 1.e4, alpha = 2./7,  ka = -40., ks =  -4., kf = -1., eps=0.0, sigma_b=0.7, chai_vallis=False, optional_subscript=''):

    lat_array = dataset['lat'].values
    pfull_array = dataset['pfull'].values    

    delh=model_params_internal['delh']
    delv=model_params_internal['delv']    
    t_zero=model_params_internal['t_zero']    
    t_strat=model_params_internal['t_strat']        

    sin_lat = np.sin(np.deg2rad(lat_array))
    sin_lat_2 = sin_lat**2.
    cos_lat = np.cos(np.deg2rad(lat_array))
    cos_lat_2 = cos_lat**2.    

    nlat = lat_array.shape[0]
    nlev = pfull_array.shape[0]    

    t_star = np.zeros((nlat))
    tstr   = np.zeros((nlev, nlat))  
    teq    = np.zeros((nlev, nlat))
    teq_final    = np.zeros((nlev, nlat))


    t_star = t_zero - delh*sin_lat_2 - eps*sin_lat
    tstr   = t_strat - eps*sin_lat

    tks = -1./(86400*ks)
    tka = -1./(86400*ka)


    tcoeff = (tks-tka)/(1.0-sigma_b)
    pref = P00
    rps  = 1./P00

    for k in range(nlev):

        p_norm = pfull_array[k]/pref
        if chai_vallis:
            the    = t_star - delv*np.log(p_norm)            
        else:
            the    = t_star - delv*cos_lat_2*np.log(p_norm)
        teq[k, :] = the*(p_norm)**model_params_internal['kappa']
        teq_final[k, :] = np.maximum( teq[k,:], tstr[:] )

    dataset[f'teq_no_min{optional_subscript}'] = (('pfull', 'lat'), teq)
    dataset[f'teq{optional_subscript}'] = (('pfull', 'lat'), teq_final)

    theta = dataset[f'teq{optional_subscript}']*((pref/dataset['pfull'])**model_params_internal['kappa'])
    dataset[f'theta{optional_subscript}'] = (theta.dims, theta.values)


    brunt_vas_freq(dataset, model_params_internal)

def diffz(data, vcoord, axis=None):
  
    dquantity_dp = np.zeros_like(data)

    dquantity_dp[0,:] = (data[1,:] - data[0,:]) / (vcoord[1]- vcoord[0]) #forward difference at the top

    dquantity_dp[-1,:] = (data[-1,:] - data[-2,:]) / (vcoord[-1]- vcoord[-2]) #backward difference at the bottom

    for n_int in range(1, len(vcoord)-1):
        dquantity_dp[n_int,:] = (data[n_int+1,:] - data[n_int-1,:]) / (vcoord[n_int+1]- vcoord[n_int-1])
    
    return dquantity_dp

def brunt_vas_freq(dataset, model_params, temp_name='teq', name_out='nsqd', theta_name='theta'):

    grav = model_params['g']

    #Find relevant axes:
    pfull_loc = [x=='pfull' for x in dataset[theta_name].dims]

    pfull_idx = np.where(pfull_loc)[0][0]

    d_theta_dp = diffz(dataset[theta_name].values, dataset.pfull.values*100., axis=pfull_idx)

    nsqd_prefactor =  (-grav**2./dataset[theta_name])*(dataset.pfull * 100. / (model_params['rdgas']*dataset[temp_name]))

    nsqd = nsqd_prefactor * d_theta_dp

    dataset[name_out]=(nsqd.dims,nsqd)

def theta_from_nsqd(dataset, model_params, temp_name='teq', nsqd_name='nsqd_total'):

    pfull_array = dataset['pfull'].values
    lat_array = dataset['lat'].values
    nlat = lat_array.shape[0]

    theta_int = np.zeros_like(dataset['theta'])

    rho = (dataset.pfull * 100. / (model_params['rdgas']*dataset[temp_name]))

    expression_to_integrate = -1.*dataset[nsqd_name].values/((model_params['g']**2.)*rho.values)

    theta_int[0, :] = dataset['theta'][0,:].values

    for lat_idx in range(nlat):
        for n_int in range(1, pfull_array.shape[0]):
            theta_int[n_int, lat_idx] = theta_int[n_int-1, lat_idx] * np.exp(0.5*(expression_to_integrate[n_int, lat_idx]+expression_to_integrate[n_int-1, lat_idx])*100.*(pfull_array[n_int]-pfull_array[n_int-1]))

    dataset['theta_int_total'] = (('pfull', 'lat'), theta_int)

    temp_int_total = dataset['theta_int_total'] * (dataset['pfull']/model_params['pref'])**model_params['kappa']

    dataset['temp_int_total'] = (('pfull', 'lat'), temp_int_total)


def construct_nsqd(dataset):

    p_peak  = 600.
    p_water = 650.
    p_decay = 50.
    p_array = dataset['pfull'].values
    nsqd_max=0.00025

    nsqd_construct = np.zeros_like(dataset['teq'])

    for lat_idx in range(dataset['lat'].shape[0]):
        for p_idx in range(dataset['pfull'].values.shape[0]):
            p_value = p_array[p_idx]
            if p_value>=p_water:
                nsqd_construct[p_idx, lat_idx] = 0.
            elif p_value<p_water and p_value>=p_peak:
                nsqd_construct[p_idx, lat_idx] = nsqd_max * np.arctanh(np.tanh(1.)*(p_value - p_water)/(p_peak-p_water))
            elif p_value<p_peak:
                nsqd_construct[p_idx, lat_idx] = nsqd_max * np.exp((p_value-p_peak)/p_decay)

    dataset['nsqd_construct'] = (('pfull', 'lat'), nsqd_construct)
    dataset['nsqd_total'] = dataset['nsqd']+dataset['nsqd_construct']

def modify_final_temp_profile_to_fix_emission_temp(dataset, model_params, pval_chosen, pval_chosen_min_eq, pval_chosen_min_polar, chai_vallis=False):

    kappa = model_params['kappa']
    pref  = model_params['pref']
    deltav= model_params['delv']
    deltah= model_params['delh']    

    min_lat = np.min(np.abs(dataset.lat.values))
    nearest_pval_chosen = dataset['pfull'].sel(pfull=pval_chosen, method='nearest').values
    # nearest_pval_chosen_min_eq = dataset['pfull'].sel(pfull=pval_chosen_min_eq, method='nearest').values
    # nearest_pval_chosen_min_polar = dataset['pfull'].sel(pfull=pval_chosen_min_polar, method='nearest').values


    teq_new = dataset['temp_int_total'].sel(lat=min_lat).sel(pfull=nearest_pval_chosen).values

    new_t0_param = teq_new*(pref/nearest_pval_chosen)**kappa + deltav*np.log(nearest_pval_chosen/pref) + deltah*(np.sin(np.deg2rad(min_lat)))**2.

    model_params_mod = set_model_params(t_zero=new_t0_param)

    offline_hs(dataset, model_params_mod, optional_subscript='_fit', chai_vallis=chai_vallis)

    teq_final = dataset['temp_int_total'].copy(deep=True).values
    teq_final_cs = dataset['temp_int_total'].copy(deep=True).values


    where_to_apply_fit = np.where(dataset['pfull'].values<=nearest_pval_chosen)[0]

    teq_final[where_to_apply_fit,:] = dataset['teq_fit'].values[where_to_apply_fit,:]

    dataset['teq_final'] = ((dataset['temp_int_total'].dims), teq_final)

    pfull_array = dataset['pfull'].values
    lat_array = dataset['lat'].values

    for lat_idx in range(nlat):

        nearest_pval_chosen_min = pval_chosen_min_eq + (pval_chosen_min_polar-pval_chosen_min_eq)*(np.sin(np.deg2rad(lat_array[lat_idx])))**2.

        where_to_fit = np.where(dataset['pfull'].values<=nearest_pval_chosen)[0]
        where_to_train= np.where(np.logical_and(nearest_pval_chosen_min<=dataset['pfull'].values, dataset['pfull'].values<=nearest_pval_chosen))[0]

        cs_obj = CubicSpline(pfull_array[where_to_train][::-1], dataset['temp_int_total'][where_to_train, lat_idx].values[::-1])
        extended_temps = cs_obj(pfull_array[where_to_fit][::-1])
        teq_final_cs[where_to_fit, lat_idx] = np.maximum(extended_temps[::-1], model_params['t_strat'])

    dataset['teq_final_cs'] = ((dataset['temp_int_total'].dims), teq_final_cs)

    lat_plot=0.
    plt.figure()
    plt.plot(dataset['teq'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='original')
    plt.plot(dataset['temp_int_total'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='modified with bump')    
    plt.plot(dataset['teq_fit'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq')
    plt.plot(dataset['teq_final'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq final', linestyle='-')
    plt.plot(dataset['teq_final_cs'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq final cs', linestyle='-')
    plt.plot([200., 300.], [pval_chosen_min_eq, pval_chosen_min_eq])
    plt.plot([200., 300.], [nearest_pval_chosen, nearest_pval_chosen])    
    plt.legend()
    plt.ylim([1000., 0.])    


    lat_plot=90.
    plt.figure()
    plt.plot(dataset['teq'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='original')
    plt.plot(dataset['temp_int_total'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='modified with bump')    
    plt.plot(dataset['teq_fit'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq')
    plt.plot(dataset['teq_final'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq final', linestyle='-')
    plt.plot(dataset['teq_final_cs'].sel(lat=lat_plot, method='nearest').values, dataset['pfull'].values, label='fitted teq final cs', linestyle='-')
    plt.plot([200., 300.], [pval_chosen_min_polar, pval_chosen_min_polar])
    plt.plot([200., 300.], [nearest_pval_chosen, nearest_pval_chosen])        
    plt.legend()
    plt.ylim([1000., 0.])    

    plt.figure()
    plt.plot(dataset['teq_fit'].sel(lat=0., method='nearest').values, dataset['pfull'].values, label='original')
    plt.plot(dataset['teq_fit'].sel(lat=90., method='nearest').values, dataset['pfull'].values, label='original')    
    plt.ylim([1000., 0.])    


def output_isca_input_files(dataset):

    ntime=12
    npfull = dataset['pfull'].values.shape[0]
    nphalf = dataset['phalf'].values.shape[0]

    nlon = dataset['lon'].values.shape[0]
    nlat = dataset['lat'].values.shape[0]

    latitudes=dataset['lat'].values
    longitudes=dataset['lon'].values    
    latitude_bounds = dataset['latb'].values
    longitude_bounds = dataset['lonb'].values    

    time_arr_adj=np.arange(15,360,360/ntime)

    for temp_name in ['teq', 'teq_final_cs']:
        teq_array = np.zeros((ntime, npfull, nlat, nlon))

        for month_tick in range(12):
            for lon_tick in range(nlon):
                teq_array[month_tick, :, :, lon_tick]  = dataset[temp_name].values

        teq_array=teq_array[:,::-1,:,:]

        p_full=dataset['pfull'].values[::-1]
        p_half=dataset['phalf'].values[::-1]

        #Output it to a netcdf file. 
        variable_name=f'earth_chai_vallis_{temp_name}'
        file_name=f'{variable_name}.nc'


        number_dict={}
        number_dict['nlat']=nlat
        number_dict['nlon']=nlon
        number_dict['nlatb']=nlat+1
        number_dict['nlonb']=nlon+1 
        number_dict['npfull']=npfull
        number_dict['nphalf']=nphalf
        number_dict['ntime']=ntime


        time_units='days since 0000-01-01 00:00:00.0'

        cts.output_to_file(teq_array,latitudes,longitudes,latitude_bounds,longitude_bounds,p_full,p_half,time_arr_adj,time_units,file_name,variable_name,number_dict)
    

if __name__=="__main__":

    nlat=128
    nlon=256
    npfull = 80

    latitudes, latitude_bounds_2  = gg.gaussian_latitudes(int(nlat/2))
    latitude_bounds = [latitude_bound[0] for latitude_bound in latitude_bounds_2] + [latitude_bounds_2[-1][1]]

    longitudes = np.linspace(0., 360., nlon, endpoint=False)
    delta_lon = longitudes[1]-longitudes[0]
    longitude_bounds = [lon_val-(0.5*delta_lon) for lon_val in longitudes] + [np.max(longitudes)+(0.5*delta_lon)]
    time_arr_adj=np.arange(15,360,30)

    lon_array_2d, lat_array_2d = np.meshgrid(longitudes, latitudes)

    lat_array = latitudes
    # lat_array  = np.linspace(-90., 90., num=128)
    phalf_array = np.linspace(1200., 0., num=npfull+1, endpoint=True)
    pfull_array = [0.5*(phalf_array[pidx]+phalf_array[pidx-1]) for pidx in range(1,npfull+1)] 

    model_params=set_model_params()

    dataset = xar.Dataset(coords=dict(
        lon=('lon', longitudes),
        lat=('lat', lat_array),
        lonb=('lonb', longitude_bounds),
        latb=('latb', latitude_bounds),   
        phalf=('phalf', phalf_array),               
        pfull=('pfull', pfull_array)))

    for do_chai_vallis in [True]:

        #construct original basic state
        offline_hs(dataset, model_params, chai_vallis=do_chai_vallis)

        #construct nsqd anomaly
        construct_nsqd(dataset)

        #construct theta and temperatures from nsqd assuming original temperatures
        theta_from_nsqd(dataset, model_params, temp_name='teq', nsqd_name='nsqd_total')

        # plt.figure()
        # for lat_value in [0.]:
        #     plt.plot(dataset['nsqd_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=f'{lat_value} total', marker='x')            

        iteration_dict = {}

        max_iterations=100

        for iteration_step in tqdm(range(max_iterations)):
            #construct theta from nsqd assuming original temperatures
            theta_from_nsqd(dataset, model_params, temp_name='temp_int_total', nsqd_name='nsqd_total')

            #construct new nsqd using new theta to compare with nsqd we set out to reproduce
            brunt_vas_freq(dataset, model_params, temp_name='temp_int_total', name_out='nsqd_int_total', theta_name='theta_int_total')

            iteration_dict[iteration_step] = dataset['nsqd_int_total'].values

            # for lat_value in [0.]:
            #     plt.plot(dataset['nsqd_int_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=f'{lat_value} total', marker='x')      

        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()

        modify_final_temp_profile_to_fix_emission_temp(dataset, model_params, 470., 250., 400., chai_vallis=do_chai_vallis) 
        output_isca_input_files(dataset)


        plt.figure()
        plt.plot(iteration_dict[0][:,64], dataset['pfull'].values)        
        plt.plot(iteration_dict[max_iterations-1][:,64], dataset['pfull'].values)        
        plt.plot(iteration_dict[0][:,64]-iteration_dict[max_iterations-1][:,64], dataset['pfull'].values)
        plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )


        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['teq'], cmap='RdBu_r', levels=np.linspace(190., 320., num=25, endpoint=True))
        plt.ylim([1000., 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis}')

        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['teq_final_cs'], cmap='RdBu_r', levels=np.linspace(190., 320., num=25, endpoint=True))
        plt.ylim([1000., 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis} with bump added')

        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['teq_final_cs'].values-dataset['teq'].values, cmap='RdBu_r', levels=30)
        plt.ylim([1000., 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis} with bump added minus original')



        # plt.figure()
        # plt.contourf(lat_array, pfull_array, dataset['theta'], cmap='RdBu_r', levels=25)
        # plt.ylim([1000., 0.])    
        # plt.colorbar(extend='both')    
        # plt.title(f'{do_chai_vallis}')

        # plt.figure()
        # plt.contourf(lat_array, pfull_array, dataset['nsqd'], cmap='RdBu_r', levels=25)
        # plt.ylim([1000., 0.])    
        # plt.colorbar(extend='both')    
        # plt.title(f'{do_chai_vallis}')        

        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['teq'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()

        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['theta'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()

        plt.figure()
        for lat_value in [90., 45., 0.]:
            plt.plot(dataset['nsqd_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        plt.legend()
        plt.title('nsqd with idealised bump added')

        plt.figure()
        for lat_value in [90., 45., 0.]:
            plt.plot(dataset['nsqd_int_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        plt.legend()
        plt.title('nsqd with idealised bump added constructed from new temps')


        plt.figure()
        for lat_value in [0., 90.]:
            plt.plot(dataset['nsqd_int_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=f'{lat_value} int total', marker='x')
            plt.plot(dataset['nsqd_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=f'{lat_value} total', marker='x')            
        plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        plt.legend()
        plt.title('nsqd with idealised bump compared with that constructed from new temps')



        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['theta_int_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()


        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['temp_int_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()


        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['nsqd_construct'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()

        # plt.figure()
        # for lat_value in [90., 45., 0.]:
        #     plt.plot(dataset['nsqd_total'].sel(lat=lat_value, method='nearest'), dataset['pfull'].values, label=lat_value, marker='x')
        # plt.ylim([dataset['pfull'].max(), dataset['pfull'].min()] )
        # plt.legend()
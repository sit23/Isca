import numpy as np
import matplotlib.pyplot as plt
import xarray as xar
from tqdm import tqdm
import gauss_grid as gg
import create_timeseries as cts
from scipy.interpolate import CubicSpline
from jupiter_hs_bump_vert_levels import calc_merged_levels
import pdb

def set_model_params_earth(kappa = 2./7., g=9.81, rdgas=287.04, delh=60., t_zero=315., t_strat=200., delv=10., rvgas=461.50):

    model_params= {'g':g, 'kappa':kappa, 'rdgas':rdgas, 'rvgas':rvgas, 'pref':1000., 
    'delh':delh, 'delv':delv, 't_zero':t_zero , 't_strat':t_strat }

    return model_params

def set_model_params_jupiter(kappa = 2./7., g=26.0, rdgas=3605.38, delh=15., t_zero=170., t_strat=110., delv=5., rvgas=461.50):

    model_params= {'g':g, 'kappa':kappa, 'rdgas':rdgas, 'rvgas':rvgas, 'pref':1000., 
    'delh':delh, 'delv':delv, 't_zero':t_zero , 't_strat':t_strat }

    return model_params    

def set_model_params_small_jupiter(kappa = 2./7., g=7.77, rdgas=3605.38, delh=15., t_zero=170., t_strat=110., delv=5., rvgas=461.50):

    model_params= {'g':g, 'kappa':kappa, 'rdgas':rdgas, 'rvgas':rvgas, 'pref':1000., 
    'delh':delh, 'delv':delv, 't_zero':t_zero , 't_strat':t_strat }

    return model_params        

def offline_lian_showman(dataset, model_params_internal, P00 = 1.e3, p_trop  = 150., alpha = 2./7,  ka = -40., ks =  -4., kf = -1., eps=0.0, sigma_b=0.7, chai_vallis=False, optional_subscript='', mode='chai_vallis', perturb_mode='B', fixed_temp=490., fixed_temp_pressure_bar=30., perturbation_mode_b_amplitude=2.6):

    lat_array = dataset['lat']
    pfull_array = dataset['pfull'] 

    delh=model_params_internal['delh']
    delv=model_params_internal['delv']    
    t_zero=model_params_internal['t_zero']    
    t_strat=model_params_internal['t_strat']        

    sin_lat = np.sin(np.deg2rad(lat_array))
    sin_lat_2 = sin_lat**2.
    cos_lat = np.cos(np.deg2rad(lat_array))
    cos_lat_2 = cos_lat**2.    

    cos_8_lat = np.cos(8.*np.deg2rad(lat_array))
    cos_8_lat_2 = cos_8_lat**2.

    nlat = lat_array.shape[0]
    nlev = pfull_array.shape[0]    

    background_t_ref = xar.zeros_like(pfull_array)


    if mode=='chai_vallis':
        G_of_p = xar.zeros_like(pfull_array)

        G_of_p = 1./(1.+(pfull_array/p_trop)**2.)

        background_t_ref = (t_strat*G_of_p) + ((1.-G_of_p)*(t_zero-(delv*np.log(pfull_array/p_trop))))*(pfull_array/p_trop)**model_params_internal['kappa']  
              
        delta_t_horiz = delh*((1./3.)-sin_lat_2)
    elif mode=='lian_showman': 
        G_of_p = xar.zeros_like(pfull_array)

        use_tanh=False
        interp_temp=True
        interp_ptemp=False

        if use_tanh:
            G_of_p = 0.5*(1.+np.tanh((p_trop-pfull_array)/(0.1*p_trop)))
            G_of_p = 0.5*(1.+np.tanh((950.-pfull_array)/(0.1*550.)))


            background_t_ref = (t_strat*G_of_p) + ((1.-G_of_p)*(600.))*(pfull_array/(100.*1000.))**model_params_internal['kappa']  
            background_t_ref_fit = xar.zeros_like(pfull_array)+background_t_ref_fit
        elif interp_temp:
            strat_pfull_points = pfull_array[np.where(pfull_array<400.)[0]].values
            trop_pfull_points = pfull_array[np.where(pfull_array>1500.)[0]].values

            fitting_points = pfull_array[np.where(np.logical_not(np.logical_or(pfull_array<400., pfull_array>1500.)))[0]]

            strat_temps = np.zeros_like(strat_pfull_points) + t_strat
            trop_temps = 600.*(trop_pfull_points/(100.*1000.))**model_params_internal['kappa']  


            combined_pfull = np.concatenate((trop_pfull_points, strat_pfull_points))
            combined_temp = np.concatenate((trop_temps, strat_temps))


            cs_obj = CubicSpline(np.log(combined_pfull[::-1]), combined_temp[::-1])
            background_t_ref_fit = cs_obj(np.log(fitting_points[::-1]))

            background_t_ref = xar.zeros_like(pfull_array)

            background_t_ref[np.where(pfull_array<400.)[0]] = strat_temps
            background_t_ref[np.where(pfull_array>1500.)[0]] = trop_temps
            background_t_ref[np.where(np.logical_not(np.logical_or(pfull_array<400., pfull_array>1500.)))[0]] = background_t_ref_fit[::-1]
        elif interp_ptemp:
            strat_pfull_points = pfull_array[np.where(pfull_array<400.)[0]].values
            trop_pfull_points = pfull_array[np.where(pfull_array>1500.)[0]].values

            fitting_points = pfull_array[np.where(np.logical_not(np.logical_or(pfull_array<400., pfull_array>1500.)))[0]]

            strat_temps =  t_strat*((100.*1000.)/strat_pfull_points)**model_params_internal['kappa']  
            trop_temps = np.zeros_like(trop_pfull_points) + 600.

            combined_pfull = np.concatenate((trop_pfull_points, strat_pfull_points))
            combined_temp = np.concatenate((trop_temps, strat_temps))


            cs_obj = CubicSpline(np.log(combined_pfull[::-1]), np.log(combined_temp[::-1]))
            background_t_ref_fit = np.exp(cs_obj(np.log(fitting_points[::-1])))

            background_t_ref = xar.zeros_like(pfull_array)

            background_t_ref[np.where(pfull_array<400.)[0]] = strat_temps*(strat_pfull_points/(100.*1000.))**model_params_internal['kappa']  
            background_t_ref[np.where(pfull_array>1500.)[0]] = trop_temps*(trop_pfull_points/(100.*1000.))**model_params_internal['kappa']  
            background_t_ref[np.where(np.logical_not(np.logical_or(pfull_array<400., pfull_array>1500.)))[0]] = background_t_ref_fit[::-1]*(fitting_points/(100.*1000.))**model_params_internal['kappa']  

    elif mode=='stevo': 
        G_of_p = xar.zeros_like(pfull_array)

        background_t_ref_trop = (fixed_temp)*(pfull_array/(fixed_temp_pressure_bar*1000.))**model_params_internal['kappa']  
        background_strat_t = xar.zeros_like(pfull_array) + t_strat

        background_t_ref = np.maximum(background_t_ref_trop, background_strat_t)

        # Use mode B and specified sphum in HS to look at impact. Nsqd varies with belts and zones?!
        # 

    if perturb_mode=='A':
        delta_t_horiz = 8.*(cos_lat_2)
    elif perturb_mode=='B':
        delta_t_horiz = perturbation_mode_b_amplitude*(cos_8_lat_2)
        delta_t_horiz[np.where(np.abs(lat_array)>78.75)] = 0.
    elif perturb_mode=='C':
        delta_t_horiz = 8.*(1.-cos_lat_2)        

    teq_background = (background_t_ref*(1.+xar.zeros_like(lat_array)))
    teq_perturbations = (delta_t_horiz*(1.+xar.zeros_like(pfull_array)))
    teq =  teq_background + teq_perturbations

    dataset['teq'] = (teq.dims, teq.values)
    dataset['teq_background'] = (teq_background.dims, teq_background.values)
    dataset['teq_perturbations'] = (teq_perturbations.dims, teq_perturbations.values)


    theta = dataset[f'teq']*((P00/dataset['pfull'])**model_params_internal['kappa'])
    dataset[f'theta'] = (theta.dims, theta.values)

    brunt_vas_freq(dataset, model_params_internal)

    # plt.figure()
    # plt.plot(G_of_p, pfull_array/1000.)
    # ax = plt.gca()
    # ax.set_yscale('log')
    # plt.ylim([10., 0.00001])

    # plt.figure()
    # plt.plot(background_t_ref, pfull_array/1000.)
    # ax = plt.gca()
    # ax.set_yscale('log')
    # plt.ylim([10., 0.00001])

    # plt.figure()
    # plt.plot(dataset[f'theta'].sel(lat=0., method='nearest'), pfull_array/1000.)
    # ax = plt.gca()
    # ax.set_yscale('log')
    # plt.ylim([10., 0.00001])

    # plt.figure()
    # plt.plot(dataset[f'nsqd'].sel(lat=0., method='nearest'), pfull_array/1000.)
    # ax = plt.gca()
    # ax.set_yscale('log')
    # plt.ylim([10., 0.00001])


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

    try:
        theta_array = dataset[theta_name]
    except:
        theta = dataset[temp_name]*((model_params['pref']/dataset['pfull'])**model_params['kappa'])
        dataset[theta_name] = (theta.dims, theta.values)
        theta_array = dataset[theta_name]

    grav = model_params['g']

    #Find relevant axes:
    pfull_loc = [x=='pfull' for x in theta_array.dims]

    pfull_idx = np.where(pfull_loc)[0][0]

    d_theta_dp = diffz(theta_array.values, dataset.pfull.values*100., axis=pfull_idx)

    nsqd_prefactor =  (-grav**2./theta_array)*(dataset.pfull * 100. / (model_params['rdgas']*dataset[temp_name]))

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




def output_isca_input_files(dataset, list_of_vars_to_output=['teq'], planet='earth', resolution_name=''):

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

    for temp_name in list_of_vars_to_output:
        teq_array = np.zeros((ntime, npfull, nlat, nlon))

        for month_tick in range(12):
            for lon_tick in range(nlon):
                teq_array[month_tick, :, :, lon_tick]  = dataset[temp_name].values

        teq_array=teq_array[:,::-1,:,:]

        p_full=dataset['pfull'].values[::-1]
        p_half=dataset['phalf'].values[::-1]

        #Output it to a netcdf file. 
        if planet=='small_jupiter':
            variable_name=f's_j_spec_teq_{temp_name}{resolution_name}'
            file_name=f'{variable_name}.nc'
        elif planet=='jupiter':
            variable_name=f'j_spec_teq_{temp_name}{resolution_name}'
            file_name=f'{variable_name}.nc'            
        else:
            variable_name=f'{planet}_spec_teq_{temp_name}{resolution_name}'
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

def rh_offline(dataset, model_params, t_name='temp', sphum_name='sphum'):


    epsilon = model_params['rdgas']/ model_params['rvgas']
    one_minus_epsilon = 1.-epsilon
    
    es_array = saturation_vapour_press(dataset[t_name].load(), model_params)
    
    rh_numerator = dataset[sphum_name]
    rh_denominator = (epsilon * es_array)/(dataset.pfull - one_minus_epsilon * es_array)
    
    rh_denom_denom = (dataset.pfull - one_minus_epsilon * es_array)
    
    #Nb don't need pfull in Pa as the constant we've used for e_0 is in hPa, so appropriate pressure units are hPa. Model will be using pfull in Pa and so es in Pa too.
    
    rh = 100.*rh_numerator / rh_denominator
        
    dataset['rh_offline'] = (dataset[sphum_name].dims, rh)

def saturation_specific_humidity(dataset, model_params, t_name='temp'):

    es_array = saturation_vapour_press(dataset[t_name].load(), model_params)

    dataset['sat_vap_pres'] = (dataset[t_name].dims, es_array)


    epsilon = model_params['rdgas']/ model_params['rvgas']

    sat_specific_hum = (epsilon * dataset['sat_vap_pres']) / (dataset['pfull'] - ((1.-epsilon)*dataset['sat_vap_pres']))

    dataset['sat_sphum'] = (sat_specific_hum.dims, sat_specific_hum)


def saturation_vapour_press(temp, model_params, L_vapour = 2.44e6):

    e_0 = 6.12
    T_0 = 273.0
    
    es = e_0 * np.exp((L_vapour/model_params['rvgas']) * ( (1./T_0) - (1./temp)))
    
    return es

def construct_sphum(dataset, model_params, p_transition=500., lat_varying=True, chosen_lat=0., t_name_sphum='teq'):

    saturation_specific_humidity(dataset, model_params, t_name=t_name_sphum)

    if lat_varying:
        sat_spec_shum_at_p_transition = dataset['sat_sphum'].sel(pfull=p_transition, method='nearest')
    else:
        sat_spec_shum_at_p_transition_single = dataset['sat_sphum'].sel(pfull=p_transition, method='nearest').sel(lat=chosen_lat, method='nearest').values
        sat_spec_shum_at_p_transition = xar.zeros_like(dataset.lat)+sat_spec_shum_at_p_transition_single

    constant_value_below_p_transition = sat_spec_shum_at_p_transition * 0.75

    sphum_array = xar.zeros_like(dataset[t_name_sphum])
    sphum_profile_lower = xar.zeros_like(dataset.pfull)
    sphum_profile_upper = xar.zeros_like(dataset.pfull)    
    sphum_profile_lower[np.where(sphum_profile_lower.pfull>=p_transition)] = 1.0
    sphum_profile_upper[np.where(sphum_profile_upper.pfull<p_transition)] = 1.0

    #first add the part where the q value is constant below the transition level
    sphum_array = constant_value_below_p_transition * sphum_profile_lower

    if lat_varying:
        sphum_array = sphum_array + (dataset['sat_sphum'] * sphum_profile_upper*0.75)
    else:
        upper_sphum_change = (dataset['sat_sphum'].sel(lat=chosen_lat, method='nearest').squeeze())*(xar.zeros_like(dataset.lat)+1.)
        sphum_array = sphum_array + (upper_sphum_change* sphum_profile_upper*0.75)
        pdb.set_trace()

    sphum_array = sphum_array.transpose('pfull', 'lat')

    dataset['sphum'] = (sphum_array.dims, sphum_array)

    rh_offline(dataset, model_params, t_name='teq')

def virt_temp(dataset,model_params, t_name='temp', sphum_name='sphum'):

    virt_fact = (model_params['rvgas']/model_params['rdgas']) - 1.0

    virt_temp = dataset[t_name]*(1.+virt_fact*dataset[sphum_name])
    
    dataset['virt_temp'] = (dataset[t_name].dims, virt_temp)

if __name__=="__main__":

    resolution='T85_small'
    npfull = 80
    planet='jupiter'
    merged_levels=True
    nhigh_merged = 61
    nlow_merged=21

    if merged_levels:
        merged_level_str = f'_{nhigh_merged}_{nlow_merged}'
    else:
        merged_level_str=''

    if resolution=='T85':
        nlat=128
        nlon=256
    elif resolution=='T213':
        nlat=320
        nlon=1024       
    elif resolution=='T213_small':
        nlat=320
        nlon=1
    elif resolution=='T85_small':
        nlat=128
        nlon=1    

    latitudes, latitude_bounds_2  = gg.gaussian_latitudes(int(nlat/2))
    latitude_bounds = [latitude_bound[0] for latitude_bound in latitude_bounds_2] + [latitude_bounds_2[-1][1]]

    if nlon!=1:
        longitudes = np.linspace(0., 360., nlon, endpoint=False)
        delta_lon = longitudes[1]-longitudes[0]
        longitude_bounds = [lon_val-(0.5*delta_lon) for lon_val in longitudes] + [np.max(longitudes)+(0.5*delta_lon)]
    else:
        longitudes = np.asarray([180.])
        longitude_bounds = np.asarray([0., 360.])

    time_arr_adj=np.arange(15,360,30)

    lon_array_2d, lat_array_2d = np.meshgrid(longitudes, latitudes)

    lat_array = latitudes
    # lat_array  = np.linspace(-90., 90., num=128)
    if planet=='earth':
        phalf_array = np.linspace(1200., 0., num=npfull+1, endpoint=True)
        pfull_array = [0.5*(phalf_array[pidx]+phalf_array[pidx-1]) for pidx in range(1,npfull+1)]      
    elif planet=='jupiter':
        if merged_levels:
            phalf_array, pfull_array, zhalf_array, zfull_array = calc_merged_levels(nlow_merged, nhigh_merged)            
        else:
            zhalf_array = np.linspace(0., 5., num=npfull, endpoint=True) #only use npfull so that 0. can be added at the end to make nphalf
            phalf_array = 10.*1200.*np.exp(-zhalf_array)
            phalf_array = np.concatenate((phalf_array, np.array([0.])))        
            zhalf_array = np.concatenate((zhalf_array, np.array([np.inf])))   
            pfull_array = [0.5*(phalf_array[pidx]+phalf_array[pidx-1]) for pidx in range(1,npfull+1)] 
            zfull_array = [0.5*(zhalf_array[pidx]+zhalf_array[pidx-1]) for pidx in range(1,npfull+1)]                     
    elif planet=='small_jupiter':
        if merged_levels:
            phalf_array, pfull_array, zhalf_array, zfull_array = calc_merged_levels(nlow_merged, nhigh_merged)                        
        else:
            zhalf_array = np.linspace(0., 5., num=npfull, endpoint=True) #only use npfull so that 0. can be added at the end to make nphalf
            phalf_array = 10.*1200.*np.exp(-zhalf_array)
            phalf_array = np.concatenate((phalf_array, np.array([0.])))  
            zhalf_array = np.concatenate((zhalf_array, np.array([np.inf])))   
            pfull_array = [0.5*(phalf_array[pidx]+phalf_array[pidx-1]) for pidx in range(1,npfull+1)] 
            zfull_array = [0.5*(zhalf_array[pidx]+zhalf_array[pidx-1]) for pidx in range(1,npfull+1)]             



    if planet=='earth':
        model_params=set_model_params_earth()
    elif planet=='jupiter':
        model_params=set_model_params_jupiter()
    elif planet=='small_jupiter':
        model_params=set_model_params_small_jupiter()        

    dataset = xar.Dataset(coords=dict(
        lon=('lon', longitudes),
        lat=('lat', lat_array),
        lonb=('lonb', longitude_bounds),
        latb=('latb', latitude_bounds),   
        phalf=('phalf', phalf_array),               
        pfull=('pfull', pfull_array)))

    for do_chai_vallis in [True]:

        #construct original basic state
        if planet=='earth':
            offline_hs(dataset, model_params, chai_vallis=do_chai_vallis)
            construct_sphum(dataset, model_params, p_transition=500.)

        elif planet=='jupiter':
            offline_lian_showman(dataset, model_params, mode='stevo', perturb_mode='A')
            construct_sphum(dataset, model_params, p_transition=4000., t_name_sphum='teq_background')

        elif planet=='small_jupiter':
            offline_lian_showman(dataset, model_params, mode='stevo')
            construct_sphum(dataset, model_params, p_transition=4000., t_name_sphum='teq_background')            

        virt_temp(dataset,model_params, t_name='teq')        
        brunt_vas_freq(dataset, model_params, temp_name='virt_temp', name_out='nsqd_virt', theta_name='theta_virt')        
        # output_isca_input_files(dataset, list_of_vars_to_output=['teq', 'virt_temp'], planet=planet, resolution_name=f'_{resolution.lower()}{merged_level_str}_A')
        # output_isca_input_files(dataset, list_of_vars_to_output=['teq', 'virt_temp'], planet=planet, resolution_name=f'_a')

        nsqd_diff = dataset['nsqd_virt'] - dataset['nsqd']
        theta_diff = dataset['theta_virt'] - dataset['theta']
        temp_diff = dataset['virt_temp'] - dataset['teq']


        plt.figure()
        # dataset['nsqd_virt'].sel(lat=0., method='nearest').plot.line(label='virt')
        # dataset['nsqd'].sel(lat=0., method='nearest').plot.line(label='normal')       
        nsqd_diff.sel(lat=0., method='nearest').plot.line(label='diff', marker='x')             
        # plt.legend()

        plt.figure()    
        plt.plot(zfull_array, nsqd_diff.sel(lat=0., method='nearest').values, label='diff', marker='x')             
        plt.legend()


        # plt.figure()
        # dataset['theta_virt'].sel(lat=0., method='nearest').plot.line(label='virt')
        # dataset['theta'].sel(lat=0., method='nearest').plot.line(label='normal')       
        # theta_diff.sel(lat=0., method='nearest').plot.line(label='diff')             
        # plt.legend()        

        # plt.figure()
        # dataset['virt_temp'].sel(lat=0., method='nearest').plot.line(label='virt')
        # dataset['teq'].sel(lat=0., method='nearest').plot.line(label='normal')       
        # temp_diff.sel(lat=0., method='nearest').plot.line(label='diff')             
        # plt.legend()             

        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['teq'], cmap='RdBu_r', levels=25)
        plt.ylim([dataset.pfull.max(), 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis}')

        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['virt_temp'], cmap='RdBu_r', levels=25)
        plt.ylim([dataset.pfull.max(), 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis} with bump added')

        plt.figure()
        plt.contourf(lat_array, pfull_array, dataset['virt_temp'].values-dataset['teq'].values, cmap='RdBu_r', levels=30)
        plt.ylim([dataset.pfull.max(), 0.])    
        plt.colorbar(extend='both')    
        plt.title(f'{do_chai_vallis} with bump added minus original')        

        # plt.figure()
        # plt.contourf(lat_array, pfull_array, dataset['nsqd_virt'].values-dataset['nsqd'].values, cmap='RdBu_r', levels=30)
        # plt.ylim([dataset.pfull.max(), 0.])    
        # plt.colorbar(extend='both')    
        # plt.title(f'nsqd {do_chai_vallis} with bump added minus original')    

        # plt.figure()
        # plt.contourf(lat_array, pfull_array, dataset['sat_sphum'], cmap='RdBu_r', levels=25)
        # plt.ylim([dataset.pfull.max(), 0.])    
        # plt.colorbar(extend='both')    
        # plt.title(f'{do_chai_vallis}')        


        # background_t_ref_trop_lian_showman_2008 = (600.)*(dataset['pfull']/(100.*1000.))**model_params['kappa']  
        # background_t_ref_trop_lian_sugiyama2014 = (490.)*(dataset['pfull']/(30.*1000.))**model_params['kappa']  
        # background_t_ref_trop_lian_sugiyama2014_alt = (160.)*(dataset['pfull']/(0.6*1000.))**model_params['kappa']  


        # plt.figure()
        # plt.plot(background_t_ref_trop_lian_showman_2008, dataset['pfull'], label='lian 2008')
        # plt.plot(background_t_ref_trop_lian_sugiyama2014, dataset['pfull'], label='sugiyama 2014')     
        # plt.plot(background_t_ref_trop_lian_sugiyama2014_alt, dataset['pfull'], label='sugiyama 2014 alt')        

        # plt.ylim([dataset.pfull.max(), 0.])    
        # plt.legend()

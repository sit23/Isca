import os

import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 16
NUM_LEVELS = 25

base_dir = os.path.dirname(os.path.realpath(__file__))
cb = SocratesCodeBase.from_directory(GFDL_BASE)
ALBEDO = 0.1
#exp = Experiment('soc_test_aquaplanet_with_clouds_linear_clr_albedo_'+str(ALBEDO)+'_qcl_with_temp', codebase=cb)
exp = Experiment('soc_test_aquaplanet_with_clouds_linear_clr_albedo_'+str(ALBEDO)+'_qcl_two_paras', codebase=cb)
exp.clear_rundir()

inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

diag.add_field('socrates', 'soc_tdt_lw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_rad', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True)
diag.add_field('socrates', 'soc_toa_swup', time_avg=True)

diag.add_field('socrates', 'soc_olr_clr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw_clr', time_avg=True)
diag.add_field('socrates', 'soc_toa_swup_clr', time_avg=True)
diag.add_field('socrates', 'soc_flux_lw_clr', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw_clr', time_avg=True)

diag.add_field('cloud_simple', 'cf', time_avg=True)
diag.add_field('cloud_simple', 'reff_rad', time_avg=True)
diag.add_field('cloud_simple', 'frac_liq', time_avg=True)
diag.add_field('cloud_simple', 'qcl_rad', time_avg=True)
diag.add_field('cloud_simple', 'simple_rhcrit', time_avg=True)
diag.add_field('cloud_simple', 'rh_in_cf', time_avg=True)


# additional output options commented out 
diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_co2', time_avg=True)
#diag.add_field('socrates', 'soc_ozone', time_avg=True)
#diag.add_field('socrates', 'soc_coszen', time_avg=True) 
#diag.add_field('socrates', 'soc_spectral_olr', time_avg=True)

exp.diag_table = diag
exp.inputfiles = inputfiles

coeff_fn = '/home/links/ql260/Documents/Exps_Analysis/clouds/data/paras_a_b.npy'
arr = np.load(coeff_fn)
arr_tmp = np.ones((100,2), dtype='double')*-999
arr_tmp[0:len(arr),:] = arr
coeff_arr = list(np.reshape(arr_tmp, np.size(arr_tmp), order='F')) # Fortran order

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':1370.,
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'dt_rad':3600,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
        'solday':90
    }, 
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,            
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection            
        'do_cloud_simple': True,
        'do_clear_sky_pass': True,
    },

    'cloud_simple_nml': {
        'simple_cca':0.0,
        'rhcsfc': 0.95,
        'rhc700': 0.7,
        'rhc200': 0.3,
        'do_simple_rhcrit': False,
        'do_read_scm_rhcrit': False,
        'do_linear_diag': True,
        'scm_linear_coeffs': coeff_arr, 
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,  
        'depth': 2.5,
        'albedo_value': ALBEDO,
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
           'do_simple':True,
           'construct_table_wrt_liq_and_ice':True
       },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,      
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels': NUM_LEVELS,      #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },

})

#Lets do a run!
if __name__=="__main__":
        cb.compile(debug=False)
        exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)#, run_idb=True)

        for i in range(2, 25):
            exp.run(i, num_cores=NCORES)

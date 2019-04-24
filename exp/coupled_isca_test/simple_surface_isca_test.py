import os

import numpy as np

from isca import IscaMOMCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 2
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaMOMCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('simple_surface_frierson_test_dry_no_add_phys_t_surf_rad_35', codebase=cb)

exp.inputfiles = [os.path.join(base_dir,'input/grid_spec.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_daily', 1, 'days', time_units='days')
# diag.add_file('ocean_daily', 1, 'days', time_units='days')
# diag.add_file('ice_daily', 1, 'days', time_units='days')


#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'bk', files = ['atmos_daily'])
diag.add_field('dynamics', 'pk', files = ['atmos_daily'])
diag.add_field('atmosphere', 'precipitation', time_avg=True, files = ['atmos_daily'])
diag.add_field('atmosphere', 't_surf_in', time_avg=True, files = ['atmos_daily'])

diag.add_field('simple_surface', 't_surf', time_avg=True, files = ['atmos_daily'])
diag.add_field('simple_surface', 'shflx', time_avg=True, files = ['atmos_daily'])
diag.add_field('simple_surface', 'evap', time_avg=True, files = ['atmos_daily'])
diag.add_field('simple_surface', 'lwflx', time_avg=True, files = ['atmos_daily'])
diag.add_field('simple_surface', 'swflx', time_avg=True, files = ['atmos_daily'])
diag.add_field('simple_surface', 'tdt_surf', time_avg=True, files = ['atmos_daily'])

diag.add_field('dynamics', 'sphum', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'ucomp', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'vcomp', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'temp', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'vor', time_avg=True, files = ['atmos_daily'])
diag.add_field('dynamics', 'div', time_avg=True, files = ['atmos_daily'])

diag.add_field('two_stream', 'net_lw_surf', time_avg=True, files = ['atmos_daily'])
diag.add_field('two_stream', 'swdn_sfc', time_avg=True, files = ['atmos_daily'])
diag.add_field('two_stream', 'swdn_toa', time_avg=True, files = ['atmos_daily'])
diag.add_field('two_stream', 'coszen', time_avg=True, files = ['atmos_daily'])
diag.add_field('two_stream', 'tdt_rad', time_avg=True, files = ['atmos_daily'])
diag.add_field('two_stream', 'tdt_solar', time_avg=True, files = ['atmos_daily'])



# diag.add_field('ocean_model', 'salt', time_avg=True, files = ['ocean_daily'])
# diag.add_field('ocean_model', 'temp', time_avg=True, files = ['ocean_daily'])
# diag.add_field('ocean_model', 'u', time_avg=True, files = ['ocean_daily'])
# diag.add_field('ocean_model', 'v', time_avg=True, files = ['ocean_daily'])

# diag.add_field('ice_model', 'CN', time_avg=True, files = ['ice_daily'])
# diag.add_field('ice_model', 'EXT', time_avg=True, files = ['ice_daily'])

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'coupler_nml':{
     'days'   : 5,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'dt_ocean':600,     
     'dt_cpld':600,          
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day',
     'do_land' : False,
     'do_ice' : False,
     'do_ocean': False,
     'do_flux': False,
     'print_s_messages': False,

    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':False,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': True,     #Use grey radiation
        'convection_scheme': 'None', #Use the simple Betts Miller convection scheme from Frierson
    },

    'moist_processes_nml': {
        'print_s_messages': False,
    },

    'physics_driver_nml': {
        'print_s_messages': False,
        'diffusion_smooth': False,
    },

    'atmos_model_nml': {
        'print_s_messages': False,
    },

    'atmosphere_nml': {
        'print_s_messages': False,
    },

    'simple_surface_nml': {
        'print_s_messages': False,
        'heat_capacity': 4.e07,
        'Tm': -275.,
        'evaporation' : False,
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

#     'atmosphere_nml': {
#         'idealized_moist_model': True
#     },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 275.,
        'delta_T': 0.,        
        'evaporation':False,   
        'depth': 2.5,                          #Depth of mixed layer used
        'albedo_value': 0.30,                  #Albedo value used          
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },

    'betts_miller_nml': {
       'rhbm': .7   , 
       'do_simp': False, 
       'do_shallower': True, 
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  5000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': False,                #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
        'atm_abs': 0.2,                      # default: 0.0        
        # 'ir_tau_eq': 0.,
        # 'ir_tau_pole': 0.,        
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
        'num_levels':25,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.],
        'vert_coord_option':'input', #Use the vertical levels from Frierson 2006
        'surf_res':0.5,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },
    'vert_coordinate_nml': {
        'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
        'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
       },

    'ocean_barotropic_nml':{
            'barotropic_time_stepping_B':True,
       },
    'ocean_tracer_nml': {

            'frazil_heating_after_vphysics':True,
       },
    'ocean_density_nml': {
            'eos_linear':True,
       },
})
#Verified as giving exactly the same result as Frierson test case.
#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    # exp.run(2, num_cores=NCORES, overwrite_data=False)

    # for i in range(2,25):
    #     exp.run(i, num_cores=NCORES)

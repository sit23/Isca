import os

import numpy as np
import pdb
from isca import IscaCodeBase, GreyCodeBase, ColumnCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from ntfy import notify
import datetime

# notify('Job running __'+str(os.path.basename(__file__))+'__ on Isca started at '+str(datetime.datetime.now().time()), 'Isca cpu update')

NCORES = 16

base_dir = os.path.dirname(os.path.realpath(__file__))

# a CodeBase can be a directory on the computer,
# useful for iterative development
# cb = ColumnCodeBase.from_directory(GFDL_BASE)
cb = GreyCodeBase.from_directory(GFDL_BASE)
#cb = GreyCodeBase.from_repo(repo='git@github.com:sit23/Isca.git', commit='c73439f')

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
exp = Experiment('dry_giant_planet_3d_60_levels_with_conv_with_diff_mk4_no_deep_velocity', codebase=cb)

exp.inputfiles = [os.path.join(base_dir,'input/jupiter_column_with_diff_input_file_4300.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_3monthly', 90, 'days', time_units='days')
# diag.add_file('atmos_daily', 1, 'days', time_units='days')
# diag.add_file('atmos_10_mins', 600, 'seconds', time_units='days')


#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
#diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'dt_qg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_convection_moist', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_convection_dry', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
#diag.add_field('atmosphere', 'dt_qg_condensation', time_avg=True)
#diag.add_field('atmosphere', 'dt_tg_condensation', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'snow', time_avg=True)
diag.add_field('atmosphere', 'klzb', time_avg=True)
diag.add_field('atmosphere', 'plzb', time_avg=True)
diag.add_field('atmosphere', 'plcl', time_avg=True)
#diag.add_field('atmosphere', 't_ref', time_avg=True)
#diag.add_field('atmosphere', 'q_ref', time_avg=True)
#diag.add_field('atmosphere', 'convflag', time_avg=True)
#diag.add_field('atmosphere', 'dtemp_bb_surface_flux', time_avg=True)
#diag.add_field('atmosphere', 'dsphum_bb_surface_flux', time_avg=True)
diag.add_field('atmosphere', 'cin', time_avg=True)
diag.add_field('atmosphere', 'cape', time_avg=True)
diag.add_field('atmosphere', 'cin_dry', time_avg=True)
# diag.add_field('atmosphere', 'cape_dry', time_avg=True)
diag.add_field('two_stream', 'tdt_rad', time_avg=True)
diag.add_field('two_stream', 'tdt_solar', time_avg=True)
diag.add_field('vert_turb', 'z_pbl', time_avg=True)
diag.add_field('vert_turb', 'diff_t', time_avg=True)
diag.add_field('vert_turb', 'diff_m', time_avg=True)
diag.add_field('rayleigh_bottom_drag', 'udt_rd', time_avg=True)
diag.add_field('rayleigh_bottom_drag', 'vdt_rd', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#s Namelist changes from default values
exp.namelist = namelist = Namelist({

    'main_nml': {
     'days'   : 90,	
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,			
     'dt_atmos':1800,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'no_calendar'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :True,
        'do_simple': False,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,     
        'two_stream_gray':  True, #Use grey radiation
        'do_rrtm_radiation':  False, #Don't use RRTM
        'convection_scheme':  'dry', #Use the dry convection scheme of Schneider & Walker
        'gp_surface':  True, #Use the giant-planet option for the surface, meaning it's not a mixed layer, and applies energy conservation and a bottom-boundary heat flux
        'mixed_layer_bc':  False, #Don't use the mixed-layer surface
        'remove_negative_sphum':False,  
        'do_lscale_cond':False,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': False,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': False,
        # Following three options are new and should match the vertical diffusion included in Peter and Roland's new 2018 paper's appendix. Needs testing though - does it make any difference?? Tried in column model but highly unstable. Try somewhere with velocities...
        'free_atm_diff':True,
        'young_read_jupiter_diffusivity':True,
        'young_read_timescale':21600.0,
        'rich_crit_diff':3.0,
        'do_pbl_diff':False,
        'rich_prandtl':0.,
        # 'fixed_depth':True,
        # 'depth_0': 0.0
    },

    'surface_flux_nml': {
        'use_virtual_temp': True,
        'do_simple': False,
        'old_dtaudv': True, 
        'diabatic_acce':  1.0, #Parameter to artificially accelerate the diabatic processes during spinup. 1.0 performs no such acceleration, >1.0 performs acceleration        
        'gp_deep_water_source': False,
        'q_deep_sphum': 0.0,
        # 'flux_heat_gp':0.0,
    },
#         'q_deep_sphum': 6.1e-3,

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,        
    },
    
    'lscale_cond_nml': {
        'do_simple':False,
        'do_evap':False
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':False,
        'tcmin':  -223, #Make sure low temperature limit of saturation vapour pressure is low enough that it doesn't cause an error (note that this giant planet has no moisture anyway, so doesn't directly affect calculation.        
        'tcmax': 350.,
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  16000.,
        'do_conserve_energy': True,         
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'Schneider', #Use the Schneider & Liu option for the grey scheme
        'do_seasonal': False,               #Don't use seasonally-varying radiation 
        'atm_abs': 0.2,                      # default: 0.0
        'solar_constant':  50.7, #Change solar constant
        'diabatic_acce':  1.0, #Parameter to artificially accelerate the diabatic processes during spinup. 1.0 performs no such acceleration, >1.0 performs acceleration        
        'do_normal_integration_method':False,  

    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
         'domains_stack_size': 5000000 #Setting size of stack available to model, which needs to be higher than the default when running at high spatial resolution
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'water_correction_limit': 200.e2,
        'valid_range_t':[50.,800.],
        'vert_coord_option':'uneven_sigma', #Use equally-spaced sigma levels
        'surf_res':1.0, #Parameter for setting vertical distribution of sigma levels
        'scale_heights' : 5.0, #Number of vertical scale-heights to include
        'exponent':7.0,#Parameter for setting vertical distribution of sigma levels
        'robert_coeff':0.03,
        'num_fourier':  85, #Number of Fourier modes
        'num_spherical':  86, #Number of spherical harmonics in triangular truncation
        'lon_max':  256, #Lon grid points
        'lat_max':  128, #Lat grid points
        'num_levels':  60, #Number of vertical levels
        'do_water_correction':  False, #Turn off enforced water conservation as model is dry
        'damping_option':  'exponential_cutoff', #Use the high-wavenumber filter option
        'damping_order':  4,
        'damping_coeff':  1.3889e-04,
        'cutoff_wn':  40,
        'initial_sphum':  0.0909, #No initial specific humidity   
        'reference_sea_level_press':  15.0e5,
        'initial_state_option':'input',
        'use_virtual_temperature':True,
    },


    # 'column_nml': {
    #     'lon_max': 2, # number of columns in longitude, default begins at lon=0.0
    #     'lat_max': 128, # number of columns in latitude, precise 
    #                   # latitude can be set in column_grid_nml if only 1 lat used. 
    #     'num_levels': 60,  # number of levels 
    #     'initial_sphum': 0.0, 
    #     'valid_range_t':[10.,700.], 
    #     'reference_sea_level_press':15.0e5,
    #     'vert_coord_option':'uneven_sigma', #Use equally-spaced sigma levels
    #     'surf_res':1.0, #Parameter for setting vertical distribution of sigma levels
    #     'scale_heights' : 5.0, #Number of vertical scale-heights to include
    #     'exponent':7.0,#Parameter for setting vertical distribution of sigma levels        
    # },

    # set initial condition, NOTE: currently there is not an option to read in initial condition from a file. 
    # 'column_init_cond_nml': {
    #     'initial_temperature': 264., # initial atmospheric temperature 
    #     'surf_geopotential': 0.0, # applied to all columns 
    #     'surface_wind': 5. # as described above 
    # },

    'ic_from_external_file_nml': {
        'file_name':'INPUT/jupiter_column_with_diff_input_file_4300.nc',    
    },    

    'constants_nml': {
#Set Jupiter constants
        'radius':  69860.0e3,
        'grav':  26.0,
        'omega':  1.7587e-4,
        'orbital_period':  4332.589*86400.,
        'PSTD':  15.0E+06,
        'PSTD_MKS':  15.0E+05,
        'rdgas':  3605.38,    
    },

#Set parameters for near-surface Rayleigh drag
    'rayleigh_bottom_drag_nml':{
    	'kf_days':10.0,
            'do_drag_at_surface':True,
            'variable_drag': False,
            'zero_eq_drag':False,
            'do_energy_conserv_ray':True,
	},

#Set parameters for dry convection scheme
    'qe_moist_convection_nml': {
        'tau_bm': 21600.,
        'rhbm':0.75,
        'k_surface_set': 25,
        'tmin': 150.,
        'tmax':620
    },

    'dry_convection_nml': {
        'tau' : 21600., 
        'gamma' : 1.0,
        'use_inertial_timescale':False,
#         'max_height_level_set':24,
    },

      
})


#exp.run(1441, multi_node=False, use_restart=True, restart_file = base_dir+'/input/moist_no_drag_sbm_res1440.tar.gz', num_cores=NCORES)
exp.run(1, multi_node=False, use_restart=False, num_cores=NCORES)
for i in range(2,5001):
   exp.run(i, num_cores=NCORES, multi_node=False)

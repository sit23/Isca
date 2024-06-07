import os

import numpy as np

from isca import ShallowCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import pdb

NCORES = 32
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = ShallowCodeBase.from_directory(GFDL_BASE)

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

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_50days', 50, 'days', time_units='days')
# diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('shallow_diagnostics', 'ucomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vcomp', time_avg=True)
diag.add_field('shallow_diagnostics', 'vor', time_avg=True)
diag.add_field('shallow_diagnostics', 'div', time_avg=True)
diag.add_field('shallow_diagnostics', 'h', time_avg=True)
diag.add_field('shallow_diagnostics', 'pv', time_avg=True)
diag.add_field('shallow_diagnostics', 'pv_corrected', time_avg=True)

diag.add_field('shallow_diagnostics', 'stream', time_avg=True)
diag.add_field('shallow_diagnostics', 'deep_geopot', time_avg=True)
# diag.add_field('shallow_diagnostics', 'trs', time_avg=True)
# diag.add_field('shallow_diagnostics', 'tr', time_avg=True)
diag.add_field('stirring_mod', 'stirring', time_avg=True)
# diag.add_field('stirring_mod', 'stirring_amp', time_avg=True)
# diag.add_field('stirring_mod', 'stirring_sq', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_kin', time_avg=True)
diag.add_field('shallow_diagnostics', 'kinetic_energy_density', time_avg=True)
diag.add_field('shallow_diagnostics', 'h_sqd_mean', time_avg=True)
diag.add_field('physics_mod', 'dt_hg_physical_forcing', time_avg=True)
diag.add_field('physics_mod', 'dt_hg_rad_forcing', time_avg=True)

diag_for_movie = DiagTable()
diag_for_movie.add_file('atmos_daily', 1, 'days', time_units='days')
diag_for_movie.add_field('shallow_diagnostics', 'vor', time_avg=True)
diag_for_movie.add_field('shallow_diagnostics', 'h', time_avg=True)
diag_for_movie.add_field('shallow_diagnostics', 'ucomp', time_avg=True)
# diag_for_movie.add_field('stirring_mod', 'stirring', time_avg=True)

#Empty the run directory ready to run

#Define values for the 'core' namelist
namelist = Namelist({
    'main_nml':{
     'days'   : 120,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos': 400,
     'calendar': 'no_calendar',
    },

 'atmosphere_nml':{
   'print_interval': 86400,
    },

'fms_io_nml':{
   'threading_write' :'single',
   'fileset_write': 'single'
    },

 'fms_nml':{
   'print_memory_usage':False,
   'domains_stack_size': 200000,
    },

 'shallow_dynamics_nml':{
   'num_lon'             : 2048,
   'num_lat'             : 1024,
   'num_fourier'         : 682,
   'num_spherical'       : 683,
   'fourier_inc'         : 1,
   'damping_option'      : 'resolution_dependent',
   'damping_order'       : 4,
   'damping_coeff'       : 1.e-04,
   'h_0'                 : 3.e04,
   'grid_tracer'         : False,
   'spec_tracer'         : False,
   'robert_coeff'        : 0.04,
   'robert_coeff_tracer' : 0.04,
    },

 'shallow_physics_nml':{
   'fric_damp_time'  :  0.0,
   'h_0'             :  3.e04, #Scott and Polvani begin with L_D_polar = 10 radii in size, i.e. depth of 9,910 metres.
   'h_amp'           :  0.,
   'h_itcz'          :  0.,
   },

 'stirring_nml': {
   'B':0.0,
   },

})

#Let's do a run!
if __name__=="__main__":


    exp = Experiment(f'showman_2007_A1_mk3_single', codebase=cb)

    exp.diag_table = diag 
    exp.namelist = namelist 
    exp.clear_rundir()             

    gh = 6e4
    taurad = 0.
    omega = 1.74e-4
    radius = 7.14e7
    dt = 400.
    total_time = 100


    exp.update_namelist({
        'shallow_dynamics_nml':{
            'h_0': gh,
            'u_deep_mag'   : 0.,
            # 'damping_coeff': eta
        },
        'shallow_physics_nml': {
            'h_0': gh, 
            'therm_damp_time' : 0.5*taurad, #Thermal damping time in Scott and Polvani is 1 rotation period (v_l = 1)
            'fric_damp_time': 0., 
            'h_width_storm': 0.7,
            'storm_strength_param': 0.333,
            'storm_interval': 1.25e4, 
            'storm_length' : 1e5,
        },
        'constants_nml': {
            'omega': omega,
            'radius': radius,
        },
        'stirring_nml': {
            # 'decay_time':10.*rotation_period, #Scott and Polvani use decorrelation time of 10 planetary rotations - i.e. a long time for Jupiter. 
            'amplitude':0.,
            # 'n_total_forcing_max': 85,
            # 'n_total_forcing_min': 79,
            # 'delta_stirring_lat' : delta_stirring_lat, 
            # 'do_ampy_localize': True,
            # 'do_localize': True,
        },
        'main_nml':{
            'dt_atmos':dt,   
            'days': total_time
        },
    })

    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2,51):
        exp.run(i, num_cores=NCORES)

    # exp.diag_table = diag_for_movie
    # exp.run(21, num_cores=NCORES)

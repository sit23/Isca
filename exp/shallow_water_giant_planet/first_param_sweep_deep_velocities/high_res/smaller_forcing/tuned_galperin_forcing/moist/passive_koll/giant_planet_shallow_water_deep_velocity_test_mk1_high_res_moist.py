import os

# import numpy as np

from isca import ShallowCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

import pdb

NCORES = 2
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
# diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily', 1, 'days', time_units='days')

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
diag.add_field('shallow_diagnostics', 'tr', time_avg=True)
diag.add_field('shallow_diagnostics', 'evap', time_avg=True)
diag.add_field('shallow_diagnostics', 'precip', time_avg=True)
diag.add_field('shallow_diagnostics', 'rh', time_avg=True)


diag.add_field('stirring_mod', 'stirring', time_avg=True)
# diag.add_field('stirring_mod', 'stirring_amp', time_avg=True)
diag.add_field('stirring_mod', 'stirring_sqr', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_kin', time_avg=True)
diag.add_field('shallow_diagnostics', 'h_sqd_mean', time_avg=True)

diag.add_field('shallow_diagnostics', 'e_kin_density', time_avg=True)
diag.add_field('shallow_diagnostics', 'eq_geopot', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_kin_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_pot_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'e_tot_real_units', time_avg=True)
diag.add_field('shallow_diagnostics', 'u_rms', time_avg=True)

#Empty the run directory ready to run

#Define values for the 'core' namelist
namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos': 300,
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
   'num_lon'             : 128,
   'num_lat'             : 64,
   'num_fourier'         : 42,
   'num_spherical'       : 43,
   'fourier_inc'         : 1,
   'damping_option'      : 'resolution_dependent',
   'damping_order'       : 4,
   'damping_coeff'       : 1.e-04,
   'h_0'                 : 2e5,
   'grid_tracer'         : True,
   'spec_tracer'         : False,
   'robert_coeff'        : 0.04,
   'robert_coeff_tracer' : 0.04,
   'sat_range_initial'   : 0.,
   'small_gridpoint_amp': 2.e-5,
    },

 'shallow_physics_nml':{
   'fric_damp_time'  :  0.0,
   'h_0'             :  2e5, #Scott and Polvani begin with L_D_polar = 10 radii in size, i.e. depth of 9,910 metres.
   'h_amp'           :  0.,
   'h_itcz'          :  0.,
   'sat_constant'    : 1.0e-5,
   'precip_timescale' : 1.0e7,
   'evap_prefactor'  : 0.0,
   'latent_heat_prefactor' : 1.0,
   },

 'stirring_nml': {
   'B':0.0,
   'do_localize': False,
   },

})

#Lets do a run!
if __name__=="__main__":

    for precip_timescale in [1e3, ]:

        for damping_time in [ 0.]:
            u_deep_mag_val = 0.

            if u_deep_mag_val!=0.:
                u_deep_merid_arr = [27]
            else:
                u_deep_merid_arr = [3]

            for u_deep_merid in u_deep_merid_arr:

                ld_value = 0.025
                exp = Experiment('passive_tidal_lock_0.5_delta_h_with_sp_drag_fix_gh_mk2', codebase=cb)

                exp.diag_table = diag 
                exp.namelist = namelist 
                exp.clear_rundir()

                rotation_period = 2.3*86400.
                omega = 2.*3.1415/ rotation_period
                radius = 8.2e7
                grav = 20.0
                number_ld_in_radius_units = ld_value

                # Model uses geopotential as its height co-ordinate. So depth is h_0/grav.

                equilibrium_geopotential = (2.*number_ld_in_radius_units*omega*radius)**2.
                equilibrium_depth = equilibrium_geopotential/grav

                exp.update_namelist({
                    'shallow_dynamics_nml':{
                        'h_0': 2e5*20,
                        'u_deep_mag'   : u_deep_mag_val,
                        'n_merid_deep_flow': u_deep_merid,  
                        'valid_range_v': [-1e4, 1e4]
                    },
                    'shallow_physics_nml': {
                        'h_0': 2e5*20, 
                        'therm_damp_time' : 0.1*86400.,
                        'fric_damp_time': 10.0*86400.,
                        'precip_timescale': precip_timescale, 
                        'latent_heat_prefactor': 0.,
                        'h_amp_tidal_locked': 0.5*2e5*20,
                        'showman_polvani_drag':True,
                    },
                    'constants_nml': {
                        'omega': omega,
                        'radius': radius,
                    },

                })

                exp.run(1, use_restart=False, num_cores=NCORES, multi_node=False)
                for i in range(2,25):
                    exp.run(i, num_cores=NCORES, multi_node=False)

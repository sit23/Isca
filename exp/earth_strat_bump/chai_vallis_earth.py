import numpy as np
import os
from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


NCORES = 16
RESOLUTION = 'T42', 25  # T42 horizontal resolution, 25 levels in pressure

base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

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
diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vcomp_temp', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('hs_forcing', 'teq', time_avg=True)


# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 600,
        'days': 30,
        'calendar': 'thirty_day',
        'current_date': [1,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    }
})



#Lets do a run!
if __name__ == '__main__':

    for forcing_file in ['earth_chai_vallis_teq', 'earth_chai_vallis_teq_final_cs']:
        exp_name = f'{forcing_file}_t42'

        exp = Experiment(exp_name, codebase=cb)
        exp.clear_rundir()
        exp.diag_table = diag
        exp.namelist = namelist
        exp.set_resolution(*RESOLUTION)
        exp.inputfiles = [os.path.join(base_dir,f'input/{forcing_file}.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]        

        exp.update_namelist({
            'hs_forcing_nml': { 
                'equilibrium_t_option': 'from_file',
                'equilibrium_t_file':f'{forcing_file}'
            }
        }) 

        exp.run(1, num_cores=NCORES, use_restart=False, overwrite_data=True)
        for i in range(2, 121):
            exp.run(i, num_cores=NCORES, overwrite_data=True)  # use the restart i-1 by default
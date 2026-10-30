"""Check run for the present-day tuning: the SimCloud atmosphere of realistic_simcloud_mom6.py
(same physics, forcing and TUNE settings) over prescribed climatological SST and sea ice
(exp/test_cases/simple_clouds/input: sst_clim_amip, siconc_clim_amip), no ocean model.
Groups that tune coupled models check that such a run has a net TOA flux of about
+0.5 to +1 W/m2 (Schmidt et al. 2017); it is a check, not a tuning target.

Usage: EXP_NAME=SA1 NMONTHS=60 NCORES=16 python realistic_simcloud_amip.py
"""
import os

from isca import Experiment

import realistic_simcloud_mom6 as cpl

exp = Experiment(os.environ.get('EXP_NAME', 'realistic_simcloud_amip'), codebase=cpl.cb)
exp.diag_table = cpl.diag
exp.inputfiles = [cpl.land_file, os.path.join(cpl.GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc'),
                  os.path.join(cpl.sc_dir, 'input/sst_clim_amip.nc'),
                  os.path.join(cpl.sc_dir, 'input/siconc_clim_amip.nc')]
exp.namelist = cpl.exp.namelist.copy()
exp.namelist['idealized_moist_phys_nml']['surface_tiles'] = False
exp.namelist['mixed_layer_nml'].update({
    'do_read_sst': True, 'do_sc_sst': True, 'sst_file': 'sst_clim_amip',
    'specify_sst_over_ocean_only': True, 'update_albedo_from_ice': True,
    'update_land_mask_from_ice': True, 'ice_albedo_value': 0.7, 'ice_albedo_method': 'ramp_function'})

if __name__ == '__main__':
    for i in range(int(os.environ.get('FIRST_MONTH', 1)), cpl.NMONTHS + 1):
        exp.run(i, use_restart=(i > 1), num_cores=cpl.NCORES, overwrite_data=True,
                mpirun_opts=os.environ.get('MPIRUN_OPTS', '--bind-to none'))

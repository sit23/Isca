"""Isca with realistic continents (T42, RRTM) coupled to the Samudra ML ocean emulator.

Based on exp/test_cases/realistic_continents, but the q-flux slab ocean is
replaced by Samudra (1-degree, 5-day steps) through OASIS3-MCT. Points under
the prescribed AMIP sea-ice climatology keep Isca's slab with ice albedo, and
the ocean beneath them receives no fluxes.

Needs an OASIS-enabled Isca build and a Python environment with Samudra
(SAMUDRA_PYTHON). Usage:
    NMONTHS=12 python realistic_samudra.py
"""
import os

import f90nml
import netCDF4
import numpy as np

from isca import IscaCodeBase, DiagTable, GFDL_BASE
from isca.coupling import CoupledExperiment, SamudraOcean

NCORES = int(os.environ.get('NCORES', 16))
NMONTHS = int(os.environ.get('NMONTHS', 1))
DAYS = int(os.environ.get('DAYS', 30))
EXP_NAME = os.environ.get('EXP_NAME', 'realistic_samudra')
SAMUDRA_PYTHON = os.environ.get('SAMUDRA_PYTHON', '/disca/share/sit204/venvs/samudra/bin/python')
SAMUDRA_CONFIG = os.environ.get('SAMUDRA_CONFIG', '/disca/share/sit204/samudra/m1_eval/eval.yaml')
SAMUDRA_CKPT = os.environ.get('SAMUDRA_CKPT',
                              '/disca/share/sit204/samudra/v2_1deg_20261008/ema_ckpt.pt')

rc_dir = os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents')
land_file = os.path.join(GFDL_BASE, 'input/land_masks/era_land_t42.nc')
with netCDF4.Dataset(land_file) as nc:
    land = np.asarray(nc['land_mask'][:]) > 0.
    zsurf = np.asarray(nc['zsurf'][:])

cb = IscaCodeBase.from_directory(GFDL_BASE)
ocean = SamudraOcean(SAMUDRA_PYTHON, SAMUDRA_CONFIG, SAMUDRA_CKPT, threads=8)
exp = CoupledExperiment(EXP_NAME, codebase=cb, ocean=ocean, coupling_period=5 * 86400,
                        land_mask=land, zsurf=zsurf)
exp.inputfiles = [land_file, os.path.join(GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc'),
                  os.path.join(rc_dir, 'input/siconc_clim_amip.nc')]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
for mod, field in [('dynamics', 'ps'), ('dynamics', 'bk'), ('dynamics', 'pk'), ('dynamics', 'zsurf'),
                   ('atmosphere', 'precipitation'), ('mixed_layer', 't_surf'),
                   ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
                   ('mixed_layer', 'coupled_ocean_frac'), ('mixed_layer', 'cpl_net_heat'),
                   ('mixed_layer', 'ice_conc'), ('atmosphere', 'cpl_runoff'),
                   ('dynamics', 'ucomp'), ('dynamics', 'vcomp'), ('dynamics', 'temp'),
                   ('dynamics', 'sphum'), ('rrtm_radiation', 'olr'), ('rrtm_radiation', 'toa_sw')]:
    diag.add_field(mod, field, time_avg=True)
exp.diag_table = diag

exp.namelist = f90nml.read(os.path.join(rc_dir, 'namelist_basefile.nml'))
exp.namelist['main_nml']['days'] = DAYS
# The ocean model replaces the q-flux; ice points keep the slab with ice albedo.
exp.namelist['mixed_layer_nml'].update({'do_qflux': False, 'load_qflux': False})

if __name__ == '__main__':
    cb.compile()
    first = int(os.environ.get('FIRST_MONTH', 1))
    for i in range(first, NMONTHS + 1):
        exp.run(i, use_restart=(i > 1), num_cores=NCORES, overwrite_data=True,
                mpirun_opts='--bind-to none')

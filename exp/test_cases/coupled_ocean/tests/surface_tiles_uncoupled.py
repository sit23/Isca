"""Fractional surface tiles in an uncoupled run (the bucket_hydrology test case: RRTM,
continents from its land file, land bucket).

  CASE=off     no tiles
  CASE=binary  surface_tiles with the same 0/1 land: should reproduce 'off' bit for bit
  CASE=half    surface_tiles with every land cell half land, half ocean

Usage: CASE=binary DAYS=10 python surface_tiles_uncoupled.py
"""
import os
import sys

import netCDF4
import numpy as np

from isca import GFDL_BASE, Experiment

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp', 'test_cases', 'bucket_hydrology'))
import bucket_model_test_case as bm   # noqa: E402

CASE = os.environ.get('CASE', 'binary')
DAYS = int(os.environ.get('DAYS', 10))
NCORES = int(os.environ.get('NCORES', 4))

exp = Experiment('tiles_uncoupled_%s' % CASE, codebase=bm.cb)
exp.namelist = bm.namelist.copy()
exp.diag_table = bm.diag
if os.environ.get('HOURLY'):
    # hourly snapshots of the surface exchange, to find where tiles first differ
    exp.diag_table.add_file('atmos_hourly', 1, 'hours', time_units='hours')
    for mod, field in [('mixed_layer', 't_surf'), ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
                       ('atmosphere', 'bucket_depth'), ('atmosphere', 'precipitation'),
                       ('dynamics', 'temp'), ('dynamics', 'sphum'), ('dynamics', 'ucomp')]:
        exp.diag_table.add_field(mod, field, time_avg=False, files=['atmos_hourly'])
exp.namelist['main_nml']['days'] = DAYS
if 'SECONDS' in os.environ:   # a few time steps, to find where two cases first differ
    exp.namelist['main_nml']['days'] = 0
    exp.namelist['main_nml']['seconds'] = int(os.environ['SECONDS'])
land_src = os.path.join(GFDL_BASE, 'exp', 'test_cases', 'bucket_hydrology', 'input', 'land.nc')
inputs = [f for f in bm.exp.inputfiles if not f.endswith('land.nc')]
if CASE == 'off':
    inputs.append(land_src)
else:
    exp.namelist['idealized_moist_phys_nml']['surface_tiles'] = True
    land = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'land_tiles_%s.nc' % CASE)
    with netCDF4.Dataset(land_src) as src, netCDF4.Dataset(land, 'w') as dst:
        for name, dim in src.dimensions.items():
            dst.createDimension(name, len(dim))
        for name, var in src.variables.items():
            v = dst.createVariable(name, var.dtype, var.dimensions)
            v.setncatts({a: var.getncattr(a) for a in var.ncattrs()})
            data = var[:]
            if name == 'land_mask' and CASE == 'half':
                data = 0.5 * (np.asarray(data) > 0)
            v[:] = data
    inputs.append(land)
exp.inputfiles = inputs
if CASE != 'off':
    # the model reads INPUT/<basename>; point it at the generated file's name
    exp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/' + os.path.basename(land)
    exp.namelist['mixed_layer_nml']['land_option'] = 'input'

if __name__ == '__main__':
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)

"""Start a THERMO_ICE=1 run from the restart of a run made without ice.

Copies SRC_EXP's restart archive for month N to DST_EXP's restarts directory,
adding the fields the ice code reads: h_thermo_ice = 0, t_ml = t_surf and
albedo = ALBEDO in mixed_layer.res.nc, and FRAZIL_OCN = 0 in the OASIS restart
of the fields the ocean sends. Then run DST_EXP with FIRST_MONTH=N+1.

Usage: python seed_ice_restart.py SRC_EXP N DST_EXP [ALBEDO]
"""
import os
import shutil
import subprocess
import sys
import tempfile

import netCDF4

D = os.environ['GFDL_DATA']
src, n, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
albedo = float(sys.argv[4]) if len(sys.argv) > 4 else 0.31
archive = 'res%04d.tar.gz' % n

tmp = tempfile.mkdtemp()
subprocess.check_call(['tar', 'xzf', os.path.join(D, src, 'restarts', archive)], cwd=tmp)

with netCDF4.Dataset(os.path.join(tmp, 'mixed_layer.res.nc'), 'a') as nc:
    ts = nc['t_surf']
    for name, value in [('h_thermo_ice', 0.), ('t_ml', ts[:]), ('albedo', albedo)]:
        if name not in nc.variables:
            v = nc.createVariable(name, 'f8', ts.dimensions)
            v.long_name, v.units = name, 'none'
            v[:] = value

with netCDF4.Dataset(os.path.join(tmp, 'ocn_rst.nc'), 'a') as nc:
    if 'FRAZIL_OCN' not in nc.variables:
        ny, nx = nc['SST_OCN'].shape
        nc.createDimension('FRAZIL_OCN_nx', nx)
        nc.createDimension('FRAZIL_OCN_ny', ny)
        nc.createVariable('FRAZIL_OCN', 'f8', ('ny', 'nx'))[:] = 0.

outdir = os.path.join(D, dst, 'restarts')
os.makedirs(outdir, exist_ok=True)
subprocess.check_call(['tar', 'czf', os.path.join(outdir, archive), '.'], cwd=tmp)
shutil.rmtree(tmp)
print('wrote', os.path.join(outdir, archive), '; run with FIRST_MONTH=%d' % (n + 1))

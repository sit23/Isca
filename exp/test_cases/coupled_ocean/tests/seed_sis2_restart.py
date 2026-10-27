"""Start a SIS2=1 run from the restart of a MOM6 run without SIS2.

Copies SRC_EXP's restart archive for month N to DST_EXP's restarts directory, adding the
OASIS restart fields the SIS2 coupling needs: ICE_FRAC = 0, ICE_TSRF = SST and
ICE_ALB = 0 in the fields the ocean sends, and zero sea-ice fluxes (ISCA_I*) in the
fields Isca sends. SIS2 itself starts with no ice (there is no ice restart). Then run
DST_EXP with SIS2=1 FIRST_MONTH=N+1.

Usage: python seed_sis2_restart.py SRC_EXP N DST_EXP
"""
import os
import shutil
import subprocess
import sys
import tempfile

import netCDF4

from isca.coupling.mom6 import ISCA_TO_SIS2

D = os.environ['GFDL_DATA']
src, n, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
archive = 'res%04d.tar.gz' % n

tmp = tempfile.mkdtemp()
subprocess.check_call(['tar', 'xzf', os.path.join(D, src, 'restarts', archive)], cwd=tmp)


def add(nc, name, values):
    if name in nc.variables:
        return
    ny, nx = values.shape
    for dim, size in (('%s_nx' % name, nx), ('%s_ny' % name, ny)):
        if dim not in nc.dimensions:
            nc.createDimension(dim, size)
    nc.createVariable(name, 'f8', ('ny', 'nx'))[:] = values


with netCDF4.Dataset(os.path.join(tmp, 'ocn_rst.nc'), 'a') as nc:
    sst = nc['SST_OCN'][:]
    add(nc, 'ICE_FRAC', 0. * sst)
    add(nc, 'ICE_TSRF', sst)
    add(nc, 'ICE_ALB', 0. * sst)

with netCDF4.Dataset(os.path.join(tmp, 'isca_rst.nc'), 'a') as nc:
    shape = nc['ISCA_TAUX'][:].shape
    for isca_name, _ in ISCA_TO_SIS2:
        if isca_name not in nc.variables:
            nc.createVariable(isca_name, 'f8', ('ny', 'nx'))[:] = 0.

outdir = os.path.join(D, dst, 'restarts')
os.makedirs(outdir, exist_ok=True)
subprocess.check_call(['tar', 'czf', os.path.join(outdir, archive), '.'], cwd=tmp)
shutil.rmtree(tmp)
print('wrote', os.path.join(outdir, archive), '; run with SIS2=1 FIRST_MONTH=%d' % (n + 1))

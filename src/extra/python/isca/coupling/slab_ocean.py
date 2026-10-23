#!/usr/bin/env python3
"""A slab ocean that runs as a separate OASIS3-MCT component.

This is the simplest possible ocean for Isca's OASIS coupling. It exists to
test the coupling plumbing (grids, partitions, signs, units, time averaging,
restarts) against Isca's own built-in slab ocean in mixed_layer.F90, before a
real ocean model such as MOM6 is involved.

It lives on the Isca grid itself (same grid name in the namcouple, so no
remapping), runs on one MPI process and steps

    rho_cp * depth * dSST/dt = SW_net + LW_net - SH - L * E

explicitly over each coupling period, using the fluxes averaged by OASIS.

Its only state is the SST, and that is kept in the OASIS coupling restart
file for SST_OCN, which OASIS reads at the start of a run and rewrites at the
end. On a cold start (no restart file) the file is created here, before
oasis_enddef, from Isca's prescribe_initial_dist formula.

Launched by isca.coupling.CoupledExperiment as the second program of an
MPMD run, e.g.
    mpirun -np 16 isca.x : -np 1 python slab_ocean.py --nlon 128 --nlat 64 ...
"""

import argparse
import os

import netCDF4
import numpy as np
from mpi4py import MPI
import pyoasis
from pyoasis import OASIS

RHO_CP = 1.035e3 * 3989.24495292815   # Isca constants_mod RHO0*CP_OCEAN
HLV = 2.500e6

FLUX_FIELDS = ['OCN_TAUX', 'OCN_TAUY', 'OCN_SWNET', 'OCN_LWNET', 'OCN_SENS',
               'OCN_EVAP', 'OCN_PRECIP', 'OCN_PSURF']


def gaussian_latitudes(nlat):
    """Gaussian latitudes [degrees], south to north, as used by Isca."""
    sin_lat, _ = np.polynomial.legendre.leggauss(nlat)
    return np.degrees(np.arcsin(sin_lat))


def write_cold_start_restart(fname, nlon, nlat, sst, ofrac):
    """Write an OASIS coupling restart file holding the initial SST and ocean fraction.

    OASIS restart fields are (nx, ny) in Fortran order, i.e. (ny, nx) as
    seen from netCDF4-python.
    """
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nx', nlon)
        nc.createDimension('ny', nlat)
        for name, data in (('SST_OCN', sst), ('FRAC_OCN', ofrac)):
            v = nc.createVariable(name, 'f8', ('ny', 'nx'))
            v[:] = data.reshape(nlat, nlon)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--nlon', type=int, required=True)
    p.add_argument('--nlat', type=int, required=True)
    p.add_argument('--runtime', type=int, required=True, help='run length [s]')
    p.add_argument('--coupling-period', type=int, required=True, help='[s]')
    p.add_argument('--depth', type=float, default=2.5, help='slab depth [m]')
    p.add_argument('--restart-file', default='ocn_rst.nc')
    p.add_argument('--tconst', type=float, default=285.,
                   help='cold start SST at the equator-ish, as Isca tconst')
    p.add_argument('--delta-t', type=float, default=40.,
                   help='cold start equator-pole contrast, as Isca delta_T')
    p.add_argument('--max-lat', type=float, default=90.,
                   help='ocean covers |lat| <= max-lat; poleward is left to Isca (tests masking)')
    args = p.parse_args()

    comm = MPI.COMM_WORLD
    comp = pyoasis.Component('slab_ocean', True, comm)
    if comp.localcomm.size != 1:
        raise SystemExit('slab_ocean.py runs on exactly one MPI process')

    nlon, nlat = args.nlon, args.nlat
    npts = nlon * nlat
    lat = np.broadcast_to(gaussian_latitudes(nlat)[:, None], (nlat, nlon)).ravel()
    ofrac = np.where(np.abs(lat) <= args.max_lat, 1.0, 0.0)

    if not os.path.exists(args.restart_file):
        sin_lat = np.sin(np.radians(lat))
        sst0 = args.tconst - args.delta_t * (3. * sin_lat**2 - 1.) / 3.
        write_cold_start_restart(args.restart_file, nlon, nlat, sst0, ofrac)
        print('slab_ocean: cold start, wrote %s' % args.restart_file, flush=True)

    # The SST state is whatever is in the restart file OASIS will send at t=0.
    with netCDF4.Dataset(args.restart_file) as nc:
        sst = np.array(nc['SST_OCN'][:], dtype=np.float64).ravel()

    partition = pyoasis.SerialPartition(npts)
    flux_vars = {name: pyoasis.Var(name, partition, OASIS.IN) for name in FLUX_FIELDS}
    sst_var = pyoasis.Var('SST_OCN', partition, OASIS.OUT)
    frac_var = pyoasis.Var('FRAC_OCN', partition, OASIS.OUT)
    comp.enddef()

    heat_capacity = RHO_CP * args.depth
    fields = {name: pyoasis.asarray(np.zeros(npts)) for name in FLUX_FIELDS}
    dt = args.coupling_period
    for date in range(0, args.runtime, dt):
        for name, var in flux_vars.items():
            var.get(date, fields[name])
        q_net = (fields['OCN_SWNET'] + fields['OCN_LWNET'] - fields['OCN_SENS']
                 - HLV * fields['OCN_EVAP'])
        sst = np.where(ofrac > 0., sst + q_net * dt / heat_capacity, sst)
        # SST at date+dt; the LAG in the namcouple delivers it to Isca at date+dt.
        sst_var.put(date, pyoasis.asarray(sst))
        frac_var.put(date, pyoasis.asarray(ofrac))

    print('slab_ocean: finished, global mean SST %.4f K' % np.mean(sst[ofrac > 0]), flush=True)
    del comp


if __name__ == '__main__':
    main()

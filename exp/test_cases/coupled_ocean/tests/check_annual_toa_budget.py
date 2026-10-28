"""Annual energy budget from the top of the atmosphere, for runs without surface radiation diagnostics.

Over a whole year the atmosphere stores almost no energy, so the global net downward flux at
the TOA should equal the change in MOM6's heat content (ocean.stats) plus the change in
SIS2's ice enthalpy (seaice.stats, if present), up to what land and Isca's slab cells store.
Global W/m2 of planet, per year of 12 runs.

Usage: python check_annual_toa_budget.py EXP
"""
import glob
import os
import sys

import numpy as np
import xarray as xr

from check_sis2_budget import ice_enthalpy

D = os.environ['GFDL_DATA']
R = 6376.0e3
At = 4 * np.pi * R**2


def toa_net(a):
    # RRTM only: the grey scheme does not output the reflected SW, which sea ice makes uneven
    # (use check_sis2_budget.py, the surface budget, for grey runs)
    return a.toa_sw - a.olr


def heat(run):
    st = xr.open_dataset(os.path.join(run, 'ocean', 'ocean.stats.nc'), decode_times=False)
    return st.Time.values, st.Heat.values


exp = sys.argv[1]
runs = sorted(glob.glob(os.path.join(D, exp, 'run[0-9]*')))
print('year   TOA    dOHC/dt  dIce/dt  residual  [W/m2 of planet]')
for y in range(len(runs) // 12):
    yr = runs[12 * y:12 * y + 12]
    a = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in yr], decode_times=False,
                          combine='by_coords').mean('time')
    _, w = np.polynomial.legendre.leggauss(a.sizes['lat'])
    toa = float((toa_net(a).mean('lon') * w).sum() / w.sum())
    d0, h0 = heat(yr[0]); d1, h1 = heat(yr[-1])
    dt = (d1[-1] - d0[0]) * 86400.
    docn = (h1[-1] - h0[0]) / dt / At
    dice = 0.
    if os.path.exists(os.path.join(yr[0], 'ocean', 'seaice.stats')):
        i0 = ice_enthalpy(os.path.join(yr[0], 'ocean', 'seaice.stats'))[0]
        i1 = ice_enthalpy(os.path.join(yr[-1], 'ocean', 'seaice.stats'))[-1]
        dice = (i1[1] - i0[1]) / dt / At
    print('%4d %7.2f %8.2f %8.2f %9.2f' % (y + 1, toa, docn, dice, toa - docn - dice))

"""Heat and water sent by Isca's ocean tiles against what MOM6 receives, per month.

With surface_tiles the cpl_* diagnostics are the ocean tile's fluxes per unit area of
ocean, so Isca sends cpl_* x coupled_ocean_frac x cell area; MOM6's net_heat_coupler
and its water fluxes are what arrived. Global totals in TW and 1e6 kg/s.

Usage: python check_tile_conservation.py EXP
"""
import glob
import os
import sys

import numpy as np
import xarray as xr

D = os.environ['GFDL_DATA']
R = 6376.0e3
exp = sys.argv[1]
for run in sorted(glob.glob(os.path.join(D, exp, 'run[0-9]*'))):
    a = xr.open_dataset(os.path.join(run, 'atmos_monthly.nc'), decode_times=False).isel(time=0)
    o = xr.open_dataset(os.path.join(run, 'ocean', 'ocean_month.nc'), decode_times=False).isel(time=0)
    ao = xr.open_dataset(os.path.join(run, 'ocean', 'ocean_static.nc'), decode_times=False).area_t.values
    nlat, nlon = a.sizes['lat'], a.sizes['lon']
    _, w = np.polynomial.legendre.leggauss(nlat)
    area = (2 * np.pi * R**2 * w / nlon)[:, None] * np.ones(nlon)
    f = a.coupled_ocean_frac.values
    heat_sent = np.sum(a.cpl_net_heat.values * f * area) / 1e12
    heat_recv = np.nansum(o.net_heat_coupler.values * ao) / 1e12
    pme_sent = np.sum((a.cpl_precip.values - a.cpl_evap.values) * f * area) / 1e6
    runoff_sent = np.sum(a.cpl_runoff.values * area) / 1e6 if 'cpl_runoff' in a else 0.
    water_recv = np.nansum(o.PRCmE.values * ao) / 1e6 if 'PRCmE' in o else np.nan
    print('%s heat sent %9.3f recv %9.3f TW (%.2e) | P-E sent %8.3f + runoff %7.3f = %8.3f, recv %8.3f'
          % (os.path.basename(run), heat_sent, heat_recv, (heat_recv - heat_sent) / abs(heat_sent),
             pme_sent, runoff_sent, pme_sent + runoff_sent, water_recv))

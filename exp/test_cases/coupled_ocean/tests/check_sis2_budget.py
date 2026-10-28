"""Monthly energy budget of a run with SIS2 sea ice (SIS2=1).

What Isca's surface takes from the atmosphere over the cells the ocean model covers
(net SW + LW down - sensible - latent, the cell means over every tile, so exact only where
cells are wholly ocean-model ocean) against the change in MOM6's heat content
(ocean.stats) plus the change in SIS2's ice enthalpy (seaice.stats, relative to liquid
water at 0 C, so it includes the latent heat of the ice). Global W/m2 of planet.

Usage: python check_sis2_budget.py EXP
"""
import glob
import os
import sys

import numpy as np
import xarray as xr

D = os.environ['GFDL_DATA']
R = 6376.0e3
At = 4 * np.pi * R**2


def ice_enthalpy(path):
    """(day, total enthalpy [J]) of each line of SIS2's seaice.stats."""
    out = []
    for line in open(path):
        if 'Enth' not in line:
            continue
        day = float(line.split(',')[1])
        e = line.split('Enth')[1].split(',')[0].split()
        out.append((day, float(e[0]) + float(e[1])))
    return out


exp = sys.argv[1]
print('run      sfc_net  dOHC/dt  dIce/dt  residual  [W/m2 of planet]')
for run in sorted(glob.glob(os.path.join(D, exp, 'run[0-9]*'))):
    o = os.path.join(run, 'ocean')
    if not os.path.exists(os.path.join(o, 'seaice.stats')):
        continue
    a = xr.open_dataset(os.path.join(run, 'atmos_monthly.nc'), decode_times=False).isel(time=0)
    nlat, nlon = a.sizes['lat'], a.sizes['lon']
    _, w = np.polynomial.legendre.leggauss(nlat)
    area = (2 * np.pi * R**2 * w / nlon)[:, None] * np.ones(nlon)
    f = a.coupled_ocean_frac.values
    net = (a.swdn_sfc + a.lwdn_sfc - a.lwup_sfc - a.flux_t - a.flux_lhe).values
    sfc = np.sum(np.where(f > 0, net * f, 0.) * area) / At
    st = xr.open_dataset(os.path.join(o, 'ocean.stats.nc'), decode_times=False)
    days = st.Time.values if 'Time' in st else np.arange(st.sizes[list(st.dims)[0]])
    heat = st.Heat.values
    dt = (days[-1] - days[0]) * 86400.
    docn = (heat[-1] - heat[0]) / dt / At
    ie = ice_enthalpy(os.path.join(o, 'seaice.stats'))
    dice = (ie[-1][1] - ie[0][1]) / ((ie[-1][0] - ie[0][0]) * 86400.) / At if len(ie) > 1 else np.nan
    print('%s %8.2f %8.2f %8.2f %9.2f' % (os.path.basename(run), sfc, docn, dice, sfc - docn - dice))

"""Energy budget of a coupled run with thermodynamic sea ice (THERMO_ICE=1), per year.

Global means [W/m2 of planet]:
  Isca:  TOA net down; surface net down (what the atmosphere gives the surface)
  MOM6:  heat it received from Isca (net_heat_coupler); its total surface heat flux
         (hfds), whose excess over the received flux is the frazil heat plus the heat
         carried by the water fluxes; the frazil heat Isca received (cpl_frazil), which
         should now match MOM6's frazil and be turned into ice
  ice:   latent heat change of the ice, -rho_i L dh/dt, from the change in annual-mean
         thickness (approximate); slab caps: heat capacity * d(t_ml)/dt (approximate)
Surface net down should equal ocean storage + ice + slab-cap storage.

Usage: python check_ice_budget.py EXP
"""
import glob
import os
import sys

import numpy as np
import xarray as xr

D = os.environ['GFDL_DATA']
exp = sys.argv[1]
R = 6376.0e3
L_ICE = 917. * 3.34e5           # J/m3 (DENS_ICE * HLF_WATER)
SLAB_CAP = 20. * 1035. * 3989.24495292815   # J/m2/K (20 m caps; see frierson_mom6_aquaplanet.py)
YEAR_S = 360 * 86400.

runs = sorted(glob.glob(os.path.join(D, exp, 'run[0-9]*')))
ny = len(runs) // 12
At = 4 * np.pi * R**2
prev = None
print('year   TOA   sfc_net | MOM6 recv  hfds  hfds-recv  frazil(Isca) | dOHC/dt  ice  caps | residual sfc-(storage)')
for y in range(ny):
    rr = runs[12 * y:12 * y + 12]
    a = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in rr], decode_times=False,
                          combine='by_coords', data_vars='minimal', coords='minimal',
                          compat='override').mean('time').load()
    o = xr.open_mfdataset([os.path.join(r, 'ocean', 'ocean_month.nc') for r in rr], decode_times=False,
                          combine='by_coords', data_vars='minimal', coords='minimal',
                          compat='override')[['hfds', 'net_heat_coupler']].mean('time').load()
    ao = xr.open_dataset(os.path.join(rr[0], 'ocean', 'ocean_static.nc'), decode_times=False).area_t.values
    nlat, nlon = a.sizes['lat'], a.sizes['lon']
    _, w = np.polynomial.legendre.leggauss(nlat)
    area = (2 * np.pi * R**2 * w / nlon)[:, None] * np.ones(nlon)
    g = lambda x: float(np.nansum(np.asarray(x) * area) / At)
    alb = 0.31
    toa = g(a.swdn_toa - a.swdn_sfc / (1 - alb) * alb - a.olr)
    sfc = g(a.swdn_sfc + a.lwdn_sfc - a.lwup_sfc - a.flux_t - a.flux_lhe)
    recv = float(np.nansum(o.net_heat_coupler.values * ao) / At)
    hfds = float(np.nansum(o.hfds.values * ao) / At)
    frz = g(a.cpl_frazil) if 'cpl_frazil' in a else np.nan   # W/m2 of Isca cell
    # ocean heat content from ocean.stats at the end of the year and of the previous one
    st = xr.open_dataset(os.path.join(rr[-1], 'ocean', 'ocean.stats.nc'), decode_times=False)
    heat = float(st.Heat[-1])
    h_mean = (a.h_thermo_ice * area).values
    caps = ((a.t_ml * (a.coupled_ocean_frac < 0.5)) * area).values
    if prev is not None:
        docn = (heat - prev[0]) / YEAR_S / At
        dice = -L_ICE * (np.nansum(h_mean) - prev[1]) / YEAR_S / At
        dcap = SLAB_CAP * (np.nansum(caps) - prev[2]) / YEAR_S / At
        print('%4d %6.2f %7.2f | %8.2f %6.2f %8.2f %10.2f | %7.2f %5.2f %5.2f | %8.2f'
              % (y + 1, toa, sfc, recv, hfds, hfds - recv, frz, docn, dice, dcap, sfc - (docn + dice + dcap)))
    else:
        print('%4d %6.2f %7.2f | %8.2f %6.2f %8.2f %10.2f | (storage terms need two years)'
              % (y + 1, toa, sfc, recv, hfds, hfds - recv, frz))
    prev = (heat, np.nansum(h_mean), np.nansum(caps))

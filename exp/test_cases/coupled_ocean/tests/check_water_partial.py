"""Water/heat budget on a non-matched grid: Isca sent (weighted by ocean fraction, plus
routed runoff) vs MOM6 received. Usage: python check_water_partial.py EXP MONTH"""
import os, sys
import numpy as np, xarray as xr
D = os.environ['GFDL_DATA']; r = '%s/%s/run%04d' % (D, sys.argv[1], int(sys.argv[2]))
a = xr.open_dataset(r + '/atmos_monthly.nc', decode_times=False).isel(time=0)
o = xr.open_dataset(r + '/ocean/ocean_month.nc', decode_times=False).isel(time=0)
ao = xr.open_dataset(r + '/ocean/ocean_static.nc', decode_times=False).area_t.values
nlat, nlon = a.sizes['lat'], a.sizes['lon']; _, w = np.polynomial.legendre.leggauss(nlat)
area = (2 * np.pi * 6376e3**2 * w / nlon)[:, None] * np.ones(nlon); frac = a.coupled_ocean_frac.values
pme = (a.precipitation - a.flux_lhe / 2.5e6).values
print('water [1e6 kg/s]: Isca global P-E %.3f | sent (P-E)*frac + runoff %.3f | MOM6 PRCmE %.3f' % (
    np.nansum(pme * area) / 1e6,
    np.nansum(((a.cpl_precip - a.cpl_evap).values * frac + a.cpl_runoff.values) * area) / 1e6,
    np.nansum(o.PRCmE.values * ao) / 1e6))
print('heat [TW]: sent*frac %.3f | MOM6 net_heat_coupler %.3f' % (
    np.nansum(a.cpl_net_heat.values * area * frac) / 1e12, np.nansum(o.net_heat_coupler.values * ao) / 1e12))
print('runoff [1e6 kg/s]: Isca cpl_runoff %.3f | MOM6 friver %.3f' % (
    np.nansum(a.cpl_runoff.values * area) / 1e6, np.nansum(o.friver.values * ao) / 1e6))

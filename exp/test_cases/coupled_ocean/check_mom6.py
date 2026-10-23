"""Sanity checks of an Isca-MOM6 coupled aquaplanet run.

1. The SST Isca uses equals the MOM6 SST (the regridding is the identity).
2. The coupled-ocean fraction is 1 in the MOM6 domain and 0 in the polar caps.
3. Wind stress has the right sign: MOM6 taux follows Isca's surface westerlies.
4. Energy: heat flux MOM6 receives vs the net surface flux Isca diagnoses,
   and MOM6 heat content change vs its time-integrated surface heat flux.

Usage: python check_mom6.py EXP_NAME first_month last_month [outdir]
"""
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = os.environ['GFDL_DATA']
exp, m0, m1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
outdir = sys.argv[4] if len(sys.argv) > 4 else '.'
runs = [os.path.join(DATA, exp, 'run%04d' % i) for i in range(m0, m1 + 1)]

atm = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in runs],
                        decode_times=False, combine='by_coords')
ocn_files = [os.path.join(r, 'ocean', 'ocean_month.nc') for r in runs]
ocn = xr.open_mfdataset(ocn_files, decode_times=False, combine='by_coords') if all(
    os.path.exists(f) for f in ocn_files) else None

a = atm.mean('time').load()
za = a.mean('lon')
print('coupled_ocean_frac zonal mean: min %.3f max %.3f' % (float(za.coupled_ocean_frac.min()),
                                                           float(za.coupled_ocean_frac.max())))
print('  rows with frac>0.5: lat %.1f to %.1f' % (float(za.lat.where(za.coupled_ocean_frac > 0.5).min()),
                                                  float(za.lat.where(za.coupled_ocean_frac > 0.5).max())))
isca_net = a.swdn_sfc + a.lwdn_sfc - a.lwup_sfc - a.flux_t - a.flux_lhe

fig, axs = plt.subplots(2, 2, figsize=(12, 8))
if ocn is not None:
    o = ocn.mean('time').load()
    lat_name = [d for d in o.tos.dims if d.startswith('y')][0]
    lon_name = [d for d in o.tos.dims if d.startswith('x')][0]
    zo = o.mean(lon_name)
    olat = o[lat_name].values
    # 1. SST
    t_isca = za.t_surf.interp(lat=olat)
    print('SST: max |Isca t_surf - (MOM6 tos+273.15)| (zonal-mean, ocean rows) = %.4f K'
          % float(abs(t_isca.values - (zo.tos.values + 273.15)).max()))
    axs[0, 0].plot(za.lat, za.t_surf, label='Isca t_surf')
    axs[0, 0].plot(olat, zo.tos + 273.15, '--', label='MOM6 tos')
    axs[0, 0].set_title('SST [K]'); axs[0, 0].legend()
    # 3. stress sign
    usurf = za.ucomp.isel(pfull=-1).interp(lat=olat)
    # taux is on u-points (xq, yh): same latitudes as the tracer points
    taux = o.taux.mean([d for d in o.taux.dims if d.startswith('x')]).values
    corr = np.corrcoef(usurf.values, taux)[0, 1]
    print('corr(Isca lowest-level u, MOM6 taux) = %.3f  (should be strongly positive)' % corr)
    ax = axs[0, 1]
    ax.plot(olat, usurf / abs(usurf).max(), label='Isca u (lowest level), scaled')
    ax.plot(olat, taux / abs(taux).max(), "--", label="MOM6 taux, scaled")
    ax.set_title('wind vs stress'); ax.legend()
    # 4. energy
    hf_isca = isca_net.mean('lon').interp(lat=olat)
    w = np.cos(np.radians(olat))
    print('ocean-domain mean net surface heat flux: Isca diag %.3f W/m2, MOM6 hfds %.3f W/m2'
          % (float((hf_isca * w).sum() / w.sum()), float((zo.hfds * w).sum() / w.sum())))
    if 'frazil_heat_tend_2d' in zo:
        print('  MOM6 frazil heating (ocean mean) %.3f W/m2'
              % float((zo.frazil_heat_tend_2d * w).sum() / w.sum()))
    axs[1, 0].plot(olat, hf_isca, label='Isca net sfc flux')
    axs[1, 0].plot(olat, zo.hfds, '--', label='MOM6 hfds')
    axs[1, 0].set_title('net surface heat flux [W/m2]'); axs[1, 0].legend()
    if 'T_ady_2d' in o:
        # northward ocean heat transport [PW] from advective + diffusive vertically integrated fluxes
        yq = [d for d in o.T_ady_2d.dims if d.startswith('y')][0]
        oht = (o.T_ady_2d + o.get('T_diffy_2d', 0.)).sum([d for d in o.T_ady_2d.dims if d.startswith('x')]) / 1e15
        axs[1, 1].plot(o[yq], oht)
        axs[1, 1].set_title('ocean heat transport [PW]')
        print('peak ocean heat transport: NH %.3f PW, SH %.3f PW' % (float(oht.max()), float(oht.min())))
fig.tight_layout()
fig.savefig(os.path.join(outdir, '%s_checks.png' % exp), dpi=110)
print('figure written to', os.path.join(outdir, '%s_checks.png' % exp))

# ---- water budget: everything Isca rains/evaporates must reach the ocean -------
HLV = 2.5e6
if ocn is not None and 'PRCmE' in o:
    static = xr.open_dataset(os.path.join(runs[0], 'ocean', 'ocean_static.nc'), decode_times=False)
    area_o = static.area_t.values
    R = 6376.0e3
    latb = np.degrees(np.arcsin(np.concatenate([[-1.], -1. + np.cumsum(np.polynomial.legendre.leggauss(a.sizes['lat'])[1])])))
    area_a = (R**2 * np.radians(360. / a.sizes['lon']) * np.diff(np.sin(np.radians(latb))))[:, None] \
             * np.ones(a.sizes['lon'])[None, :]
    pme_isca = (a.precipitation - a.flux_lhe / HLV).values             # kg/m2/s, all Isca cells
    frac = a.coupled_ocean_frac.values
    tot_isca = np.nansum(pme_isca * area_a)
    caps_isca = np.nansum(pme_isca * area_a * (frac < 0.5))
    tot_mom6 = np.nansum(o.PRCmE.values * area_o)
    riv_mom6 = np.nansum(o.friver.values * area_o) if 'friver' in o else np.nan
    print('freshwater [1e6 kg/s]: Isca global P-E %.4f | MOM6 total in (PRCmE) %.4f'
          % (tot_isca / 1e6, tot_mom6 / 1e6))
    print('                       Isca P-E over uncoupled cells %.4f | MOM6 river input %.4f'
          % (caps_isca / 1e6, riv_mom6 / 1e6))
    print('  (for scale: global precipitation %.2f)' % (np.nansum(a.precipitation.values * area_a) / 1e6))

# ---- exact budget: what Isca sent vs what MOM6 received ------------------------
if ocn is not None and 'cpl_net_heat' in a:
    cpl = frac > 0.5
    heat_sent = np.nansum(a.cpl_net_heat.values * area_a * cpl)
    print('heat [TW]: Isca sent (coupled cells) %.4f | MOM6 hfds %.4f | MOM6 net_heat_surface %.4f'
          % (heat_sent / 1e12, np.nansum(o.hfds.values * area_o) / 1e12,
             np.nansum(o.net_heat_surface.values * area_o) / 1e12 if 'net_heat_surface' in o else np.nan))
    for name in ('net_heat_coupler', 'heat_content_massout', 'heat_content_massin', 'frazil_heat_tend_2d'):
        if name in o:
            print('          MOM6 %-22s %.4f TW' % (name, np.nansum(o[name].values * area_o) / 1e12))
    water_sent = np.nansum(((a.cpl_precip - a.cpl_evap).values * cpl + a.cpl_runoff.values) * area_a)
    print('water [1e6 kg/s]: Isca sent P-E+runoff %.4f | MOM6 PRCmE %.4f'
          % (water_sent / 1e6, np.nansum(o.PRCmE.values * area_o) / 1e6))

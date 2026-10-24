"""T5: sensitivity of the Isca-MOM6 aquaplanet to the coupling period.

Compares year-2 annual, zonal means of runs that differ only in coupling period
with the 1-hour reference (S1, same configuration and initial state). The S1
year-2 minus year-3 difference gives a scale for drift plus internal variability.

Usage: python compare_coupling_period.py REF EXP [EXP ...] [--out fig.png]
"""
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = os.environ['GFDL_DATA']
args = [a for a in sys.argv[1:] if not a.startswith('--')]
out = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--out=')), 'T5_coupling_period.png')
ref, exps = args[0], args[1:]
R = 6376.0e3


def year(exp, y):
    rr = [os.path.join(DATA, exp, 'run%04d' % m) for m in range(12 * (y - 1) + 1, 12 * y + 1)]
    a = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in rr], decode_times=False,
                          combine='by_coords').mean('time').load()
    o = xr.open_mfdataset([os.path.join(r, 'ocean', 'ocean_month.nc') for r in rr], decode_times=False,
                          combine='by_coords').mean('time').load()
    return a, o


def summary(a, o):
    _, w = np.polynomial.legendre.leggauss(a.sizes['lat'])
    za = a.mean('lon')
    gm = lambda f: float((f.mean('lon') * w).sum() / w.sum())
    # ocean heat transport from the zonally integrated MOM6 surface heat flux (no polar caps)
    hf = o.hfds.fillna(0.)
    yname = [d for d in hf.dims if d.startswith('y')][0]
    lat = o[yname].values
    dlat = np.deg2rad(np.gradient(lat))
    zint = (hf.mean([d for d in hf.dims if d.startswith('x')]) * 2 * np.pi * R**2
            * np.cos(np.deg2rad(lat)) * dlat).values
    zint -= zint.sum() * np.cos(np.deg2rad(lat)) * dlat / (np.cos(np.deg2rad(lat)) * dlat).sum()
    oht = np.cumsum(zint) / 1e15
    return dict(lat=za.lat.values, sst=za.t_surf.values, pr=za.precipitation.values * 86400,
                u=za.ucomp.isel(pfull=-1).values, olat=lat, oht=oht,
                gm_sst=gm(a.t_surf), gm_pr=gm(a.precipitation) * 86400,
                max_oht=float(np.abs(oht).max()), itcz=float(za.precipitation.max() * 86400))


S = {ref + ' y2': summary(*year(ref, 2)), ref + ' y3': summary(*year(ref, 3))}
for e in exps:
    S[e + ' y2'] = summary(*year(e, 2))
r2 = S[ref + ' y2']
print('%-22s %8s %8s %8s %9s %7s %7s' % ('run', 'SST[K]', 'P[mm/d]', 'ITCZ', 'maxOHT', 'dSSTrms', 'dPrms'))
for k, s in S.items():
    print('%-22s %8.3f %8.3f %8.2f %9.3f %7.3f %7.3f' % (
        k, s['gm_sst'], s['gm_pr'], s['itcz'], s['max_oht'],
        np.sqrt(np.mean((s['sst'] - r2['sst'])**2)), np.sqrt(np.mean((s['pr'] - r2['pr'])**2))))

fig, axs = plt.subplots(2, 2, figsize=(12, 8))
for k, s in S.items():
    ls = '--' if k.endswith('y3') else '-'
    axs[0, 0].plot(s['lat'], s['sst'] - r2['sst'], ls, label=k)
    axs[0, 1].plot(s['lat'], s['pr'], ls, label=k)
    axs[1, 0].plot(s['lat'], s['u'], ls, label=k)
    axs[1, 1].plot(s['olat'], s['oht'], ls, label=k)
for ax, t in zip(axs.flat, ['SST minus %s y2 [K]' % ref, 'precipitation [mm/day]',
                            'lowest-level u [m/s]', 'ocean heat transport [PW]']):
    ax.set_title(t); ax.set_xlabel('latitude'); ax.grid(alpha=.3)
axs[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(out, dpi=110)
print('wrote', out)

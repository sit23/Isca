"""Compare the Phase-1 runs: Isca's built-in slab vs the same slab coupled through OASIS.

Usage: python check_phase1.py [first_month] [last_month] [outdir]
"""
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = os.environ['GFDL_DATA']
m0, m1 = int(sys.argv[1]), int(sys.argv[2])
outdir = sys.argv[3] if len(sys.argv) > 3 else '.'


def load(name):
    files = [os.path.join(DATA, name, 'run%04d' % i, 'atmos_monthly.nc') for i in range(m0, m1 + 1)]
    ds = xr.open_mfdataset(files, decode_times=False, combine='by_coords')
    return ds.mean('time').mean('lon').load(), ds.mean('lon').load()


ref, ref_t = load('coupled_ocean_phase1_reference')
cpl, cpl_t = load('coupled_ocean_phase1_slab')
w = np.cos(np.radians(ref.lat))

fields = [('t_surf', 'K'), ('precipitation', 'mm/day'), ('flux_lhe', 'W/m2'), ('flux_t', 'W/m2')]
scale = {'precipitation': 86400.}
fig, axs = plt.subplots(2, 3, figsize=(14, 8))
for ax, (f, units) in zip(axs.flat, fields):
    s = scale.get(f, 1.)
    ax.plot(ref.lat, ref[f] * s, label='built-in slab')
    ax.plot(cpl.lat, cpl[f] * s, '--', label='OASIS-coupled slab')
    d = (cpl[f] - ref[f]) * s
    # month-to-month spread of the reference as a yardstick for internal variability
    spread = (ref_t[f] * s).std('time') / np.sqrt(ref_t.sizes['time'])
    print('%-14s global mean ref %.4f cpl %.4f | max |zonal-mean diff| %.4f, typical sampling error %.4f %s'
          % (f, float((ref[f] * s * w).sum() / w.sum()), float((cpl[f] * s * w).sum() / w.sum()),
             float(abs(d).max()), float(spread.mean()), units))
    ax.set_title(f); ax.set_ylabel(units); ax.legend(fontsize=8)
for ax, f in zip(axs.flat[4:], ['ucomp', 'temp']):
    d = (cpl[f] - ref[f])
    c = ax.contourf(ref.lat, ref.pfull, d, 21, cmap='RdBu_r')
    ax.contour(ref.lat, ref.pfull, ref[f], 10, colors='k', linewidths=0.5)
    ax.invert_yaxis(); plt.colorbar(c, ax=ax); ax.set_title('%s: coupled - built-in' % f)
    print('%-14s max |zonal-mean diff| %.4f' % (f, float(abs(d).max())))
fig.suptitle('Phase 1: months %d-%d' % (m0, m1))
fig.tight_layout()
fig.savefig(os.path.join(outdir, 'phase1_slab_comparison.png'), dpi=110)
print('figure written to', os.path.join(outdir, 'phase1_slab_comparison.png'))

"""Science sanity checks of the long coupled runs (S1 aquaplanet, S3 ridge world).

* spin-up: global TOA imbalance, mean SST and volume-mean ocean temperature per year
* annual-mean zonal-mean SST, ocean heat transport (OHT) and atmospheric heat
  transport (AHT), compared with the slab-ocean reference (R1, no OHT):
  Bjerknes compensation means AHT should drop where OHT is large while the
  total stays similar (Czaja & Marshall 2006; Enderton & Marshall 2009)
* barotropic streamfunction (ridge world: gyres, western boundary currents)

Usage: python analyse_long_runs.py EXP [REF_EXP] [last_n_years] [outdir]
"""
import glob
import os
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = os.environ['GFDL_DATA']
exp = sys.argv[1]
ref = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != '-' else None
nyears = int(sys.argv[3]) if len(sys.argv) > 3 else 5
outdir = sys.argv[4] if len(sys.argv) > 4 else '.'
R = 6376.0e3


def runs(name):
    return sorted(glob.glob(os.path.join(DATA, name, 'run[0-9]*')))


def gaussian_area(nlat, nlon):
    _, w = np.polynomial.legendre.leggauss(nlat)
    return (2 * np.pi * R**2 * w / nlon)[:, None] * np.ones(nlon)[None, :]


def isca_net_sfc_down(a):
    """Net surface energy flux into the surface [W/m2] from grey-radiation diagnostics."""
    return a.swdn_sfc + a.lwdn_sfc - a.lwup_sfc - a.flux_t - a.flux_lhe


def toa_net_down(a, albedo=None):
    """Net downward TOA flux [W/m2] for the Frierson grey scheme: SW reflected at the
    surface escapes without absorption, so absorbed SW = swdn_toa - swdn_sfc_down*albedo."""
    if 'toa_sw' in a:
        return a.toa_sw - a.olr
    alb = 0.31 if albedo is None else albedo
    sw_sfc_down = a.swdn_sfc / (1. - alb)
    return a.swdn_toa - sw_sfc_down * alb - a.olr


def northward_transport(heating, area):
    """Northward transport [PW] implied by column heating (lat, lon) [W/m2], with the
    global-mean imbalance removed so the transport vanishes at both poles."""
    hz = (heating * area).sum(axis=1)
    hz = hz - hz.sum() * area.sum(axis=1) / area.sum()
    return np.concatenate([[0.], np.cumsum(hz)]) / 1e15


rr = runs(exp)
nm = len(rr)
print('%s: %d months' % (exp, nm))
atm_all = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in rr], decode_times=False,
                            combine='by_coords', data_vars='minimal', coords='minimal', compat='override')
nlat, nlon = atm_all.sizes['lat'], atm_all.sizes['lon']
area = gaussian_area(nlat, nlon)
_, wts = np.polynomial.legendre.leggauss(nlat)
latb = np.degrees(np.arcsin(np.concatenate([[-1.], -1. + np.cumsum(wts)]).clip(-1, 1)))

# ---- spin-up time series (annual means)
toa = (toa_net_down(atm_all) * area).sum(['lat', 'lon']) / area.sum()
sst = (atm_all.t_surf * area).sum(['lat', 'lon']) / area.sum()
ny = nm // 12
toa_y = toa.values[:ny * 12].reshape(ny, 12).mean(1)
sst_y = sst.values[:ny * 12].reshape(ny, 12).mean(1)
ocn_T = []
for r in rr[11::12]:
    s = os.path.join(r, 'ocean', 'ocean.stats.nc')
    if os.path.exists(s):
        st = xr.open_dataset(s, decode_times=False)
        if 'Heat' in st and 'Mass' in st:
            ocn_T.append(float(st.Heat[-1] / st.Mass[-1] / 3989.24495292815))
print('year  TOA[W/m2]  SST[K]   ocean mean T [C]')
for y in range(ny):
    print('%4d  %8.3f  %7.3f  %s' % (y + 1, toa_y[y], sst_y[y], '%.4f' % ocn_T[y] if y < len(ocn_T) else ''))

# ---- last n years: climate and transports
last = rr[-12 * min(nyears, max(ny, 1)):]
a = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in last], decode_times=False,
                      combine='by_coords', data_vars='minimal', coords='minimal', compat='override').mean('time').load()
toa_net = toa_net_down(a).values
sfc_net = isca_net_sfc_down(a).values
# total transport from TOA; atmospheric from TOA minus surface; ocean = what the surface takes up
total = northward_transport(toa_net, area)
aht = northward_transport(toa_net - sfc_net, area)
oht_implied = northward_transport(sfc_net, area)

fig, axs = plt.subplots(2, 2, figsize=(13, 9))
ax = axs[0, 0]
ax.plot(latb, total, 'k', label='total (TOA)')
ax.plot(latb, aht, 'C0', label='atmosphere')
ax.plot(latb, oht_implied, 'C1', label='ocean (implied by surface flux)')

o = None
ocn_files = [os.path.join(r, 'ocean', 'ocean_month.nc') for r in last]
if all(os.path.exists(f) for f in ocn_files):
    o = xr.open_mfdataset(ocn_files, decode_times=False, combine='by_coords').mean('time').load()
    if 'T_ady_2d' in o:
        yq = [d for d in o.T_ady_2d.dims if d.startswith('y')][0]
        xd = [d for d in o.T_ady_2d.dims if d.startswith('x')][0]
        oht = (o.T_ady_2d + (o.T_diffy_2d if 'T_diffy_2d' in o else 0.)).sum(xd) / 1e15
        ax.plot(o[yq], oht, 'C1--', label='ocean (MOM6 advective+diffusive)')
        print('OHT (MOM6) peak: NH %.2f PW at %.0f, SH %.2f PW at %.0f'
              % (float(oht.max()), float(o[yq][oht.argmax()]), float(oht.min()), float(o[yq][oht.argmin()])))
    if 'thetao' in o and 'hfds' in o:
        # Out of equilibrium the ocean stores heat unevenly, so the transport implied by the
        # surface flux must subtract the actual storage (heat content change over the period;
        # z* layers taken at their rest thickness).
        th = xr.open_mfdataset([ocn_files[0], ocn_files[-1]], decode_times=False,
                               combine='by_coords').thetao.load()
        zi = np.concatenate([[0.], 0.5 * (o.zl.values[1:] + o.zl.values[:-1]), [2 * o.zl.values[-1] - o.zl.values[-2]]])
        rho_cp = 1035. * 3989.24495292815
        hc = (th.fillna(0.) * xr.DataArray(np.diff(zi), dims='zl')).sum('zl') * rho_cp
        storage = (hc.isel(time=-1) - hc.isel(time=0)).values / ((len(ocn_files) - 1) * 30 * 86400.)
        yh, xh = o.yh.values, o.xh.values
        ocn_area = (R**2 * np.cos(np.deg2rad(yh)) * np.deg2rad(np.gradient(yh)))[:, None] * \
            np.deg2rad(np.gradient(xh))[None, :]
        oht_store = np.cumsum(((o.hfds.fillna(0.).values - storage) * ocn_area).sum(1)) / 1e15
        ax.plot(yh, oht_store, 'C3', lw=1, label='ocean (surface flux minus storage)')
        print('ocean heat uptake %.2f W/m2 of ocean, storage %.2f W/m2'
              % ((o.hfds.fillna(0.).values * ocn_area).sum() / ocn_area.sum(),
                 (storage * ocn_area).sum() / ocn_area.sum()))
print('AHT peak: NH %.2f PW, SH %.2f PW; total peak NH %.2f PW' % (aht.max(), aht.min(), total.max()))

if ref:
    rref = runs(ref)
    ar = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in rref[-12 * nyears:]],
                           decode_times=False, combine='by_coords', data_vars='minimal', coords='minimal',
                           compat='override').mean('time').load()
    aht_ref = northward_transport(toa_net_down(ar).values - isca_net_sfc_down(ar).values, area)
    ax.plot(latb, aht_ref, 'C0:', label='atmosphere, slab-ocean reference')
    print('slab reference AHT peak: NH %.2f PW, SH %.2f PW' % (aht_ref.max(), aht_ref.min()))
    j = np.argmin(abs(latb - 35))
    print('at 35N: coupled AHT %.2f + OHT %.2f = %.2f PW vs slab AHT %.2f PW'
          % (aht[j], oht_implied[j], total[j], aht_ref[j]))
    axs[0, 1].plot(ar.lat, ar.t_surf.mean('lon'), 'k:', label='slab reference')
ax.axhline(0, color='0.7', lw=0.5); ax.set_xlabel('lat'); ax.set_ylabel('PW'); ax.legend(fontsize=8)
ax.set_title('northward heat transport, last %d yr' % nyears)

ax = axs[0, 1]
ax.plot(a.lat, a.t_surf.mean('lon'), 'k', label='coupled')
ax.set_title('zonal-mean surface temperature [K]'); ax.legend(fontsize=8)

ax = axs[1, 0]
ax.plot(np.arange(1, ny + 1), toa_y, 'k')
ax.set_xlabel('year'); ax.set_ylabel('W/m2'); ax.set_title('global TOA net downward flux')
if ocn_T:
    ax2 = ax.twinx(); ax2.plot(np.arange(1, len(ocn_T) + 1), ocn_T, 'C3'); ax2.set_ylabel('ocean mean T [C]', color='C3')

ax = axs[1, 1]
if o is not None and 'vmo' in o:
    # barotropic streamfunction: zonal integral (from 0E, the ridge if there is one) of the
    # depth-summed meridional mass transport vmo [kg/s] -> Sv
    yq = [d for d in o.vmo.dims if d.startswith('y')][0]
    xh = [d for d in o.vmo.dims if d.startswith('x')][0]
    zl = [d for d in o.vmo.dims if d.startswith('z')][0]
    psi = o.vmo.sum(zl).cumsum(xh) / 1035. / 1e6
    lim = float(abs(psi).max())
    c = ax.contourf(o[xh], o[yq], psi, np.linspace(-lim, lim, 21), cmap='RdBu_r')
    plt.colorbar(c, ax=ax)
    ax.set_title('barotropic streamfunction [Sv]')
    print('barotropic streamfunction range %.1f to %.1f Sv' % (float(psi.min()), float(psi.max())))
fig.suptitle(exp)
fig.tight_layout()
fn = os.path.join(outdir, '%s_science.png' % exp)
fig.savefig(fn, dpi=110)
print('figure written to', fn)

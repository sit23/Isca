"""Lat-lon maps of the last-year-mean SST of each coupled experiment, with the land of
both models, the routed rivers and the cells where Isca's land/ocean mask and the ocean
model's coverage disagree.

SST is MOM6's tos on its own grid (shown by nearest-neighbour lookup on a 0.5 deg plotting
grid), or for Samudra, Isca's t_surf on the cells Samudra couples (zero coupling lag).

Usage: python plot_sst_maps.py [outdir] [EXP ...]
"""
import os
import sys

import netCDF4
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

from isca import GFDL_BASE
from isca.coupling import gaussian_latitudes, runoff_destinations

DATA = os.environ['GFDL_DATA']
NLAT, NLON = 64, 128
LAT = gaussian_latitudes(NLAT)
LON = np.arange(NLON) * 360. / NLON
_, W = np.polynomial.legendre.leggauss(NLAT)
EDGES = np.degrees(np.arcsin(np.concatenate([[-1.], -1. + np.cumsum(W)]).clip(-1, 1)))
AREA = (2 * np.pi * 6376e3**2 * W / NLON)[:, None] * np.ones(NLON)

with netCDF4.Dataset(os.path.join(GFDL_BASE, 'input/land_masks/era_land_t42.nc')) as nc:
    EARTH_LAND = np.asarray(nc['land_mask'][:]) > 0.
    EARTH_ZSURF = np.asarray(nc['zsurf'][:])
RIDGE = np.zeros((NLAT, NLON), dtype=bool)
RIDGE[:, [0, NLON - 1]] = True

# name: title, Isca land, Isca topography for routing, ocean domain max_lat, river mode
# [, longitude to centre the map on]
EXPS = {
    'S1_aquaplanet_mom6_50yr': ('S1 aquaplanet, MOM6 to $\\pm$76.7$^\\circ$, year 50', None, None, 76.7, None),
    'T7_maxlat70': ('T7 aquaplanet, MOM6 to $\\pm$70$^\\circ$, year 2', None, None, 70., None),
    'T6_regular2deg': ('T6 aquaplanet, MOM6 on a regular 2$^\\circ$ grid to $\\pm$80$^\\circ$, year 2', None, None, 80., None),
    'S3_ridge_mom6_50yr': ('S3 ridge world, year 50', RIDGE, None, 76.7, 'all land P$-$E routed (old code)', 0.),
    'E1b_realistic_mom6': ('E1b Earth-like MOM6, no land bucket, year 20', EARTH_LAND, EARTH_ZSURF, 78., 'none: no bucket, so nothing routed'),
    'E1e_realistic_mom6_bucket': ('E1e Earth-like MOM6 with rivers, year 20', EARTH_LAND, EARTH_ZSURF, 78., 'bucket overflow'),
    'M3_samudra_ra025': ('M3 Earth-like Samudra, year 3', EARTH_LAND, None, None, 'none: Samudra takes no freshwater'),
}


def last_year(exp):
    runs = sorted(d for d in os.listdir(os.path.join(DATA, exp)) if d.startswith('run'))
    return [os.path.join(DATA, exp, r) for r in runs[-12:]]


def domain_mask(max_lat):
    """Isca cells meant to be coupled: more than half inside max_lat (as the experiments set up)."""
    if max_lat is None:
        return np.ones((NLAT, NLON), dtype=bool)
    s, m = np.sin(np.radians(EDGES)), np.sin(np.radians(max_lat))
    covered = (np.minimum(s[1:], m) - np.maximum(s[:-1], -m)).clip(0) / np.diff(s)
    return np.broadcast_to((covered > 0.5)[:, None], (NLAT, NLON)).copy()


def next_step(coupled, zsurf):
    """Downhill neighbour (or nearest coupled cell for pits) of every uncoupled cell,
    following runoff_destinations step by step. Returns {(j, i): ((j2, i2), jumped)}."""
    lat, lon = np.radians(LAT), np.radians(LON)
    LO, LA = np.meshgrid(lon, lat)
    xyz = np.stack([np.cos(LA) * np.cos(LO), np.cos(LA) * np.sin(LO), np.sin(LA)], -1)
    oj, oi = np.nonzero(coupled)
    oxyz = xyz[oj, oi]
    step = {}
    for j in range(NLAT):
        for i in range(NLON):
            if coupled[j, i]:
                continue
            best, nxt = (zsurf[j, i] if zsurf is not None else None), None
            if zsurf is not None:
                for dj in (-1, 0, 1):
                    for di in (-1, 0, 1):
                        nj, ni = j + dj, (i + di) % NLON
                        if (dj or di) and 0 <= nj < NLAT and zsurf[nj, ni] < best:
                            best, nxt = zsurf[nj, ni], (nj, ni)
            if nxt is None:
                k = np.argmax(oxyz @ xyz[j, i])
                step[(j, i)] = ((oj[k], oi[k]), True)
            else:
                step[(j, i)] = (nxt, False)
    return step


def regrid_nearest(lon2d, lat2d, values, lon_out, lat_out):
    """Nearest-neighbour lookup of a curvilinear field onto a regular plotting grid."""
    def xyz(lo, la):
        lo, la = np.radians(lo), np.radians(la)
        return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    tree = cKDTree(xyz(lon2d.ravel(), lat2d.ravel()))
    LO, LA = np.meshgrid(lon_out, lat_out)
    _, k = tree.query(xyz(LO.ravel(), LA.ravel()).reshape(-1, 3))
    return values.ravel()[k].reshape(LO.shape)


def model_wet_fraction(lon2d, lat2d, wet, res=0.125):
    """Fraction of each Isca cell covered by the ocean model's wet cells, from a
    nearest-neighbour sampling of its wet mask on a fine grid. (Isca's own
    coupled_ocean_frac is zero on Isca land by construction, so it can't show this.)"""
    flon = np.arange(res / 2, 360., res)
    flat = np.arange(-90 + res / 2, 90., res)
    fine = regrid_nearest(lon2d, lat2d, wet.astype(float), flon, flat)
    lat_ok = (flat >= np.nanmin(lat2d) - 1.) & (flat <= np.nanmax(lat2d) + 1.)
    fine[~lat_ok, :] = 0.
    jj = np.clip(np.searchsorted(EDGES, flat) - 1, 0, NLAT - 1)
    ii = np.floor((flon + 180. / NLON) / (360. / NLON)).astype(int) % NLON
    wts = np.cos(np.radians(flat))[:, None] * np.ones(len(flon))
    num = np.zeros((NLAT, NLON)); den = np.zeros((NLAT, NLON))
    np.add.at(num, (jj[:, None] * np.ones(len(flon), int), ii[None, :] * np.ones((len(flat), 1), int)), fine * wts)
    np.add.at(den, (jj[:, None] * np.ones(len(flon), int), ii[None, :] * np.ones((len(flat), 1), int)), wts)
    return num / den


def isca_cells(ax, mask, offsets=(0.,), **kw):
    """Draw Isca grid cells where mask is true as rectangles."""
    dlon = 360. / NLON
    for off in offsets:
        for j, i in zip(*np.nonzero(mask)):
            ax.add_patch(plt.Rectangle((LON[i] - dlon / 2 + off, EDGES[j]), dlon,
                                       EDGES[j + 1] - EDGES[j], **kw))


def plot(exp, outdir):
    title, land, zsurf, max_lat, river_mode = EXPS[exp][:5]
    centre = EXPS[exp][5] if len(EXPS[exp]) > 5 else 180.
    # everything is drawn in 0-360 coordinates, repeated 360 degrees west when centred on 0
    offsets = (0., -360.) if centre < 90. else (0.,)
    runs = last_year(exp)
    atm = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in runs], decode_times=False,
                            combine='by_coords', data_vars='minimal', coords='minimal',
                            compat='override').mean('time').load()
    frac = atm.coupled_ocean_frac.values
    isca_land = land if land is not None else np.zeros((NLAT, NLON), dtype=bool)
    domain = domain_mask(max_lat)

    plon = np.arange(0.25, 360., 0.5)
    plat = np.arange(-89.75, 90., 0.5)
    if os.path.exists(os.path.join(runs[0], 'ocean', 'ocean_month.nc')):
        st = xr.open_dataset(os.path.join(runs[0], 'ocean', 'ocean_static.nc'), decode_times=False)
        tos = xr.open_mfdataset([os.path.join(r, 'ocean', 'ocean_month.nc') for r in runs],
                                decode_times=False, combine='by_coords', data_vars='minimal',
                                coords='minimal', compat='override').tos.mean('time').values
        sst = regrid_nearest(st.geolon.values % 360., st.geolat.values, tos, plon, plat)
        lat_top, lat_bot = np.nanmax(st.geolat.values), np.nanmin(st.geolat.values)
        sst[(plat > lat_top + 1.) | (plat < lat_bot - 1.), :] = np.nan
        sst_label = 'MOM6 SST (tos) [$^\\circ$C]'
        wetfrac = model_wet_fraction(st.geolon.values % 360., st.geolat.values, st.wet.values > 0)
    else:
        # Isca's t_surf is Samudra's SST only in months without the ice lid: average those
        mon = xr.open_mfdataset([os.path.join(r, 'atmos_monthly.nc') for r in runs], decode_times=False,
                                combine='by_coords', data_vars='minimal', coords='minimal',
                                compat='override')[['t_surf', 'ice_conc']].load()
        free = mon.ice_conc <= 0.5     # Isca's ice_concentration_threshold: the lid is off
        ts = (mon.t_surf.where(free).sum('time') / free.sum('time').where(free.sum('time') > 0)).values
        ts = np.where((frac > 0.5) & ~isca_land, ts - 273.15, np.nan)
        jj = np.searchsorted(EDGES, plat) - 1
        ii = np.round(plon / (360. / NLON)).astype(int) % NLON
        sst = ts[jj.clip(0, NLAT - 1)][:, ii]
        sst_label = 'Samudra SST (Isca t$_{surf}$ on coupled cells, ice-free months) [$^\\circ$C]'
        c = xr.open_dataset('/disca/share/sit204/samudra/om4_sst_clim_2003_2012.nc')
        wetfrac = model_wet_fraction(c.lon.values % 360., c.lat.values, c.mask.values)

    fig, ax = plt.subplots(figsize=(13, 6.6))
    ax.set_facecolor('#d9d9d9')       # ocean model land / outside the ocean domain
    for off in offsets:
        pm = ax.pcolormesh(plon + off, plat, sst, cmap='RdYlBu_r', vmin=-2, vmax=30, shading='nearest')
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(sst_label)

    # Isca land: outlined in black (ridge or continents)
    if land is not None:
        jj = np.clip(np.searchsorted(EDGES, plat) - 1, 0, NLAT - 1)
        ii = np.floor((plon + 180. / NLON) / (360. / NLON)).astype(int) % NLON
        for off in offsets:
            ax.contour(plon + off, plat, isca_land[jj][:, ii].astype(float), levels=[0.5], colors='k',
                       linewidths=0.9)
    # Isca slab outside the ocean domain (polar caps)
    if max_lat is not None:
        caps = ~domain[:, 0]
        for j in np.nonzero(caps)[0]:
            ax.axhspan(EDGES[j], EDGES[j + 1], fill=False, hatch='///', ec='0.4', lw=0)

    # Isca/ocean-model mask mismatches (inside the intended ocean domain)
    ocean_under_land = isca_land & (wetfrac >= 0.25)
    uncoupled_ocean = ~isca_land & domain & (frac <= 0.5)
    partial = ~isca_land & domain & (frac > 0.5) & (frac < 0.99)
    isca_cells(ax, partial, offsets, fill=False, ec='0.35', lw=0.4, ls=':')
    isca_cells(ax, ocean_under_land, offsets, fill=False, ec='m', lw=1.3)
    isca_cells(ax, uncoupled_ocean, offsets, fill=False, ec='lime', lw=1.3)

    handles = [Line2D([], [], color='k', lw=0.9, label='Isca land')] if land is not None else []
    handles += [plt.Rectangle((0, 0), 1, 1, fc='#d9d9d9', label='land / outside the ocean model')]
    if max_lat is not None:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc='w', hatch='///', ec='0.4',
                                     label='Isca slab (outside ocean domain)'))
    handles += [plt.Rectangle((0, 0), 1, 1, fill=False, ec='m', lw=1.3,
                              label='Isca land, ocean model wet over $\\geq$25%% (%d cells)' % ocean_under_land.sum()),
                plt.Rectangle((0, 0), 1, 1, fill=False, ec='lime', lw=1.3,
                              label='Isca ocean, ocean model covers $\\leq$50%% (%d)' % uncoupled_ocean.sum()),
                plt.Rectangle((0, 0), 1, 1, fill=False, ec='0.35', lw=0.4, ls=':',
                              label='Isca ocean, partly covered 50-99%% (%d)' % partial.sum())]

    # rivers: the routing paths Isca sets up, and the mean discharge arriving at each mouth
    active = river_mode is not None and river_mode.startswith(('bucket', 'all'))
    if land is not None and active:
        coupled = domain & ~isca_land
        step = next_step(coupled, zsurf)
        dest = runoff_destinations(coupled, zsurf)
        upstream = np.zeros((NLAT, NLON))
        for (j, i) in step:
            c = (j, i)
            for _ in range(NLAT * NLON):
                upstream[c] += AREA[j, i]
                if coupled[c]:
                    break
                c = step[c][0]
        for (j, i), ((j2, i2), jumped) in step.items():
            x0, x1 = LON[i], LON[i2]
            if abs(x1 - x0) > 180. or upstream[j, i] < 4 * AREA[j, i]:
                continue
            lw = 0.6 * np.sqrt(upstream[j, i] / AREA[j, i])
            for off in offsets:
                ax.plot([x0 + off, x1 + off], [LAT[j], LAT[j2]], color='navy' if active else '0.35',
                        lw=min(lw, 3.), ls=':' if jumped else '-', alpha=0.8 if active else 0.5,
                        solid_capstyle='round')
        handles.append(Line2D([], [], color='navy' if active else '0.35', lw=1.5,
                              label='river channel, $\\geq$4 cells upstream (width ~ $\\sqrt{\\mathrm{catchment}}$; dotted: pit, jump to nearest ocean)'))
        if active and 'cpl_runoff' in atm:
            run = atm.cpl_runoff.values
            spread = np.median(run[(frac > 0.99) & ~isca_land & (dest == 0)]) if 'bucket' in river_mode else 0.
            mouth = (run - spread * np.clip(frac, 0, 1)) * AREA / 1000.     # [m3/s]
            mouths = np.unique(dest[dest > 0]) - 1
            mj, mi = mouths // NLON, mouths % NLON
            q = mouth[mj, mi]
            ok = np.abs(q) > 2e2
            for off in offsets:
                ax.scatter(LON[mi[ok]] + off, LAT[mj[ok]], s=np.sqrt(np.abs(q[ok])) / 6.,
                           c=np.where(q[ok] > 0, 'gold', 'crimson'), edgecolors='k',
                           linewidths=0.5, zorder=5)
            handles.append(Line2D([], [], marker='o', color='w', mfc='gold', mec='k', ms=8,
                                  label='river mouth, area ~ discharge (total %.2f$\\times10^6$ m$^3$/s)'
                                  % (q[q > 0].sum() / 1e6)))
            if (q < -2e2).any():
                handles.append(Line2D([], [], marker='o', color='w', mfc='crimson', mec='k', ms=8,
                                      label='negative "runoff" (E > P routed)'))
    ax.set_title('%s\nrivers: %s' % (title, river_mode) if river_mode else title, fontsize=11)

    ax.set_xlim(centre - 180., centre + 180.)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(centre - 180., centre + 181., 60.))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xlabel('longitude [$^\\circ$E]')
    ax.set_ylabel('latitude [$^\\circ$N]')
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8,
              frameon=False)
    fig.tight_layout()
    fn = os.path.join(outdir, '%s_sst_map.png' % exp)
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('%s: %d Isca-land cells with ocean-model water, %d Isca ocean cells uncoupled, %d partial -> %s'
          % (exp, ocean_under_land.sum(), uncoupled_ocean.sum(), partial.sum(), fn))


if __name__ == '__main__':
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    for exp in (sys.argv[2:] or EXPS):
        plot(exp, outdir)

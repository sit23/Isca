"""River routing on a prescribed real-world drainage network, for Earth-like coupled runs.

Isca's default routing (isca.coupling.runoff_destinations) sends water downhill on the
model's own T42 topography. That suits idealised land, but on Earth it gives the wrong
basins and mouths, since T42 topography does not resolve real river networks. Climate
models route on prescribed networks instead (e.g. TRIP, MOSART, DDM30, HydroSHEDS-derived
maps). Here each Isca land cell sends its runoff to the mouth of the real basin that drains
most of its area, taken from the DDM30 30-minute drainage direction map (Doll and Lehner
2002, J. Hydrol. 258, 214-231; ISIMIP3 version on the CRU land mask, CC BY 4.0,
doi:10.48364/ISIMIP.865475), stored in input/rivers/. The mouth is then joined to the nearest
coupled ocean cell. There is no travel time: runoff reaches the mouth in the same step.

Cells the network does not cover (Antarctica, and land in Isca but not in DDM30) fall back to
downhill routing. Endorheic basins (inland DDM30 sinks, e.g. the Caspian's) have nowhere to
go in the coupled model, so their water goes to the coupled ocean nearest the sink.
"""
import os

import netCDF4
import numpy as np

# DDM30 flow directions: 1 E, 2 SE, 3 S, 4 SW, 5 W, 6 NW, 7 N, 8 NE; 0 outlet or sink
_STEP = {1: (0, 1), 2: (-1, 1), 3: (-1, 0), 4: (-1, -1), 5: (0, -1), 6: (1, -1), 7: (1, 0), 8: (1, 1)}
DDM30_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..',
                          'input', 'rivers', 'ddm30_flowdir_cru_neva.nc')


def gaussian_edges(nlat):
    """Latitudes [deg] of the edges of Isca's Gaussian grid cells, south to north."""
    _, w = np.polynomial.legendre.leggauss(nlat)
    return np.degrees(np.arcsin(np.clip(np.concatenate([[-1.], -1. + np.cumsum(w)]), -1, 1)))


def _xyz(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)


def ddm30_terminals(path=DDM30_FILE):
    """lat, lon (1-D, ascending lat) and, for every DDM30 land cell, the flat index of the
    cell its water ends in (an outlet to the sea or an inland sink); -1 off the network."""
    with netCDF4.Dataset(path) as nc:
        lat = np.asarray(nc['lat'][:], dtype=float)
        lon = np.asarray(nc['lon'][:], dtype=float)
        fd = np.ma.filled(nc['flowdirection'][:].astype(float), np.nan)
    if lat[0] > lat[-1]:
        lat, fd = lat[::-1], fd[::-1]
    nj, ni = fd.shape
    nxt = np.full(nj * ni, -1, dtype=np.int64)
    land = np.isfinite(fd)
    jj, ii = np.nonzero(land)
    for k, (dj, di) in _STEP.items():
        sel = fd[jj, ii] == k
        j2, i2 = jj[sel] + dj, (ii[sel] + di) % ni
        ok = (j2 >= 0) & (j2 < nj)
        src = jj[sel][ok] * ni + ii[sel][ok]
        dst = j2[ok] * ni + i2[ok]
        # a step into the sea or off the map ends the path at the current cell
        nxt[src] = np.where(land[j2[ok], i2[ok]], dst, src)
    terminal_cells = land.ravel() & (nxt < 0)       # 0 (outlet/sink) and -1 cells end paths
    nxt[terminal_cells] = np.nonzero(terminal_cells)[0]
    # pointer jumping: follow the network until every cell points at its end
    for _ in range(64):
        new = np.where(nxt >= 0, nxt[np.maximum(nxt, 0)], -1)
        if np.array_equal(new, nxt):
            break
        nxt = new
    return lat, lon, fd, nxt.reshape(nj, ni)


def isca_river_mouths(lat_edges, nlon, path=DDM30_FILE):
    """Mouth (lat, lon) of the basin draining most of each Isca cell's area (an outlet to the
    sea, or an inland sink), NaN where no DDM30 land lies in the cell."""
    lat, lon, fd, term = ddm30_terminals(path)
    nj, ni = fd.shape
    LON, LAT = np.meshgrid(lon % 360., lat)
    land = term >= 0
    nlat = len(lat_edges) - 1
    J = np.clip(np.searchsorted(lat_edges, LAT[land]) - 1, 0, nlat - 1)
    I = np.floor((LON[land] + 180. / nlon) / (360. / nlon)).astype(int) % nlon
    cell = J * nlon + I
    t = term[land]
    w = np.cos(np.radians(LAT[land]))
    mouth_lat = np.full((nlat, nlon), np.nan)
    mouth_lon = np.full((nlat, nlon), np.nan)
    order = np.lexsort((t, cell))
    cell, t, w = cell[order], t[order], w[order]
    keys = cell.astype(np.int64) * (nj * ni) + t
    uk, start = np.unique(keys, return_index=True)
    area = np.add.reduceat(w, start)
    c_of = uk // (nj * ni)
    t_of = uk % (nj * ni)
    best = {}
    for c, tt, a in zip(c_of, t_of, area):
        if c not in best or a > best[c][1]:
            best[c] = (tt, a)
    for c, (tt, _) in best.items():
        j, i = divmod(int(c), nlon)
        tj, ti = divmod(int(tt), ni)
        mouth_lat[j, i], mouth_lon[j, i] = lat[tj], lon[ti] % 360.
    return mouth_lat, mouth_lon


def basin_runoff_destinations(coupled, lat_edges, zsurf=None, path=DDM30_FILE):
    """Like isca.coupling.runoff_destinations, but land cells send their runoff to the coupled
    ocean cell nearest the mouth of their real (DDM30) basin. Returns the 1-based destination
    index i + (j-1)*nlon, 0 for coupled cells."""
    from isca.coupling import runoff_destinations, gaussian_latitudes
    nlat, nlon = coupled.shape
    dest = runoff_destinations(coupled, zsurf)          # fallback where DDM30 has no land
    mouth_lat, mouth_lon = isca_river_mouths(lat_edges, nlon, path)
    lat_c = gaussian_latitudes(nlat)
    lon_c = np.arange(nlon) * 360. / nlon
    LON, LAT = np.meshgrid(lon_c, lat_c)
    oj, oi = np.nonzero(coupled)
    ocean_xyz = _xyz(LON[oj, oi], LAT[oj, oi])
    use = ~coupled & np.isfinite(mouth_lat)
    k = np.argmax(_xyz(mouth_lon[use], mouth_lat[use]) @ ocean_xyz.T, axis=1)
    dest[use] = oi[k] + 1 + oj[k] * nlon
    return dest

"""Tripolar ocean grids (Murray 1996), so a MOM6 ocean coupled to Isca can reach the North
Pole, as GFDL's OM4 does.

South of lat_join the grid is latitude-longitude. North of it the cap is mapped
conformally onto bipolar coordinates in a stereographic projection from the North Pole:
with w = i a cot(z/2), z = sigma + i tau, the join circle |w| = a is sigma = +-pi/2,
the two foci w = +-a (on the join circle, at longitudes lon0 and lon0 + 180) are the
grid's two northern poles, and the segment between them through the North Pole is
sigma = +-pi, the fold. Each grid longitude on the join circle keeps its tau, and sigma
runs uniformly from the join circle to the fold. The grid is orthogonal everywhere, and
the fold maps column i onto column ni - i, as MOM6's TRIPOLAR_N expects. Cells around the
two foci shrink to nothing, so they must be land (as OM4 puts them on land): see
pole_island_mask.

Metric terms (dx, dy, area, angle_dx) of the 'supergrid' MOM6 reads are computed with
spherical geometry, so they are exact for the given points.
"""
import numpy as np
import netCDF4


def _unit_vectors(lon, lat):
    lam, phi = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], axis=-1)


def great_circle(lon1, lat1, lon2, lat2, radius):
    a, b = _unit_vectors(lon1, lat1), _unit_vectors(lon2, lat2)
    return radius * np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1))


def _triangle_area(a, b, c):
    """Area of the spherical triangle with unit-vector corners a, b, c (unit sphere)."""
    num = np.abs(np.sum(a * np.cross(b, c), axis=-1))
    den = 1. + np.sum(a * b, axis=-1) + np.sum(b * c, axis=-1) + np.sum(c * a, axis=-1)
    return 2. * np.arctan2(num, den)


def quad_areas(x, y, radius):
    """Areas of the cells of a (ny+1, nx+1) array of corners, as two spherical triangles."""
    p = _unit_vectors(x, y)
    sw, se, ne, nw = p[:-1, :-1], p[:-1, 1:], p[1:, 1:], p[1:, :-1]
    return radius**2 * (_triangle_area(sw, se, ne) + _triangle_area(sw, ne, nw))


def tripolar_points(lon_e, lat_e_south, lat_join, n_cap):
    """Supergrid points (x, y) [degrees], shape (2*nj+1, 2*ni+1).

    lon_e: the ni+1 cell-edge longitudes, uniform, starting at the first focus' longitude
        and symmetric about lon_e[0] + 180 (ni even);
    lat_e_south: cell-edge latitudes from the southern boundary up to, not including,
        lat_join; n_cap: number of cell rows from lat_join to the fold.
    """
    lon_e = np.asarray(lon_e, float)
    ni = len(lon_e) - 1
    if ni % 2:
        raise ValueError('a tripolar grid needs an even number of longitudes')
    xs = np.empty(2 * ni + 1); xs[0::2] = lon_e; xs[1::2] = 0.5 * (lon_e[:-1] + lon_e[1:])
    lat_s = np.append(np.asarray(lat_e_south, float), lat_join)
    ys_s = np.empty(2 * (len(lat_s) - 1) + 1); ys_s[0::2] = lat_s
    ys_s[1::2] = 0.5 * (lat_s[:-1] + lat_s[1:])
    xsouth, ysouth = np.meshgrid(xs, ys_s)

    # the cap, in the stereographic plane of the unit sphere: w = 2 tan(colat/2) e^(i lon)
    a = 2. * np.tan(np.radians(90. - lat_join) / 2.)
    lam = np.radians(xs - lon_e[0])           # longitude measured from the first focus
    with np.errstate(divide='ignore', invalid='ignore'):   # infinite at the foci
        z0 = 2. * np.arctan(1j * a / (a * np.exp(1j * lam)))   # z on the join circle
    tau = z0.imag
    s0 = np.sign(z0.real)                      # +pi/2 east of the focus line, -pi/2 west
    s0[s0 == 0] = 1.
    m = np.arange(1, 2 * n_cap + 1) / (2. * n_cap)          # supergrid rows above the join
    sigma = s0[None, :] * (np.pi / 2.) * (1. + m[:, None])
    with np.errstate(divide='ignore', invalid='ignore'):
        w = 1j * a / np.tan((sigma + 1j * tau[None, :]) / 2.)
    # the foci (tau infinite): every row of that column is the focus itself
    focus = ~np.isfinite(tau)
    w[:, focus] = a * np.sign(np.cos(lam[focus]))[None, :]
    r, theta = np.abs(w), np.angle(w)
    colat = 2. * np.arctan(r / 2.)
    xcap = np.degrees(theta)[...] + lon_e[0]
    ycap = 90. - np.degrees(colat)
    # keep longitudes continuous along each row (no jumps of 360)
    xcap = np.degrees(np.unwrap(np.radians(xcap), axis=1))
    xcap += 360. * np.round((xs[None, :] - xcap) / 360.)
    x = np.vstack([xsouth, xcap]); y = np.vstack([ysouth, ycap])
    return x, y


def tripolar_metrics(x, y, radius):
    """dx (nyp, nx), dy (ny, nxp), area (ny, nx) and angle_dx (nyp, nxp) of a supergrid."""
    dx = great_circle(x[:, :-1], y[:, :-1], x[:, 1:], y[:, 1:], radius)
    dy = great_circle(x[:-1, :], y[:-1, :], x[1:, :], y[1:, :], radius)
    area = quad_areas(x, y, radius)
    # angle of the grid's i direction anticlockwise from east, from centred differences
    xe = np.pad(x, ((0, 0), (1, 1)), mode='edge'); ye = np.pad(y, ((0, 0), (1, 1)), mode='edge')
    dlon = np.radians(xe[:, 2:] - xe[:, :-2]); dlon = (dlon + np.pi) % (2 * np.pi) - np.pi
    dlat = np.radians(ye[:, 2:] - ye[:, :-2])
    angle = np.degrees(np.arctan2(dlat, dlon * np.cos(np.radians(y))))
    return dx, dy, area, angle


def write_tripolar_supergrid(fname, lon_e, lat_e_south, lat_join, n_cap, radius):
    """Write MOM6's ocean_hgrid.nc for a tripolar grid; return (ni, nj)."""
    x, y = tripolar_points(lon_e, lat_e_south, lat_join, n_cap)
    dx, dy, area, angle = tripolar_metrics(x, y, radius)
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nyp', x.shape[0]); nc.createDimension('nxp', x.shape[1])
        nc.createDimension('ny', x.shape[0] - 1); nc.createDimension('nx', x.shape[1] - 1)
        nc.createDimension('string', 255)
        for name, dims, data in (('x', ('nyp', 'nxp'), x), ('y', ('nyp', 'nxp'), y),
                                 ('dx', ('nyp', 'nx'), dx), ('dy', ('ny', 'nxp'), dy),
                                 ('area', ('ny', 'nx'), area), ('angle_dx', ('nyp', 'nxp'), angle)):
            nc.createVariable(name, 'f8', dims)[:] = data
        chars = np.zeros(255, dtype='S1'); chars[:5] = list('tile1')
        nc.createVariable('tile', 'S1', ('string',))[:] = chars
    return x.shape[1] // 2, x.shape[0] // 2


def cell_centres(x, y):
    """Cell-centre longitudes and latitudes (the odd supergrid points)."""
    return x[1::2, 1::2], y[1::2, 1::2]


def pole_island_mask(x, y, radius, min_size):
    """True for cells (ny, nx) narrower than min_size [m] in either direction: the cells
    crowding around the two northern foci, which have to be land."""
    dx = great_circle(x[1::2, :-2:2], y[1::2, :-2:2], x[1::2, 2::2], y[1::2, 2::2], radius)
    dy = great_circle(x[:-2:2, 1::2], y[:-2:2, 1::2], x[2::2, 1::2], y[2::2, 1::2], radius)
    return (dx < min_size) | (dy < min_size)

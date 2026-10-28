"""MOM6 as the ocean component of a coupled Isca run.

MOM6 runs as its own program (mom6_isca.x, built from MOM6's FMS_cap and the
driver in isca_ocean_coupling/mom6_isca) in the 'ocean' subdirectory of the
run directory, because its FMS reads input.nml, diag_table etc. from its
working directory just as Isca's does.

The idealised configurations are generated here rather than read from files:
the horizontal grid is MOM6 'mosaic' supergrid whose cell edges are exactly
the Isca Gaussian-grid cell edges (so the conservative remapping between the
two is the identity wherever both exist), cut off at max_lat because MOM6
cannot include the poles; poleward of that Isca keeps its own slab ocean.
"""

import os
import shutil

import netCDF4
import numpy as np

from isca.coupling import OceanComponent, OCEAN_RESTART, OCEAN_TO_ISCA, ISCA_TO_OCEAN

ISCA_RADIUS = 6376.0e3     # constants_mod RADIUS
ISCA_CP_OCEAN = 3989.24495292815

# Parameters shared by every idealised configuration. The numerics follow
# MOM6-examples ocean_only/global_ALE (z* coordinate, ePBL, neutral diffusion);
# the lateral viscosity/diffusivities are chosen for coarse (2-3 degree) grids.
BASE_MOM_INPUT = {
    'REENTRANT_X': True,
    'TRIPOLAR_N': False,
    'GRID_CONFIG': '"mosaic"',
    'GRID_FILE': '"ocean_hgrid.nc"',
    'INPUTDIR': '"INPUT"',
    'RAD_EARTH': ISCA_RADIUS,
    'C_P': ISCA_CP_OCEAN,
    'RHO_0': 1035.0,
    'DT': 1800.0,
    'DT_THERM': 3600.0,
    'THICKNESSDIFFUSE': True,
    'THICKNESSDIFFUSE_FIRST': True,
    'DIABATIC_FIRST': True,
    'FRAZIL': True,
    'BOUND_SALINITY': True,
    'MIN_SALINITY': 0.01,
    'EQN_OF_STATE': '"WRIGHT"',
    'MAXIMUM_DEPTH': 4000.0,
    'MINIMUM_DEPTH': 0.5,
    'TOPO_CONFIG': '"flat"',
    'USE_REGRIDDING': True,
    'BULKMIXEDLAYER': False,
    'REGRIDDING_COORDINATE_MODE': '"Z*"',
    'ALE_COORDINATE_CONFIG': '"FILE:vgrid.nc,dz"',
    'REMAPPING_SCHEME': '"PPM_H4"',
    'COORD_CONFIG': '"ALE"',
    'INIT_LAYERS_FROM_Z_FILE': True,
    'TEMP_SALT_Z_INIT_FILE': '"ts_init.nc"',
    'Z_INIT_FILE_PTEMP_VAR': '"ptemp"',
    'Z_INIT_FILE_SALT_VAR': '"salt"',
    'Z_INIT_ALE_REMAPPING': True,
    'ENERGETICS_SFC_PBL': True,
    'USE_MLD_ITERATION': False,
    'HMIX_FIXED': 0.5,
    'MIXEDLAYER_RESTRAT': True,
    'FOX_KEMPER_ML_RESTRAT_COEF': 20.0,
    'USE_NEUTRAL_DIFFUSION': True,
    'USE_STORED_SLOPES': True,
    'KHTH': 800.0,
    'KHTR': 800.0,
    'KD': 2.0e-5,
    'KV': 1.0e-4,
    'HBBL': 10.0,
    'LAPLACIAN': True,
    'KH': 2.0e4,
    'KH_VEL_SCALE': 0.01,
    'BIHARMONIC': True,
    'SMAGORINSKY_AH': True,
    'SMAG_BI_CONST': 0.06,
    'AH_VEL_SCALE': 0.05,
    'MAXVEL': 6.0,
    'BOUND_CORIOLIS': True,
    'DTBT': -0.95,
    'BEBT': 0.2,
    'VAR_PEN_SW': False,
    'PEN_SW_SCALE': 15.0,
    'PEN_SW_FRAC': 0.42,
    'TRACER_ADVECTION_SCHEME': '"PPM:H3"',
    # off: Isca routes all surface water to the ocean (runoff), so the budget closes by itself
    'ADJUST_NET_FRESH_WATER_TO_ZERO': False,
    'CHECK_BAD_SURFACE_VALS': True,
    'PARALLEL_RESTARTFILES': False,
    'ENERGYSAVEDAYS': 1.0,
    'MAXTRUNC': 5000,
}


def gaussian_grid(nlat):
    """Isca's Gaussian latitudes and cell edges [degrees], south to north."""
    sin_lat, wts = np.polynomial.legendre.leggauss(nlat)
    sin_edges = np.concatenate([[-1.], -1. + np.cumsum(wts)])
    sin_edges[-1] = 1.
    return np.degrees(np.arcsin(sin_lat)), np.degrees(np.arcsin(np.clip(sin_edges, -1, 1)))


def ocean_rows(nlat, max_lat):
    """Indices of the Isca latitude rows the ocean covers: cells entirely within |lat| <= max_lat."""
    _, edges = gaussian_grid(nlat)
    return np.where((edges[:-1] >= -max_lat) & (edges[1:] <= max_lat))[0]


def write_supergrid(fname, nlon, nlat, max_lat, radius=ISCA_RADIUS):
    """MOM6 'mosaic' supergrid with the Isca Gaussian grid's cell edges and centres."""
    lat_c, lat_e = gaussian_grid(nlat)
    rows = ocean_rows(nlat, max_lat)
    # supergrid latitudes: edge, centre, edge, centre, ..., edge
    ys = np.empty(2 * len(rows) + 1)
    ys[0::2] = lat_e[rows[0]:rows[-1] + 2]
    ys[1::2] = lat_c[rows]
    # Isca's first cell is centred on 0E, so its edges start at -dlon/2
    xs = np.linspace(0., 360., 2 * nlon + 1) - 180. / nlon
    x, y = np.meshgrid(xs, ys)
    lam, phi = np.radians(xs), np.radians(ys)
    dx = radius * np.cos(np.radians(y[:, :-1])) * np.diff(lam)[None, :]
    dy = np.broadcast_to((radius * np.diff(phi))[:, None], (len(ys) - 1, len(xs))).copy()
    area = radius**2 * np.diff(lam)[None, :] * np.diff(np.sin(phi))[:, None]
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nyp', len(ys)); nc.createDimension('nxp', len(xs))
        nc.createDimension('ny', len(ys) - 1); nc.createDimension('nx', len(xs) - 1)
        nc.createDimension('string', 255)
        for name, dims, data, units in (('x', ('nyp', 'nxp'), x, 'degree_east'),
                                         ('y', ('nyp', 'nxp'), y, 'degree_north'),
                                         ('dx', ('nyp', 'nx'), dx, 'm'),
                                         ('dy', ('ny', 'nxp'), dy, 'm'),
                                         ('area', ('ny', 'nx'), area, 'm2'),
                                         ('angle_dx', ('nyp', 'nxp'), np.zeros_like(x), 'degrees')):
            v = nc.createVariable(name, 'f8', dims)
            v.units = units
            v[:] = data
        tile = nc.createVariable('tile', 'S1', ('string',))
        chars = np.zeros(255, dtype='S1')
        chars[:5] = list('tile1')
        tile[:] = chars
    return len(rows)


def regular_grid(dlat, max_lat):
    """Cell edges [degrees] of a regular grid covering -max_lat..max_lat."""
    n = int(round(2 * max_lat / dlat))
    return np.linspace(-max_lat, max_lat, n + 1)


def write_supergrid_regular(fname, dlon, dlat, max_lat, radius=ISCA_RADIUS):
    """MOM6 supergrid for a regular dlon x dlat grid (cell edges from 0E): used to test
    coupling through a non-trivial (conservative) remapping."""
    lat_e = regular_grid(dlat, max_lat)
    ys = np.empty(2 * (len(lat_e) - 1) + 1)
    ys[0::2] = lat_e
    ys[1::2] = 0.5 * (lat_e[:-1] + lat_e[1:])
    nlon = int(round(360. / dlon))
    xs = np.linspace(0., 360., 2 * nlon + 1)
    x, y = np.meshgrid(xs, ys)
    lam, phi = np.radians(xs), np.radians(ys)
    dx = radius * np.cos(np.radians(y[:, :-1])) * np.diff(lam)[None, :]
    dy = np.broadcast_to((radius * np.diff(phi))[:, None], (len(ys) - 1, len(xs))).copy()
    area = radius**2 * np.diff(lam)[None, :] * np.diff(np.sin(phi))[:, None]
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nyp', len(ys)); nc.createDimension('nxp', len(xs))
        nc.createDimension('ny', len(ys) - 1); nc.createDimension('nx', len(xs) - 1)
        nc.createDimension('string', 255)
        for name, dims, data in (('x', ('nyp', 'nxp'), x), ('y', ('nyp', 'nxp'), y),
                                 ('dx', ('nyp', 'nx'), dx), ('dy', ('ny', 'nxp'), dy),
                                 ('area', ('ny', 'nx'), area),
                                 ('angle_dx', ('nyp', 'nxp'), np.zeros_like(x))):
            nc.createVariable(name, 'f8', dims)[:] = data
        chars = np.zeros(255, dtype='S1')
        chars[:5] = list('tile1')
        nc.createVariable('tile', 'S1', ('string',))[:] = chars
    return nlon, len(lat_e) - 1


def write_supergrid_edges(fname, lon_e, lat_e, lat_c=None, radius=ISCA_RADIUS):
    """MOM6 supergrid from arbitrary cell edges [degrees] (and optional centre latitudes)."""
    lat_c = 0.5 * (lat_e[:-1] + lat_e[1:]) if lat_c is None else lat_c
    ys = np.empty(2 * len(lat_c) + 1); ys[0::2] = lat_e; ys[1::2] = lat_c
    xs = np.empty(2 * (len(lon_e) - 1) + 1); xs[0::2] = lon_e; xs[1::2] = 0.5 * (lon_e[:-1] + lon_e[1:])
    x, y = np.meshgrid(xs, ys)
    lam, phi = np.radians(xs), np.radians(ys)
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nyp', len(ys)); nc.createDimension('nxp', len(xs))
        nc.createDimension('ny', len(ys) - 1); nc.createDimension('nx', len(xs) - 1)
        nc.createDimension('string', 255)
        for name, dims, data in (
                ('x', ('nyp', 'nxp'), x), ('y', ('nyp', 'nxp'), y),
                ('dx', ('nyp', 'nx'), radius * np.cos(np.radians(y[:, :-1])) * np.diff(lam)[None, :]),
                ('dy', ('ny', 'nxp'), np.broadcast_to((radius * np.diff(phi))[:, None], (len(ys) - 1, len(xs)))),
                ('area', ('ny', 'nx'), radius**2 * np.diff(lam)[None, :] * np.diff(np.sin(phi))[:, None]),
                ('angle_dx', ('nyp', 'nxp'), np.zeros_like(x))):
            nc.createVariable(name, 'f8', dims)[:] = data
        chars = np.zeros(255, dtype='S1'); chars[:5] = list('tile1')
        nc.createVariable('tile', 'S1', ('string',))[:] = chars
    return len(xs) // 2, len(ys) // 2


def write_vgrid(fname, nk=30, depth=4000.0, dz_top=10.0):
    """Layer thicknesses growing geometrically from dz_top to fill depth."""
    # solve for the stretching factor r: dz_top*(r**nk - 1)/(r - 1) = depth
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        if dz_top * (r**nk - 1) / (r - 1) > depth:
            hi = r
        else:
            lo = r
    dz = dz_top * r**np.arange(nk)
    dz *= depth / dz.sum()
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nz', nk)
        v = nc.createVariable('dz', 'f8', ('nz',))
        v.units = 'm'
        v[:] = dz
    return dz


def write_initial_ts(fname, nlon, nlat, max_lat, sst_func, t_deep=2.0, decay=800.0, salt=35.0):
    """Initial potential temperature/salinity on z levels for INIT_LAYERS_FROM_Z_FILE.

    The surface temperature is sst_func(lat) [K], decaying exponentially with
    depth to t_deep [degC] (never below freezing).
    """
    # global 1-degree data, so it covers whatever ocean grid MOM6 is using
    # (FMS horizontal interpolation fails if the data do not cover the model grid)
    lat = np.arange(-89.5, 90., 1.)
    lon = np.arange(0.5, 360., 1.)
    nlon = len(lon)
    z = np.array([0., 5., 15., 30., 50., 75., 100., 150., 200., 300., 400., 500., 700., 1000.,
                  1500., 2000., 2500., 3000., 3500., 4000., 4500.])
    sst = np.maximum(sst_func(lat) - 273.15, -1.8)
    temp = t_deep + (sst[None, :] - t_deep) * np.exp(-z[:, None] / decay)
    temp = np.maximum(temp, -1.8)
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('lon', nlon); nc.createDimension('lat', len(lat)); nc.createDimension('depth', len(z))
        for name, data, units, axis in (('lon', lon, 'degrees_east', 'X'), ('lat', lat, 'degrees_north', 'Y'),
                                        ('depth', z, 'm', 'Z')):
            v = nc.createVariable(name, 'f8', (name,))
            v.units = units; v.axis = axis
            if axis == 'Z':
                v.positive = 'down'
            v[:] = data
        # MOM6's horizontal interpolation of initial conditions needs a missing value
        t = nc.createVariable('ptemp', 'f8', ('depth', 'lat', 'lon'), fill_value=-1.e20)
        t.missing_value = -1.e20
        t.units = 'degC'
        t[:] = np.broadcast_to(temp[:, :, None], (len(z), len(lat), nlon))
        s = nc.createVariable('salt', 'f8', ('depth', 'lat', 'lon'), fill_value=-1.e20)
        s.missing_value = -1.e20
        s.units = 'psu'
        s[:] = salt


def write_topography(fname, nlon, nlat, max_lat, land_mask, depth):
    """MOM6 TOPO_FILE: 'depth' [m, positive down] on the ocean T-points, 0 on land.

    land_mask is (nlat, nlon) on the Isca grid; the MOM6 grid is the Isca grid's
    ocean rows, so the two masks agree exactly."""
    rows = ocean_rows(nlat, max_lat)
    d = np.where(np.asarray(land_mask, dtype=bool)[rows, :], 0., depth)
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('ny', len(rows)); nc.createDimension('nx', nlon)
        v = nc.createVariable('depth', 'f8', ('ny', 'nx'))
        v.units = 'm'
        v[:] = d


def fill_dead_end_cells(depth):
    """Set to land (repeatedly) wet cells with at most one wet neighbour (zonally periodic):
    single-cell inlets have no through-flow and blow up (e.g. the Gulf of Ob at 2 degrees)."""
    depth = np.array(depth, dtype=float)
    while True:
        w = (depth > 0).astype(int)
        n = np.roll(w, 1, 1) + np.roll(w, -1, 1)
        n[1:] += w[:-1]; n[:-1] += w[1:]
        dead = (w > 0) & (n <= 1)
        if not dead.any():
            return depth
        depth[dead] = 0.


# Straits closed by regridding OM4 to 2 degrees, opened as OM4 does by hand:
# (lon, lat) of the cell centres to open and their sill depth [m].
EARTH_STRAITS_2DEG = (
    [(355., 36.8, 300.), (357., 36.8, 300.), (359., 36.8, 300.)]                  # Gibraltar
    + [(lo, 56.7, 30.) for lo in (9., 11., 13., 15.)]                              # Kattegat/Danish straits
    + [(29., 42.8, 50.), (29., 40.8, 50.), (27., 40.8, 50.), (27., 38.8, 50.)]    # Bosporus/Dardanelles
    + [(lo, 62.7, 300.) for lo in range(283, 295, 2)])                            # Hudson Strait


def connect_basins(depth, max_channel=4, max_sill=300., min_sill=20.):
    """Join every enclosed basin to the largest one through the shortest channel of at most
    max_channel land cells (4-connected, zonally periodic), as OM4 opens its straits by hand;
    the channel depth is the basin's mean depth within [min_sill, max_sill]. Basins with no
    such channel are filled. Coarse regridding closes Gibraltar, the Danish straits, the
    Bosporus and Hudson Strait, and closed basins fill or empty with P-E and runoff."""
    from collections import deque
    from scipy import ndimage
    depth = np.array(depth, dtype=float)
    nj, ni = depth.shape
    while True:
        lab, _ = ndimage.label(depth > 0)
        for j in range(nj):                       # zonal periodicity
            a, b = lab[j, 0], lab[j, -1]
            if a and b and a != b:
                lab[lab == b] = a
        ids, counts = np.unique(lab[lab > 0], return_counts=True)
        if len(ids) == 1:
            return depth
        main = ids[np.argmax(counts)]
        basin = ids[ids != main][np.argmin(counts[ids != main])]   # smallest first
        inb = lab == basin
        # breadth-first search through land from the basin to any other basin
        prev = {(j, i): None for j, i in zip(*np.where(inb))}
        queue, end = deque(prev), None
        while queue and end is None:
            j, i = queue.popleft()
            steps, p = 0, prev[(j, i)]
            while p is not None and not inb[p]:
                steps, p = steps + 1, prev[p]
            if steps >= max_channel:
                continue
            for jj, ii in ((j + 1, i), (j - 1, i), (j, (i + 1) % ni), (j, (i - 1) % ni)):
                if 0 <= jj < nj and (jj, ii) not in prev:
                    prev[(jj, ii)] = (j, i)
                    if lab[jj, ii] not in (0, basin):
                        end = (jj, ii)
                        break
                    queue.append((jj, ii))
        if end is None:
            depth[inb] = 0.
            continue
        sill = min(max(depth[inb].mean(), min_sill), max_sill)
        p = prev[end]
        while not inb[p]:
            depth[p] = sill
            p = prev[p]


def write_isca_land_file(fname, land_mask, zsurf=None):
    """Isca land file (idealized_moist_phys_nml land_option='input')."""
    nlat, nlon = land_mask.shape
    lat_c, _ = gaussian_grid(nlat)
    with netCDF4.Dataset(fname, 'w', format='NETCDF3_CLASSIC') as nc:
        nc.createDimension('lat', nlat); nc.createDimension('lon', nlon)
        nc.createVariable('lat', 'f4', ('lat',))[:] = lat_c
        nc.createVariable('lon', 'f4', ('lon',))[:] = np.arange(nlon) * 360. / nlon
        nc.createVariable('land_mask', 'f4', ('lat', 'lon'))[:] = np.asarray(land_mask, dtype=float)
        nc.createVariable('zsurf', 'f4', ('lat', 'lon'))[:] = 0. if zsurf is None else zsurf


def format_mom_input(params):
    lines = []
    for k, v in params.items():
        if isinstance(v, bool):
            v = 'True' if v else 'False'
        lines.append('%s = %s' % (k, v))
    return '\n'.join(lines) + '\n'


DIAG_TABLE = """Isca coupled MOM6
1 1 1 0 0 0
"ocean_month", 30, "days", 1, "days", "time"
"ocean_month_z", 30, "days", 1, "days", "time"
"ocean_static", -1, "days", 1, "days", "time"
"ocean_model", "geolon", "geolon", "ocean_static", "all", .false., "none", 2
"ocean_model", "geolat", "geolat", "ocean_static", "all", .false., "none", 2
"ocean_model", "area_t", "area_t", "ocean_static", "all", .false., "none", 2
"ocean_model", "wet", "wet", "ocean_static", "all", .false., "none", 2
"ocean_model", "tos", "tos", "ocean_month", "all", .true., "none", 2
"ocean_model", "sos", "sos", "ocean_month", "all", .true., "none", 2
"ocean_model", "zos", "zos", "ocean_month", "all", .true., "none", 2
"ocean_model", "mlotst", "mlotst", "ocean_month", "all", .true., "none", 2
"ocean_model", "hfds", "hfds", "ocean_month", "all", .true., "none", 2
"ocean_model", "taux", "taux", "ocean_month", "all", .true., "none", 2
"ocean_model", "tauy", "tauy", "ocean_month", "all", .true., "none", 2
"ocean_model", "frazil_heat_tend_2d", "frazil_heat_tend_2d", "ocean_month", "all", .true., "none", 2
"ocean_model", "net_heat_surface", "net_heat_surface", "ocean_month", "all", .true., "none", 2
"ocean_model", "PRCmE", "PRCmE", "ocean_month", "all", .true., "none", 2
"ocean_model", "friver", "friver", "ocean_month", "all", .true., "none", 2
"ocean_model", "net_heat_coupler", "net_heat_coupler", "ocean_month", "all", .true., "none", 2
"ocean_model", "heat_content_massout", "heat_content_massout", "ocean_month", "all", .true., "none", 2
"ocean_model", "heat_content_massin", "heat_content_massin", "ocean_month", "all", .true., "none", 2
"ocean_model", "evs", "evs", "ocean_month", "all", .true., "none", 2
"ocean_model", "prlq", "prlq", "ocean_month", "all", .true., "none", 2
"ocean_model", "T_ady_2d", "T_ady_2d", "ocean_month", "all", .true., "none", 2
"ocean_model", "T_diffy_2d", "T_diffy_2d", "ocean_month", "all", .true., "none", 2
"ocean_model", "vmo", "vmo", "ocean_month", "all", .true., "none", 2
"ocean_model", "uo", "uo", "ocean_month", "all", .true., "none", 2
"ocean_model", "vo", "vo", "ocean_month", "all", .true., "none", 2
"ocean_model", "thetao", "thetao", "ocean_month", "all", .true., "none", 2
"ocean_model", "so", "so", "ocean_month", "all", .true., "none", 2
"ocean_model_z", "thetao", "thetao", "ocean_month_z", "all", .true., "none", 2
"ocean_model_z", "vmo", "vmo", "ocean_month_z", "all", .true., "none", 2
"""


class MOM6Ocean(OceanComponent):
    """MOM6 on an idealised (aquaplanet-type) domain matched to the Isca grid."""
    grid_name = 'mom6'
    # the frazil heat MOM6 adds to keep water at freezing: Isca's thermodynamic sea ice
    # (mixed_layer_nml do_thermo_ice) grows ice from it instead
    fields_to_isca = OCEAN_TO_ISCA + [('FRAZIL_OCN', 'OCN_FRAZIL')]

    def __init__(self, executable, num_cores=4, max_lat=80., nk=30, depth=4000.,
                 layout=None, mom_input=None, sst_func=None, land_mask=None, regular_dlat=None):
        """land_mask: optional (nlat, nlon) boolean array on the Isca grid (e.g. a
        ridge); MOM6 gets zero depth there. Give the same mask to Isca
        (CoupledExperiment land_mask, and an Isca land file) so they agree."""
        self.land_mask = land_mask
        # None: MOM6 cells are exactly Isca's cells. A number: a regular grid of that
        # spacing instead (tests real remapping; land masks not supported then).
        self.regular_dlat = regular_dlat
        if regular_dlat is not None and land_mask is not None:
            raise ValueError('land_mask is only supported on the Isca-matched grid')
        self.executable = executable
        self.num_cores = num_cores
        self.max_lat = max_lat
        self.nk = nk
        self.depth = depth
        self.layout = layout
        self.mom_input = dict(BASE_MOM_INPUT)
        self.mom_input['MAXIMUM_DEPTH'] = depth
        self.mom_input.update(mom_input or {})
        # default initial SST: Isca's prescribe_initial_dist with tconst=285, delta_T=40
        self.sst_func = sst_func or (lambda lat: 285. - 40. * (3. * np.sin(np.radians(lat))**2 - 1.) / 3.)

    def ocean_grid_shape(self, nlon, nlat):
        if self.regular_dlat is not None:
            return int(round(360. / self.regular_dlat)), len(regular_grid(self.regular_dlat, self.max_lat)) - 1
        return nlon, len(ocean_rows(nlat, self.max_lat))

    def isca_ocean_mask(self, nlon, nlat):
        mask = np.zeros((nlat, nlon), dtype=bool)
        if self.regular_dlat is not None:
            # Isca cells more than half covered by the ocean (as Isca decides at run time)
            _, edges = gaussian_grid(nlat)
            s = np.sin(np.radians(edges))
            m = np.sin(np.radians(self.max_lat))
            covered = (np.minimum(s[1:], m) - np.maximum(s[:-1], -m)).clip(0) / np.diff(s)
            mask[covered > 0.5, :] = True
            return mask
        mask[ocean_rows(nlat, self.max_lat), :] = True
        return mask

    def mpmd_options(self, exp):
        return '-wdir ocean'

    def command(self, exp, runtime):
        return [self.executable]

    def write_grid_files(self, indir, nlon, nlat, params):
        """Write the horizontal grid, topography and initial T/S into indir; update
        params accordingly; return the ocean grid size (ni, nj)."""
        if self.regular_dlat is not None:
            ni, nj = write_supergrid_regular(os.path.join(indir, 'ocean_hgrid.nc'), self.regular_dlat,
                                             self.regular_dlat, self.max_lat)
        else:
            ni = nlon
            nj = write_supergrid(os.path.join(indir, 'ocean_hgrid.nc'), nlon, nlat, self.max_lat)
        write_initial_ts(os.path.join(indir, 'ts_init.nc'), nlon, nlat, self.max_lat, self.sst_func)
        if self.land_mask is not None:
            write_topography(os.path.join(indir, 'topog.nc'), nlon, nlat, self.max_lat,
                             self.land_mask, self.depth)
            params.update({'TOPO_CONFIG': '"file"', 'TOPO_FILE': '"topog.nc"', 'TOPO_VARNAME': '"depth"'})
        return ni, nj

    def ocean_dir(self, exp):
        return os.path.join(exp.rundir, 'ocean')

    def prepare_rundir(self, exp, cold_start):
        nlon, nlat = exp.isca_grid_shape()
        odir = self.ocean_dir(exp)
        indir = os.path.join(odir, 'INPUT')
        os.makedirs(indir, exist_ok=True)
        os.makedirs(os.path.join(odir, 'RESTART'), exist_ok=True)

        write_vgrid(os.path.join(indir, 'vgrid.nc'), self.nk, self.depth)
        params = dict(self.mom_input)
        ni, nj = self.write_grid_files(indir, nlon, nlat, params)
        params.update({'NIGLOBAL': ni, 'NJGLOBAL': nj, 'NK': self.nk})
        if self.layout:
            params['LAYOUT'] = '%d, %d' % tuple(self.layout)
        with open(os.path.join(odir, 'MOM_input'), 'w') as f:
            f.write(format_mom_input(params))
        open(os.path.join(odir, 'MOM_override'), 'w').close()
        with open(os.path.join(odir, 'diag_table'), 'w') as f:
            f.write(DIAG_TABLE)

        main = exp.namelist['main_nml']
        m = {k: int(main.get(k, 0)) for k in ('days', 'hours', 'minutes', 'seconds')}
        date = list(main.get('current_date', [1, 1, 1, 0, 0, 0]))
        with open(os.path.join(odir, 'input.nml'), 'w') as f:
            f.write("&mom6_isca_nml\n  calendar = '%s'\n  current_date = %s\n"
                    "  days = %d, hours = %d, minutes = %d, seconds = %d\n  dt_cpld = %d\n"
                    "  oasis_restart_file = '%s'\n/\n"
                    % (main.get('calendar', 'thirty_day').lower(), ', '.join(str(d) for d in date),
                       m['days'], m['hours'], m['minutes'], m['seconds'], exp.coupling_period,
                       OCEAN_RESTART))
            f.write("&MOM_input_nml\n  output_directory = './'\n  input_filename = '%s'\n"
                    "  restart_input_dir = 'INPUT/'\n  restart_output_dir = 'RESTART/'\n"
                    "  parameter_filename = 'MOM_input', 'MOM_override'\n/\n"
                    % ('n' if cold_start else 'r'))
            f.write("&fms_nml\n  domains_stack_size = 2000000\n/\n")
            f.write("&diag_manager_nml\n/\n")

        # OASIS files are shared with Isca's directory (OASIS writes the grid files one
        # model at a time); the restart file the ocean sends from is moved here too.
        for name in ('namcouple', 'grids.nc', 'masks.nc', 'areas.nc'):
            link = os.path.join(odir, name)
            if not os.path.lexists(link):
                os.symlink(os.path.join('..', name), link)
        rst = os.path.join(exp.rundir, OCEAN_RESTART)
        if os.path.exists(rst):
            shutil.move(rst, os.path.join(odir, OCEAN_RESTART))

    def save_restart(self, exp, resdir):
        odir = self.ocean_dir(exp)
        dest = os.path.join(resdir, 'ocean')
        os.makedirs(dest, exist_ok=True)
        for f in os.listdir(os.path.join(odir, 'RESTART')):
            shutil.copy(os.path.join(odir, 'RESTART', f), dest)
        rst = os.path.join(odir, OCEAN_RESTART)
        if os.path.exists(rst):
            shutil.copy(rst, os.path.join(resdir, OCEAN_RESTART))

    def save_output(self, exp, outdir):
        odir = self.ocean_dir(exp)
        dest = os.path.join(outdir, 'ocean')
        os.makedirs(dest, exist_ok=True)
        for f in sorted(os.listdir(odir)):
            if (f.startswith('ocean_') and '.nc' in f) or f in ('ocean.stats', 'ocean.stats.nc',
                                                                 'MOM_parameter_doc.all', 'input.nml',
                                                                 'MOM_input', 'MOM_override'):
                shutil.copy(os.path.join(odir, f), dest)

    def load_restart(self, exp, indir):
        src = os.path.join(indir, 'ocean')
        dest = os.path.join(self.ocean_dir(exp), 'INPUT')
        os.makedirs(dest, exist_ok=True)
        if os.path.isdir(src):
            for f in os.listdir(src):
                shutil.copy(os.path.join(src, f), dest)


class MOM6TripolarOcean(MOM6Ocean):
    """MOM6 on a tripolar grid (see tripolar.py) that reaches the North Pole, for aquaplanet
    and ridge configurations without northern polar caps. South of lat_join the rows are
    Isca's Gaussian rows from lat_south (the ocean's southern wall); north of it, n_cap rows
    of the bipolar cap, with the two northern poles at 0E and 180E on lat_join. Cells
    narrower than island_size around those poles are land (two small islands, as OM4 puts
    its poles on land). land_func(lon, lat) -> bool, if given, adds land (e.g. a ridge)
    by cell centre. Isca and this grid never match, so use Isca's surface tiles."""

    def __init__(self, executable, lat_south=-85., lat_join=65., n_cap=26, island_size=20.e3,
                 land_func=None, **kwargs):
        super(MOM6TripolarOcean, self).__init__(executable, **kwargs)
        self.lat_south, self.lat_join, self.n_cap = lat_south, lat_join, n_cap
        self.island_size, self.land_func = island_size, land_func
        self.mom_input['TRIPOLAR_N'] = True

    def _edges(self, nlon, nlat):
        _, e = gaussian_grid(nlat)
        join = e[np.argmin(abs(e - self.lat_join))]       # the Gaussian edge nearest lat_join
        south = e[(e >= self.lat_south) & (e < join)]
        return np.linspace(0., 360., nlon + 1), south, join

    def ocean_grid_shape(self, nlon, nlat):
        _, south, _ = self._edges(nlon, nlat)
        return nlon, len(south) + self.n_cap

    def isca_ocean_mask(self, nlon, nlat):
        _, e = gaussian_grid(nlat)
        mask = np.zeros((nlat, nlon), dtype=bool)
        mask[e[:-1] >= self.lat_south, :] = True
        return mask

    def write_grid_files(self, indir, nlon, nlat, params):
        from isca.coupling.tripolar import write_tripolar_supergrid, tripolar_points, pole_island_mask
        lon_e, south, join = self._edges(nlon, nlat)
        ni, nj = write_tripolar_supergrid(os.path.join(indir, 'ocean_hgrid.nc'), lon_e, south, join,
                                          self.n_cap, ISCA_RADIUS)
        x, y = tripolar_points(lon_e, south, join, self.n_cap)
        land = pole_island_mask(x, y, ISCA_RADIUS, self.island_size)
        if self.land_func is not None:
            land |= np.asarray(self.land_func(x[1::2, 1::2] % 360., y[1::2, 1::2]), dtype=bool)
        with netCDF4.Dataset(os.path.join(indir, 'topog.nc'), 'w') as nc:
            nc.createDimension('ny', nj); nc.createDimension('nx', ni)
            v = nc.createVariable('depth', 'f8', ('ny', 'nx')); v.units = 'm'
            v[:] = np.where(land, 0., self.depth)
        params.update({'TOPO_CONFIG': '"file"', 'TOPO_FILE': '"topog.nc"', 'TOPO_VARNAME': '"depth"'})
        write_initial_ts(os.path.join(indir, 'ts_init.nc'), nlon, nlat, 90., self.sst_func)
        return ni, nj


class MOM6EarthOcean(MOM6Ocean):
    """MOM6 with Earth geometry derived from GFDL OM4 data (as regridded by M2LiNES for
    Samudra): bathymetry from the per-level ocean fractions, initial temperature and
    salinity from an OM4 5-day mean. The grid is the OM4-derived lat-lon grid itself,
    equatorward of max_lat. Sea ice is not modelled (use Isca's prescribed ice lid).

    om4_zarr: path/URL of an om4_*deg OM4.zarr; init_date: date of the initial state.
    """

    def __init__(self, executable, om4_zarr, init_date='2014-10-10', max_lat=78., **kwargs):
        super(MOM6EarthOcean, self).__init__(executable, max_lat=max_lat, depth=6500., **kwargs)
        self.om4_zarr, self.init_date = om4_zarr, init_date
        self._grid = None

    def _om4(self):
        if self._grid is None:
            import xarray as xr
            ds = xr.open_zarr(self.om4_zarr)
            lat_b = np.asarray(ds.lat_b.values)[:, 0]
            lon_b = np.asarray(ds.lon_b.values)[0, :]
            lat_c = np.asarray(ds.y.values)
            rows = np.where((lat_b[:-1] >= -self.max_lat) & (lat_b[1:] <= self.max_lat))[0]
            frac = ds.ocean_fraction.isel(y=rows).values           # (lev, ny, nx)
            dz = np.asarray(ds.dz.values)
            depth = (frac * dz[:, None, None]).sum(0)
            depth = np.where(frac[0] >= 0.5, np.maximum(depth, 10.), 0.)
            depth = fill_dead_end_cells(depth)
            lon_c = 0.5 * (lon_b[:-1] + lon_b[1:]) % 360.
            for lo, la, sill in EARTH_STRAITS_2DEG:
                j = np.argmin(abs(lat_c[rows] - la)); i = np.argmin(abs(lon_c - lo))
                depth[j, i] = max(depth[j, i], sill)
            depth = fill_dead_end_cells(connect_basins(depth, max_channel=2))
            self._grid = dict(ds=ds, rows=rows, lat_e=lat_b[rows[0]:rows[-1] + 2], lat_c=lat_c[rows],
                              lon_e=lon_b, depth=depth)
        return self._grid

    def ocean_grid_shape(self, nlon, nlat):
        g = self._om4()
        return len(g['lon_e']) - 1, len(g['lat_c'])

    def isca_ocean_mask(self, nlon, nlat):
        # Isca cells more than half inside max_lat (the land mask is removed separately)
        _, edges = gaussian_grid(nlat)
        s = np.sin(np.radians(edges)); m = np.sin(np.radians(self.max_lat))
        covered = (np.minimum(s[1:], m) - np.maximum(s[:-1], -m)).clip(0) / np.diff(s)
        return np.broadcast_to((covered > 0.5)[:, None], (nlat, nlon)).copy()

    def write_grid_files(self, indir, nlon, nlat, params):
        # generated once from the (remote) OM4 data, then reused for every run
        tag = '%s_%s_%g_v3' % (os.path.basename(os.path.dirname(self.om4_zarr.rstrip('/'))),
                            self.init_date, self.max_lat)
        cache = os.path.join(os.environ.get('GFDL_WORK', '.'), 'mom6_earth_cache', tag)
        files = ('ocean_hgrid.nc', 'topog.nc', 'ts_init.nc')
        if not all(os.path.exists(os.path.join(cache, f)) for f in files):
            os.makedirs(cache, exist_ok=True)
            self._write_grid_files(cache)
        for f in files:
            shutil.copy(os.path.join(cache, f), indir)
        params.update({'TOPO_CONFIG': '"file"', 'TOPO_FILE': '"topog.nc"', 'TOPO_VARNAME': '"depth"',
                       'MINIMUM_DEPTH': 9.5, 'MASKING_DEPTH': 9.0})
        with netCDF4.Dataset(os.path.join(indir, 'topog.nc')) as nc:
            nj, ni = nc['depth'].shape
        return ni, nj

    def _write_grid_files(self, indir):
        g = self._om4()
        ni, nj = write_supergrid_edges(os.path.join(indir, 'ocean_hgrid.nc'), g['lon_e'], g['lat_e'], g['lat_c'])
        with netCDF4.Dataset(os.path.join(indir, 'topog.nc'), 'w') as nc:
            nc.createDimension('ny', nj); nc.createDimension('nx', ni)
            nc.createVariable('depth', 'f8', ('ny', 'nx'))[:] = g['depth']
        # initial state: OM4 thetao/so on its 19 levels (land filled by MOM6's extrapolation)
        ds = g['ds'].sel(time=slice(self.init_date, None)).isel(time=0)   # first record from init_date (cftime axis)
        lev = np.asarray(ds.lev.values)
        # float64 first: a float32 -1e20 is not equal to the float64 missing_value
        T = np.stack([ds['thetao_%d' % k].values for k in range(len(lev))]).astype(np.float64)
        S = np.stack([ds['so_%d' % k].values for k in range(len(lev))]).astype(np.float64)
        lat = np.asarray(g['ds'].y.values); lon = 0.5 * (g['lon_e'][:-1] + g['lon_e'][1:])
        with netCDF4.Dataset(os.path.join(indir, 'ts_init.nc'), 'w') as nc:
            nc.createDimension('lon', len(lon)); nc.createDimension('lat', len(lat)); nc.createDimension('depth', len(lev))
            for name, data, axis in (('lon', lon, 'X'), ('lat', lat, 'Y'), ('depth', lev, 'Z')):
                v = nc.createVariable(name, 'f8', (name,)); v.axis = axis; v[:] = data
                if axis == 'Z':
                    v.positive = 'down'
            for name, data in (('ptemp', T), ('salt', S)):
                v = nc.createVariable(name, 'f8', ('depth', 'lat', 'lon'), fill_value=-1.e20)
                v.missing_value = -1.e20
                v[:] = np.where(np.isfinite(data), data, -1.e20)


# Sea-ice fields Isca sends with SIS2 (Isca name -> ocean name): fluxes over its sea-ice
# tile at the ice surface temperature it was sent, with their derivatives
ISCA_TO_SIS2 = [('ISCA_ITAUX', 'ICE_TAUX'), ('ISCA_ITAUY', 'ICE_TAUY'),
                ('ISCA_ISWDN', 'ICE_SWDN'), ('ISCA_ILWDN', 'ICE_LWDN'),
                ('ISCA_ISENS', 'ICE_SENS'), ('ISCA_IDHDT', 'ICE_DHDT'),
                ('ISCA_IEVAP', 'ICE_EVAP'), ('ISCA_IDEDT', 'ICE_DEDT'),
                ('ISCA_IRAIN', 'ICE_RAIN'), ('ISCA_ISNOW', 'ICE_SNOW')]
# ... and the ice state SIS2 sends (ocean name -> Isca name)
SIS2_TO_ISCA = [('ICE_FRAC', 'OCN_ICEFR'), ('ICE_TSRF', 'OCN_ICETS'), ('ICE_ALB', 'OCN_ICEAL')]

SIS2_DIAG_TABLE = """
"ice_month", 1, "months", 1, "days", "time"
"ice_model", "siconc", "siconc", "ice_month", "all", .true., "none", 2
"ice_model", "sithick", "sithick", "ice_month", "all", .true., "none", 2
"ice_model", "EXT", "EXT", "ice_month", "all", .true., "none", 2
"ice_model", "HI", "HI", "ice_month", "all", .true., "none", 2
"ice_model", "TS", "TS", "ice_month", "all", .true., "none", 2
"ice_model", "ALB", "ALB", "ice_month", "all", .true., "none", 2
"ice_model", "FRAZIL", "FRAZIL", "ice_month", "all", .true., "none", 2
"ice_model", "SH", "SH", "ice_month", "all", .true., "none", 2
"ice_model", "LW", "LW", "ice_month", "all", .true., "none", 2
"ice_model", "SW", "SW", "ice_month", "all", .true., "none", 2
"ice_model", "BHEAT", "BHEAT", "ice_month", "all", .true., "none", 2
"""

# SIS2 parameters (beyond the grid and topography, which are copied from MOM_input)
BASE_SIS_INPUT = {
    'CONSTANT_COSZEN_IC': 0.5,   # no diurnal cycle in the idealised setups
    'ADD_DIURNAL_SW': False,
    'ICE_STATS_INTERVAL': 30.0,
    'MERGED_CONTINUITY': True,   # needed when ice and ocean steps are interspersed
}
# the MOM_input entries SIS2 needs too
SIS_FROM_MOM = ('REENTRANT_X', 'REENTRANT_Y', 'TRIPOLAR_N', 'NIGLOBAL', 'NJGLOBAL', 'NIHALO', 'NJHALO',
                'LAYOUT', 'INPUTDIR', 'GRID_CONFIG', 'GRID_FILE', 'TOPO_CONFIG', 'TOPO_FILE',
                'TOPO_VARNAME', 'MAXIMUM_DEPTH', 'MINIMUM_DEPTH', 'MASKING_DEPTH')


class SIS2Mixin(object):
    """MOM6 with GFDL's SIS2 sea ice in the same executable (mom6_sis2_isca.x), stepped
    together by SIS2's combined ice-ocean driver. Isca must run with surface_tiles: its
    sea-ice tile sends fluxes over the ice (ISCA_TO_SIS2), and SIS2 sends back the ice
    fraction, surface temperature and albedo (SIS2_TO_ISCA)."""
    fields_to_isca = OCEAN_TO_ISCA + SIS2_TO_ISCA
    fields_from_isca = ISCA_TO_OCEAN + ISCA_TO_SIS2
    sis_input = {}

    def prepare_rundir(self, exp, cold_start):
        super(SIS2Mixin, self).prepare_rundir(exp, cold_start)
        odir = self.ocean_dir(exp)
        mom = {}
        with open(os.path.join(odir, 'MOM_input')) as f:
            for line in f:
                if '=' in line:
                    k, v = line.split('=', 1)
                    mom[k.strip()] = v.split('!')[0].strip()
        sis = {k: mom[k] for k in SIS_FROM_MOM if k in mom}
        sis.update(BASE_SIS_INPUT)
        sis['DT_ICE_DYNAMICS'] = float(exp.coupling_period)
        sis.update(self.sis_input)
        with open(os.path.join(odir, 'SIS_input'), 'w') as f:
            f.write(format_mom_input(sis))
        open(os.path.join(odir, 'SIS_override'), 'w').close()
        # the combined ice-ocean driver: ice slow thermodynamics, then ocean thermodynamics,
        # then dynamics of both, within each coupling step
        with open(os.path.join(odir, 'CIOD_input'), 'w') as f:
            # USE_INTERSPERSE_BUG (default True) never passes the ocean state, frazil included,
            # to the ice after the first step
            f.write(format_mom_input({'INTERSPERSE_ICE_OCEAN': True, 'SINGLE_MOM_CALL': False,
                                      'USE_INTERSPERSE_BUG': False}))
        with open(os.path.join(odir, 'diag_table'), 'a') as f:
            f.write(SIS2_DIAG_TABLE)
        with open(os.path.join(odir, 'input.nml'), 'a') as f:
            # the ice starts from nothing if there is no ice restart (e.g. a run without
            # SIS2 continued with it)
            ice_restart = os.path.exists(os.path.join(odir, 'INPUT', 'ice_model.res.nc'))
            f.write("&SIS_input_nml\n  output_directory = './'\n  input_filename = '%s'\n"
                    "  restart_input_dir = 'INPUT/'\n  restart_output_dir = 'RESTART/'\n"
                    "  parameter_filename = 'SIS_input', 'SIS_override'\n/\n"
                    % ('r' if ice_restart else 'n'))
            f.write("&ice_ocean_driver_nml\n  output_directory = './'\n"
                    "  parameter_filename = 'CIOD_input'\n/\n")
            f.write("&ice_model_nml\n/\n")

    def save_output(self, exp, outdir):
        super(SIS2Mixin, self).save_output(exp, outdir)
        odir = self.ocean_dir(exp)
        dest = os.path.join(outdir, 'ocean')
        for f in sorted(os.listdir(odir)):
            if (f.startswith('ice_') and '.nc' in f) or f in ('SIS_parameter_doc.all', 'SIS_input',
                                                               'seaice.stats', 'seaice.stats.nc'):
                shutil.copy(os.path.join(odir, f), dest)


class MOM6SIS2Ocean(SIS2Mixin, MOM6Ocean):
    """MOM6Ocean (idealised domains) with SIS2 sea ice."""

    def __init__(self, executable=None, **kwargs):
        executable = executable or os.path.expanduser(
            '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build_sis2/mom6_sis2_isca.x')
        super(MOM6SIS2Ocean, self).__init__(executable, **kwargs)


class MOM6SIS2EarthOcean(SIS2Mixin, MOM6EarthOcean):
    """MOM6EarthOcean with SIS2 sea ice."""


class MOM6SIS2TripolarOcean(SIS2Mixin, MOM6TripolarOcean):
    """MOM6TripolarOcean with SIS2 sea ice."""

    def __init__(self, executable=None, **kwargs):
        executable = executable or os.path.expanduser(
            '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build_sis2/mom6_sis2_isca.x')
        super(MOM6SIS2TripolarOcean, self).__init__(executable, **kwargs)

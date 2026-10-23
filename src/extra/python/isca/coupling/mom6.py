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

from isca.coupling import OceanComponent, OCEAN_RESTART

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
    'ADJUST_NET_FRESH_WATER_TO_ZERO': True,
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
    lat_c, _ = gaussian_grid(nlat)
    rows = ocean_rows(nlat, max_lat)
    lat = lat_c[rows]
    lon = np.arange(nlon) * 360. / nlon
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

    def __init__(self, executable, num_cores=4, max_lat=80., nk=30, depth=4000.,
                 layout=None, mom_input=None, sst_func=None):
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
        return nlon, len(ocean_rows(nlat, self.max_lat))

    def mpmd_options(self, exp):
        return '-wdir ocean'

    def command(self, exp, runtime):
        return [self.executable]

    def ocean_dir(self, exp):
        return os.path.join(exp.rundir, 'ocean')

    def prepare_rundir(self, exp, cold_start):
        nlon, nlat = exp.isca_grid_shape()
        odir = self.ocean_dir(exp)
        indir = os.path.join(odir, 'INPUT')
        os.makedirs(indir, exist_ok=True)
        os.makedirs(os.path.join(odir, 'RESTART'), exist_ok=True)

        nj = write_supergrid(os.path.join(indir, 'ocean_hgrid.nc'), nlon, nlat, self.max_lat)
        write_vgrid(os.path.join(indir, 'vgrid.nc'), self.nk, self.depth)
        write_initial_ts(os.path.join(indir, 'ts_init.nc'), nlon, nlat, self.max_lat, self.sst_func)

        params = dict(self.mom_input)
        params.update({'NIGLOBAL': nlon, 'NJGLOBAL': nj, 'NK': self.nk})
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

    def load_restart(self, exp, indir):
        src = os.path.join(indir, 'ocean')
        dest = os.path.join(self.ocean_dir(exp), 'INPUT')
        os.makedirs(dest, exist_ok=True)
        if os.path.isdir(src):
            for f in os.listdir(src):
                shutil.copy(os.path.join(src, f), dest)

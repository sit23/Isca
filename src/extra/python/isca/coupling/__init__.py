"""Running Isca coupled to an external ocean model through OASIS3-MCT.

Isca and the ocean are separate programs started by one MPMD mpirun. Isca's
side of the exchange is src/coupler/ocean_coupler.F90 (switched on with
mixed_layer_nml: ocean_coupling = 'oasis'); this package writes the OASIS
configuration (namcouple), prepares the ocean component's files, launches
the run and carries the OASIS coupling restart files from one run to the next
inside Isca's normal restart archives.

    from isca.coupling import CoupledExperiment, SlabOcean
    ocean = SlabOcean(depth=2.5)
    exp = CoupledExperiment('my_exp', codebase=cb, ocean=ocean, coupling_period=3600)
"""

import os
import shutil
import sys

import netCDF4
import numpy as np

from isca.experiment import Experiment

# Fields Isca sends (Isca name -> ocean name) and receives (ocean name -> Isca name).
# The Isca names must match src/coupler/ocean_coupler.F90.
ISCA_TO_OCEAN = [('ISCA_TAUX', 'OCN_TAUX'), ('ISCA_TAUY', 'OCN_TAUY'),
                 ('ISCA_SWNET', 'OCN_SWNET'), ('ISCA_LWNET', 'OCN_LWNET'),
                 ('ISCA_SENS', 'OCN_SENS'), ('ISCA_EVAP', 'OCN_EVAP'),
                 ('ISCA_PRECIP', 'OCN_PRECIP'), ('ISCA_PSURF', 'OCN_PSURF'),
                 ('ISCA_RUNOFF', 'OCN_RUNOFF')]
OCEAN_TO_ISCA = [('SST_OCN', 'OCN_SST'), ('FRAC_OCN', 'OCN_FRAC')]

ISCA_GRID = 'isca'
ISCA_GRID_RUNOFF = 'iscr'   # the Isca grid without its land mask (runoff only)
ISCA_RESTART = 'isca_rst.nc'    # OASIS restart file for the fluxes Isca sends
OCEAN_RESTART = 'ocn_rst.nc'    # OASIS restart file for the fields the ocean sends


class OceanComponent(object):
    """Interface for an ocean program that Isca can be coupled to."""
    grid_name = ISCA_GRID   # OASIS grid name used by the ocean
    fields_to_isca = OCEAN_TO_ISCA   # (ocean name, Isca name) of the fields the ocean sends
    fields_from_isca = ISCA_TO_OCEAN # (Isca name, ocean name) of the fields Isca sends
    num_cores = 1
    # Lag of the ocean-to-Isca fields [s]. None: one coupling period (concurrent
    # coupling, both models run at once). 0: sequential, Isca waits for the
    # ocean's step and uses its SST without lag (no ocean OASIS restart needed).
    sst_lag = None

    def ocean_grid_shape(self, nlon, nlat):
        """(nx, ny) of the ocean grid, given the Isca grid size."""
        return nlon, nlat

    def isca_ocean_mask(self, nlon, nlat):
        """Boolean (nlat, nlon) array: Isca cells the ocean model covers (before land is removed)."""
        return np.ones((nlat, nlon), dtype=bool)

    def command(self, exp, runtime):
        """Program and arguments to launch (a list of strings)."""
        raise NotImplementedError

    def mpmd_options(self, exp):
        """Extra mpirun options for the ocean's MPMD context (e.g. -wdir)."""
        return ''

    def prepare_rundir(self, exp, cold_start):
        """Put any files the ocean needs into exp.rundir before a run."""
        pass

    def save_restart(self, exp, resdir):
        """Copy the ocean model's own restart files into Isca's RESTART dir."""
        pass

    def load_restart(self, exp, indir):
        """Bring the ocean's restart files out of the extracted restart archive."""
        pass

    def save_output(self, exp, outdir):
        """Copy the ocean's diagnostic output into the run's output directory."""
        pass


class SlabOcean(OceanComponent):
    """The Python slab ocean in slab_ocean.py, on the Isca grid. For testing the coupling."""

    def __init__(self, depth=2.5, tconst=285., delta_t=40., max_lat=90.):
        self.depth = depth
        self.tconst = tconst
        self.delta_t = delta_t
        self.max_lat = max_lat

    def isca_ocean_mask(self, nlon, nlat):
        lat = gaussian_latitudes(nlat)
        return np.broadcast_to((np.abs(lat) <= self.max_lat)[:, None], (nlat, nlon)).copy()

    def command(self, exp, runtime):
        nlon, nlat = exp.isca_grid_shape()
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'slab_ocean.py')
        return [sys.executable, script, '--nlon', nlon, '--nlat', nlat,
                '--runtime', runtime, '--coupling-period', exp.coupling_period,
                '--depth', self.depth, '--tconst', self.tconst,
                '--delta-t', self.delta_t, '--max-lat', self.max_lat,
                '--restart-file', OCEAN_RESTART]


def gaussian_latitudes(nlat):
    """Isca's Gaussian latitudes [degrees], south to north."""
    return np.degrees(np.arcsin(np.polynomial.legendre.leggauss(nlat)[0]))


def runoff_destinations(coupled, zsurf=None):
    """Where the runoff of every Isca cell goes.

    coupled: boolean (nlat, nlon), cells whose surface is the coupled ocean.
    zsurf:   optional surface height (nlat, nlon) [m]; water follows the
             steepest descent (8 neighbours, periodic in longitude) until it
             reaches a coupled cell. Cells with no downhill path (pits, flat
             land, or no zsurf at all) send their water to the nearest coupled
             cell by great-circle distance.
    Returns int (nlat, nlon): 1-based global index i + (j-1)*nlon of the
    destination cell, 0 for coupled cells.
    """
    nlat, nlon = coupled.shape
    if not coupled.any():
        raise ValueError('no coupled ocean cells to route runoff to')
    lat = np.radians(gaussian_latitudes(nlat))
    lon = np.radians(np.arange(nlon) * 360. / nlon)
    LON, LAT = np.meshgrid(lon, lat)
    xyz = np.stack([np.cos(LAT) * np.cos(LON), np.cos(LAT) * np.sin(LON), np.sin(LAT)], -1)
    ocean_j, ocean_i = np.nonzero(coupled)
    ocean_xyz = xyz[ocean_j, ocean_i]

    def nearest_ocean(j, i):
        k = np.argmax(ocean_xyz @ xyz[j, i])
        return ocean_j[k], ocean_i[k]

    dest = np.zeros((nlat, nlon), dtype=np.int32)
    for j in range(nlat):
        for i in range(nlon):
            if coupled[j, i]:
                continue
            jj, ii = j, i
            for _ in range(nlat * nlon):
                if coupled[jj, ii]:
                    break
                step = None
                if zsurf is not None:
                    best = zsurf[jj, ii]
                    for dj in (-1, 0, 1):
                        for di in (-1, 0, 1):
                            nj, ni = jj + dj, (ii + di) % nlon
                            if (dj or di) and 0 <= nj < nlat and zsurf[nj, ni] < best:
                                best, step = zsurf[nj, ni], (nj, ni)
                if step is None:
                    jj, ii = nearest_ocean(jj, ii)
                else:
                    jj, ii = step
            dest[j, i] = ii + 1 + jj * nlon
    return dest


def write_runoff_routing(fname, dest):
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('lat', dest.shape[0])
        nc.createDimension('lon', dest.shape[1])
        v = nc.createVariable('dest', 'i4', ('lat', 'lon'))
        v.long_name = '1-based global index (i + (j-1)*nlon) of the cell receiving this cell runoff; 0 = none'
        v[:] = dest


def write_namcouple(fname, runtime, coupling_period, atmos_dt, isca_shape, ocean_shape,
                    ocean_grid, ocean_coupling_period=None, nlogprt=1, sst_lag=None,
                    ocean_to_isca=OCEAN_TO_ISCA, surface_tiles=False, isca_to_ocean=ISCA_TO_OCEAN):
    """Write the OASIS3-MCT namcouple file.

    Fluxes from Isca are averaged by OASIS over the coupling period (LOCTRANS
    AVERAGE) and lagged by one atmospheric step; ocean fields are lagged by one
    ocean coupling period. Both components therefore run concurrently, with
    the ocean seeing the fluxes of the previous coupling period and Isca the
    SST at the end of the previous ocean step (the usual concurrent coupling).
    If the ocean grid differs from Isca's, fields are remapped conservatively.

    With surface_tiles (Isca's fractional surface tiles) Isca sends the fluxes of its
    ocean tile, per unit area of ocean, from every cell with ocean in it. They are then
    all remapped like runoff, from the unmasked grid by overlap area (DESTAREA), so the
    ocean gets exactly what each cell's ocean tile gave up.
    """
    ocean_coupling_period = ocean_coupling_period or coupling_period
    sst_lag = ocean_coupling_period if sst_lag is None else sst_lag
    same_grid = (ocean_grid == ISCA_GRID)
    nx, ny = isca_shape
    ox, oy = ocean_shape
    lines = ['$NFIELDS', '  %d' % (len(isca_to_ocean) + len(ocean_to_isca)), '$END',
             '$RUNTIME', '  %d' % runtime, '$END',
             '$NLOGPRT', '  %d 0' % nlogprt, '$END',
             '$STRINGS']
    conserv = 'CONSERV LR SCALAR LATLON 1 FRACAREA FIRST'
    # ocean cells under Isca land (coastlines differ) take the nearest Isca ocean value
    conserv_to_ocean = 'CONSERV LR SCALAR LATLON 1 FRACNNEI FIRST'
    for src, tgt in isca_to_ocean:
        ntrans = 1 if same_grid else 2
        # Runoff arrives in coastal cells that the ocean may barely cover. Isca has already
        # divided it by the covered fraction, so it must be remapped from the unmasked grid
        # by overlap area (DESTAREA); FRACNNEI would hand a sliver's value to whole ocean
        # cells lying under Isca land and multiply the water delivered. With surface tiles
        # every flux goes this way.
        by_overlap = (src == 'ISCA_RUNOFF') or surface_tiles
        lines += ['%s %s 1 %d %d %s EXPORTED' % (src, tgt, coupling_period, ntrans, ISCA_RESTART),
                  '%d %d %d %d %s %s LAG=+%d' % (nx, ny, ox, oy,
                                                 ISCA_GRID_RUNOFF if by_overlap and not same_grid else ISCA_GRID,
                                                 ocean_grid, atmos_dt),
                  'P 0 P 0',
                  'LOCTRANS' if same_grid else 'LOCTRANS SCRIPR',
                  'AVERAGE']
        if not same_grid:
            lines.append(conserv.replace('FRACAREA', 'DESTAREA') if by_overlap else conserv_to_ocean)
    for src, tgt in ocean_to_isca:
        if same_grid:
            lines += ['%s %s 1 %d 0 %s EXPORTED' % (src, tgt, ocean_coupling_period, OCEAN_RESTART),
                      '%d %d %d %d %s %s LAG=+%d' % (ox, oy, nx, ny, ocean_grid, ISCA_GRID,
                                                     sst_lag),
                      'P 0 P 0']
        else:
            lines += ['%s %s 1 %d 1 %s EXPORTED' % (src, tgt, ocean_coupling_period, OCEAN_RESTART),
                      '%d %d %d %d %s %s LAG=+%d' % (ox, oy, nx, ny, ocean_grid, ISCA_GRID,
                                                     sst_lag),
                      'P 0 P 0', 'SCRIPR',
                      # the ocean fraction must be normalised by the whole Isca cell
                      # area to give the covered fraction, and frazil heat too, so that
                      # every Isca cell gets the frazil of all the ocean under it
                      # (DESTAREA conserves the total); the SST by the covered area only
                      # (FRACAREA), so values are not diluted.
                      conserv.replace('FRACAREA', 'DESTAREA')
                      if src in ('FRAC_OCN', 'FRAZIL_OCN') else conserv]
    lines.append('$END')
    with open(fname, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def write_zero_flux_restart(fname, isca_shape, isca_to_ocean=ISCA_TO_OCEAN):
    """Cold-start OASIS restart file for the Isca fluxes: zero fluxes for the
    ocean's first coupling period (surface pressure set to 1e5 Pa)."""
    nx, ny = isca_shape
    with netCDF4.Dataset(fname, 'w') as nc:
        nc.createDimension('nx', nx)
        nc.createDimension('ny', ny)
        for src, _ in isca_to_ocean:
            v = nc.createVariable(src, 'f8', ('ny', 'nx'))
            v[:] = 1.e5 if src == 'ISCA_PSURF' else 0.


class CoupledExperiment(Experiment):
    """An Isca Experiment whose surface is coupled to an external ocean through OASIS3-MCT."""

    def __init__(self, name, codebase, ocean, coupling_period=3600, land_mask=None, zsurf=None,
                 **kwargs):
        """land_mask, zsurf: (nlat, nlon) arrays on the Isca grid, matching the land Isca
        is configured with, used to route runoff to the ocean (None: aquaplanet)."""
        super(CoupledExperiment, self).__init__(name, codebase, **kwargs)
        self.ocean = ocean
        self.land_mask = land_mask
        self.zsurf = zsurf
        self.coupling_period = int(coupling_period)
        self.on('run:ready', self._prepare_coupling)
        self.on('run:complete', self._save_coupling_restarts)

    def isca_grid_shape(self):
        nml = self.namelist.get('spectral_dynamics_nml', {})
        return int(nml.get('lon_max', 128)), int(nml.get('lat_max', 64))   # Isca's defaults

    def runtime_seconds(self):
        m = self.namelist['main_nml']
        return (int(m.get('days', 0)) * 86400 + int(m.get('hours', 0)) * 3600
                + int(m.get('minutes', 0)) * 60 + int(m.get('seconds', 0)))

    def _check_namelist(self):
        self.namelist.setdefault('mixed_layer_nml', {})['ocean_coupling'] = 'oasis'
        runtime, dt = self.runtime_seconds(), int(self.namelist['main_nml']['dt_atmos'])
        if self.coupling_period % dt != 0:
            raise ValueError('coupling_period (%d s) must be a multiple of dt_atmos (%d s)'
                             % (self.coupling_period, dt))
        if runtime % self.coupling_period != 0:
            raise ValueError('run length (%d s) must be a multiple of coupling_period (%d s)'
                             % (runtime, self.coupling_period))

    def run(self, i, num_cores=8, mpirun_opts='', **kwargs):
        self._check_namelist()
        runtime = self.runtime_seconds()
        program = ' '.join(str(a) for a in self.ocean.command(self, runtime))
        mpmd = ': %s -np %d %s' % (self.ocean.mpmd_options(self), self.ocean.num_cores, program)
        return super(CoupledExperiment, self).run(i, num_cores=num_cores, mpirun_opts=mpirun_opts,
                                                  mpmd_suffix=mpmd, **kwargs)

    def _prepare_coupling(self, exp, i):
        indir = os.path.join(self.rundir, 'INPUT')
        nlon, nlat = self.isca_grid_shape()
        ocean_shape = self.ocean.ocean_grid_shape(nlon, nlat)
        write_namcouple(os.path.join(self.rundir, 'namcouple'), self.runtime_seconds(),
                        self.coupling_period, int(self.namelist['main_nml']['dt_atmos']),
                        (nlon, nlat), ocean_shape, self.ocean.grid_name, sst_lag=self.ocean.sst_lag,
                        ocean_to_isca=self.ocean.fields_to_isca,
                        isca_to_ocean=self.ocean.fields_from_isca,
                        surface_tiles=bool(self.namelist['idealized_moist_phys_nml'].get('surface_tiles', False)))

        coupled = self.ocean.isca_ocean_mask(nlon, nlat)
        if self.land_mask is not None:
            coupled = coupled & ~np.asarray(self.land_mask, dtype=bool)
        write_runoff_routing(os.path.join(indir, 'runoff_routing.nc'),
                             runoff_destinations(coupled, self.zsurf))

        cold_start = True
        for rst in (ISCA_RESTART, OCEAN_RESTART):
            src = os.path.join(indir, rst)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(self.rundir, rst))
                cold_start = False
        if cold_start:
            self.log.info('Coupled cold start: zero fluxes for the first ocean coupling period')
            write_zero_flux_restart(os.path.join(self.rundir, ISCA_RESTART), (nlon, nlat),
                                    self.ocean.fields_from_isca)
        else:
            self.ocean.load_restart(self, indir)
        self.ocean.prepare_rundir(self, cold_start)

    def _save_coupling_restarts(self, exp, i):
        # OASIS rewrites its restart files in the run directory at the end of
        # the run; put them in RESTART/ so they go into Isca's restart archive.
        resdir = os.path.join(self.rundir, 'RESTART')
        for rst in (ISCA_RESTART, OCEAN_RESTART):
            src = os.path.join(self.rundir, rst)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(resdir, rst))
        self.ocean.save_restart(self, resdir)
        # Experiment.run only archives Isca's own diagnostics, and empties the
        # run directory afterwards, so the ocean's output is copied out here.
        outdir = os.path.join(self.datadir, self.runfmt % i)
        os.makedirs(outdir, exist_ok=True)
        self.ocean.save_output(self, outdir)


class SamudraOcean(OceanComponent):
    """The Samudra ML ocean emulator (samudra_ocean.py) on its 1-degree OM4 grid.

    python: interpreter of an environment with Samudra, torch, mpi4py and pyOASIS.
    config: Samudra eval-style YAML (data + model), ckpt: checkpoint.
    Couples every 5 days, sequentially (Isca uses the SST Samudra predicts for
    each 5-day interval without lag). Earth geometry only.
    extra_args: further samudra_ocean.py options, e.g. ['--ra-nu', 0.1].
    """
    grid_name = 'smdr'
    sst_lag = 0
    state_file = 'samudra_state.pt'

    def __init__(self, python, config, ckpt, threads=8, grid_shape=(360, 180), extra_args=()):
        self.python, self.config, self.ckpt = python, config, ckpt
        self.extra_args = list(extra_args)
        self.threads = threads
        self.grid_shape = grid_shape

    def ocean_grid_shape(self, nlon, nlat):
        return self.grid_shape

    def command(self, exp, runtime):
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'samudra_ocean.py')
        return ['env', 'OMP_NUM_THREADS=%d' % self.threads, self.python, script,
                '--config', self.config, '--ckpt', self.ckpt, '--runtime', runtime,
                '--coupling-period', exp.coupling_period, '--state-file', self.state_file,
                '--threads', self.threads] + self.extra_args

    def save_restart(self, exp, resdir):
        src = os.path.join(exp.rundir, self.state_file)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(resdir, self.state_file))

    def load_restart(self, exp, indir):
        src = os.path.join(indir, self.state_file)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(exp.rundir, self.state_file))

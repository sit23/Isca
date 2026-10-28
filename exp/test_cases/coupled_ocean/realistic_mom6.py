"""Isca with realistic continents (T42, RRTM) coupled to MOM6 with OM4-derived Earth geometry.

Based on exp/test_cases/realistic_continents, but the q-flux slab ocean is
replaced by MOM6 on a 2-degree grid with GFDL OM4 bathymetry (as regridded by
M2LiNES), initialised from an OM4 state, coupled hourly through OASIS3-MCT. Points under
the prescribed AMIP sea-ice climatology keep Isca's slab with ice albedo, and
the ocean beneath them receives no fluxes.

THERMO_ICE=1 replaces the prescribed ice with Isca's thermodynamic sea ice (zero-layer,
PR #266) on top of the coupled ocean, grown by MOM6's frazil. SURFACE_TILES=1 turns on
Isca's fractional surface tiles (ocean, land and slab ocean); it needs THERMO_ICE=1.
SIS2=1 instead runs GFDL's SIS2 sea ice with MOM6 (implies SURFACE_TILES=1, no THERMO_ICE).

Usage:
    NMONTHS=12 python realistic_mom6.py
"""
import os

import f90nml
import netCDF4
import numpy as np

from isca import IscaCodeBase, DiagTable, GFDL_BASE
from isca.coupling import CoupledExperiment
from isca.coupling.mom6 import MOM6EarthOcean, MOM6SIS2EarthOcean

NCORES = int(os.environ.get('NCORES', 16))
NMONTHS = int(os.environ.get('NMONTHS', 1))
DAYS = int(os.environ.get('DAYS', 30))
EXP_NAME = os.environ.get('EXP_NAME', 'realistic_mom6')
MOM6_EXE = os.environ.get('MOM6_ISCA_EXE', os.path.expanduser(
    '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build/mom6_isca.x'))
OM4 = os.environ.get('OM4_ZARR', 'https://nyu1.osn.mghpcc.org/m2lines-pubs/Samudra/v2026-09/om4_twodeg/OM4.zarr')

rc_dir = os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents')
land_file = os.path.join(GFDL_BASE, 'input/land_masks/era_land_t42.nc')
with netCDF4.Dataset(land_file) as nc:
    land = np.asarray(nc['land_mask'][:]) > 0.
    zsurf = np.asarray(nc['zsurf'][:])

THERMO_ICE = bool(int(os.environ.get('THERMO_ICE', 0)))
SIS2 = bool(int(os.environ.get('SIS2', 0)))
SURFACE_TILES = bool(int(os.environ.get('SURFACE_TILES', 0))) or SIS2
if SIS2 and THERMO_ICE:
    raise ValueError('SIS2=1 replaces Isca\'s sea ice: use THERMO_ICE=0')
if SURFACE_TILES and not (THERMO_ICE or SIS2):
    raise ValueError('SURFACE_TILES=1 needs THERMO_ICE=1 or SIS2=1 (tiles do not support the prescribed-ice lid)')

cb = IscaCodeBase.from_directory(GFDL_BASE)
if SIS2:
    ocean = MOM6SIS2EarthOcean(os.environ.get('MOM6_SIS2_EXE', os.path.expanduser(
        '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build_sis2/mom6_sis2_isca.x')),
        OM4, num_cores=int(os.environ.get('NCORES_OCEAN', 16)))
else:
    ocean = MOM6EarthOcean(MOM6_EXE, OM4, num_cores=int(os.environ.get('NCORES_OCEAN', 16)))
exp = CoupledExperiment(EXP_NAME, codebase=cb, ocean=ocean, coupling_period=3600,
                        land_mask=land, zsurf=zsurf)
exp.inputfiles = [land_file, os.path.join(GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc'),
                  os.path.join(rc_dir, 'input/siconc_clim_amip.nc')]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
for mod, field in [('dynamics', 'ps'), ('dynamics', 'bk'), ('dynamics', 'pk'), ('dynamics', 'zsurf'),
                   ('atmosphere', 'precipitation'), ('mixed_layer', 't_surf'),
                   ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
                   ('mixed_layer', 'coupled_ocean_frac'), ('mixed_layer', 'cpl_net_heat'),
                   ('atmosphere', 'cpl_runoff'),
                   ('dynamics', 'ucomp'), ('dynamics', 'vcomp'), ('dynamics', 'temp'),
                   ('dynamics', 'sphum'), ('rrtm_radiation', 'olr'), ('rrtm_radiation', 'toa_sw'),
                   ('mixed_layer', 'cpl_evap'), ('mixed_layer', 'cpl_precip'),
                   ('atmosphere', 'bucket_depth')]:
    diag.add_field(mod, field, time_avg=True)
ice_fields = (('h_thermo_ice', 't_ml', 'albedo', 'cpl_frazil', 'ice_ocean_heat') if THERMO_ICE else
              () if SIS2 else ('ice_conc',))
for field in ice_fields:
    diag.add_field('mixed_layer', field, time_avg=True)
exp.diag_table = diag

exp.namelist = f90nml.read(os.path.join(rc_dir, 'namelist_basefile.nml'))
exp.namelist['main_nml']['days'] = DAYS
# The ocean model replaces the q-flux; ice points keep the slab with ice albedo.
exp.namelist['mixed_layer_nml'].update({'do_qflux': False, 'load_qflux': False})
# Land bucket hydrology: its overflow is the river runoff routed to the river mouths.
# Without it there are no rivers, and land P-E is spread over the whole ocean.
exp.namelist['idealized_moist_phys_nml'].update({
    'bucket': True, 'init_bucket_depth_land': 0.15, 'max_bucket_depth_land': 0.15})
if THERMO_ICE:
    # ice forms and its base sits at MOM6's surface freezing point (linear, -0.054 K/psu at
    # about 34.7 psu); thermodynamic ice requires the explicit surface step
    t_freeze = 273.15 - 0.054 * 34.7
    exp.namelist['mixed_layer_nml'].update({
        'update_albedo_from_ice': False, 'update_land_mask_from_ice': False,
        'do_thermo_ice': True, 'do_explicit': True,
        'thermo_ice_albedo': float(os.environ.get('THERMO_ICE_ALBEDO', 0.55)),
        't_thermo_ice_base': t_freeze, 't_surf_freeze': t_freeze,
        'read_const_correct': False, 'read_nudge_out': False})
if SIS2:
    exp.namelist['mixed_layer_nml'].update({'update_albedo_from_ice': False,
                                            'update_land_mask_from_ice': False})
if SURFACE_TILES:
    exp.namelist['idealized_moist_phys_nml']['surface_tiles'] = True
# POLAR_ALBEDO=1: surface albedo rising towards the poles, albedo_value + (HIGHER_ALBEDO -
# albedo_value)*(|lat|/90)**ALBEDO_EXP, standing in for the polar clouds the clear-sky RRTM
# lacks. Applies to land too (albedo_choice other than 1 drops land_albedo_prefactor).
if int(os.environ.get('POLAR_ALBEDO', 0)):
    exp.namelist['mixed_layer_nml'].update({
        'albedo_choice': 4, 'higher_albedo': float(os.environ.get('HIGHER_ALBEDO', 0.55)),
        'albedo_exp': float(os.environ.get('ALBEDO_EXP', 4.))})

if __name__ == '__main__':
    if not os.environ.get('SKIP_COMPILE'):
        cb.compile()
    first = int(os.environ.get('FIRST_MONTH', 1))
    for i in range(first, NMONTHS + 1):
        exp.run(i, use_restart=(i > 1), num_cores=NCORES, overwrite_data=True,
                mpirun_opts='--bind-to none')

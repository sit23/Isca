"""Isca (Frierson grey-radiation aquaplanet, T42) coupled to a MOM6 aquaplanet ocean.

CONFIG=ridge adds a pole-to-pole ridge (two grid columns at 0E) in both models,
with a land bucket whose overflow is routed to the ocean (Ridge world, cf.
Enderton & Marshall 2009). RIDGE_WIDTH widens it to that many degrees of longitude
(rounded to whole grid columns), leaving a narrower basin as in Vallis & Farneti
(2009). MOM6_KD sets the ocean's background diapycnal diffusivity (m2/s).

THERMO_ICE=1 turns on Isca's thermodynamic sea ice (zero-layer, Lewis's PR #266):
over the coupled ocean it grows from MOM6's frazil heat and insulates the ocean,
in the polar caps it sits on Isca's slab. THERMO_ICE_ALBEDO sets its albedo (0.55).

SURFACE_TILES=1 turns on Isca's fractional surface tiles: each cell's ocean tile covers
what the ocean model covers, and all fluxes reach the ocean by overlap area.

The ocean is 4 km deep everywhere with no continents, on a grid whose cells
are exactly Isca's Gaussian-grid cells equatorward of +-80 degrees (MOM6
cannot include the poles); poleward of that Isca's own slab ocean is used.
The components exchange fields through OASIS3-MCT every hour.

Needs: an OASIS-enabled Isca build (GFDL_ENV=oasis_conda) and mom6_isca.x
(isca_ocean_coupling/mom6_isca/build.sh). Usage:
    NMONTHS=12 python frierson_mom6_aquaplanet.py
"""
import os
import sys

from isca import IscaCodeBase, DiagTable, GFDL_BASE
from isca.coupling import CoupledExperiment
from isca.coupling.mom6 import MOM6Ocean, write_isca_land_file
import numpy as np

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp', 'test_cases', 'frierson'))
from frierson_test_case import exp as frierson_exp

NCORES = int(os.environ.get('NCORES', 8))
NCORES_OCEAN = int(os.environ.get('NCORES_OCEAN', 8))
NMONTHS = int(os.environ.get('NMONTHS', 1))
DAYS = int(os.environ.get('DAYS', 30))
DT_ATMOS = 720
COUPLING_PERIOD = 3600
MOM6_EXE = os.environ.get('MOM6_ISCA_EXE', os.path.expanduser(
    '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build/mom6_isca.x'))
EXP_NAME = os.environ.get('EXP_NAME', 'coupled_ocean_mom6_aquaplanet')
CONFIG = os.environ.get('CONFIG', 'aqua')
COUPLING_PERIOD = int(os.environ.get('COUPLING_PERIOD', COUPLING_PERIOD))
NLON, NLAT = 128, 64

land = None
if CONFIG == 'ridge':
    land = np.zeros((NLAT, NLON), dtype=bool)
    # cells centred on 0E and 357.2E by default
    ncol = int(round(float(os.environ.get('RIDGE_WIDTH', 2 * 360. / NLON)) * NLON / 360.))
    land[:, np.arange(ncol) - ncol // 2] = True

cb = IscaCodeBase.from_directory(GFDL_BASE)

MAX_LAT = float(os.environ.get('MAX_LAT', 80.))
OCEAN_DLAT = float(os.environ['OCEAN_DLAT']) if 'OCEAN_DLAT' in os.environ else None
ocean = MOM6Ocean(MOM6_EXE, num_cores=NCORES_OCEAN, max_lat=MAX_LAT, land_mask=land,
                  regular_dlat=OCEAN_DLAT)
if 'MOM6_KD' in os.environ:
    ocean.mom_input['KD'] = float(os.environ['MOM6_KD'])
exp = CoupledExperiment(EXP_NAME, codebase=cb, ocean=ocean, coupling_period=COUPLING_PERIOD,
                        land_mask=land)

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
for mod, field in [('dynamics', 'ps'), ('dynamics', 'bk'), ('dynamics', 'pk'),
                   ('atmosphere', 'precipitation'), ('mixed_layer', 't_surf'),
                   ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
                   ('mixed_layer', 'coupled_ocean_frac'), ('two_stream', 'olr'),
                   ('two_stream', 'swdn_toa'), ('two_stream', 'swdn_sfc'), ('two_stream', 'lwdn_sfc'),
                   ('two_stream', 'lwup_sfc'), ('dynamics', 'sphum'), ('dynamics', 'ucomp'),
                   ('dynamics', 'vcomp'), ('dynamics', 'temp'), ('dynamics', 'height'),
                   ('mixed_layer', 'cpl_net_heat'), ('mixed_layer', 'cpl_evap'),
                   ('mixed_layer', 'cpl_precip'), ('atmosphere', 'cpl_runoff')]:
    diag.add_field(mod, field, time_avg=True)
THERMO_ICE = bool(int(os.environ.get('THERMO_ICE', 0)))
SURFACE_TILES = bool(int(os.environ.get('SURFACE_TILES', 0)))
if THERMO_ICE:
    for field in ('h_thermo_ice', 't_ml', 'albedo', 'flux_thermo_ice', 'cpl_frazil', 'ice_ocean_heat'):
        diag.add_field('mixed_layer', field, time_avg=True)
exp.diag_table = diag

exp.namelist = frierson_exp.namelist.copy()
exp.namelist['main_nml']['days'] = DAYS
exp.namelist['main_nml']['dt_atmos'] = DT_ATMOS
if SURFACE_TILES:
    exp.namelist['idealized_moist_phys_nml']['surface_tiles'] = True
# slab depth used only for the polar caps outside the MOM6 domain
exp.namelist['mixed_layer_nml']['depth'] = 20.
if THERMO_ICE:
    # ice forms and its base sits at MOM6's surface freezing point (linear, -0.054 K/psu at
    # about 34.7 psu); thermodynamic ice requires the explicit surface step
    t_freeze = 273.15 - 0.054 * 34.7
    exp.namelist['mixed_layer_nml'].update({
        'do_thermo_ice': True, 'do_explicit': True,
        'thermo_ice_albedo': float(os.environ.get('THERMO_ICE_ALBEDO', 0.55)),
        't_thermo_ice_base': t_freeze, 't_surf_freeze': t_freeze,
        'read_const_correct': False, 'read_nudge_out': False})

if land is not None:
    land_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'land_%s.nc' % CONFIG)
    write_isca_land_file(land_file, land)
    exp.inputfiles = [land_file]
    exp.namelist['idealized_moist_phys_nml'].update({
        'land_option': 'input', 'land_file_name': 'INPUT/land_%s.nc' % CONFIG,
        'bucket': True, 'init_bucket_depth_land': 0.15, 'max_bucket_depth_land': 0.15})
    exp.namelist['mixed_layer_nml'].update({
        'land_option': 'input', 'land_h_capacity_prefactor': 0.05})   # ~1 m of water on land

if __name__ == '__main__':
    if not os.environ.get('SKIP_COMPILE'):
        cb.compile()
    # several coupled jobs share this machine: let the OS place processes rather than
    # every mpirun binding its ranks to the same cores
    opts = os.environ.get('MPIRUN_OPTS', '--bind-to none')
    first = int(os.environ.get('FIRST_MONTH', 1))
    for i in range(first, NMONTHS + 1):
        exp.run(i, use_restart=(i > 1), num_cores=NCORES, overwrite_data=True, mpirun_opts=opts)

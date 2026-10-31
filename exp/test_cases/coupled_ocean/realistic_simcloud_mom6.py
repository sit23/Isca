"""Earth-like Isca with Socrates radiation and the SimCloud diagnostic clouds (Liu et al. 2021),
coupled to MOM6 and SIS2 on the OM4 2-degree grid: the configuration for tuning Isca towards
present-day Earth.

The atmosphere is the published SimCloud setup (exp/test_cases/simple_clouds/socrates_simcloud.py:
T42, 25 levels, simple Betts-Miller convection, linear cloud scheme with marine stratocumulus),
but with the prescribed SST and ice replaced by MOM6 + SIS2 as in realistic_mom6.py (SIS2=1):
fractional surface tiles, SIS2 sea ice with snow, land-bucket rivers routed to the mouths of
their real basins (DDM30 network).
There is a full seasonal cycle; the forcing is the same every year (present day):
  CO2 = CO2_PPMV (default 380, the 1991-2020 mean), CH4 and N2O on at Socrates' defaults
  (about 1.8 ppmv and 325 ppbv), 1990 ozone climatology, solar constant 1361 W/m2.

Tuning parameters can be set from the environment without editing this file, e.g.
    TUNE='large_scale_cloud_nml:linear_a_surf=40,mixed_layer_nml:albedo_value=0.1'
(namelist:name=value pairs, comma separated; values are parsed as Python literals).

Usage:
    EXP_NAME=SC1 NMONTHS=36 NCORES=8 NCORES_OCEAN=8 python realistic_simcloud_mom6.py
"""
import ast
import os

import netCDF4
import numpy as np

from isca import SocratesCodeBase, DiagTable, Namelist, GFDL_BASE
from isca.coupling import CoupledExperiment
from isca.coupling.mom6 import MOM6SIS2EarthOcean

NCORES = int(os.environ.get('NCORES', 8))
NCORES_OCEAN = int(os.environ.get('NCORES_OCEAN', 8))
NMONTHS = int(os.environ.get('NMONTHS', 1))
EXP_NAME = os.environ.get('EXP_NAME', 'realistic_simcloud_mom6')
CO2_PPMV = float(os.environ.get('CO2_PPMV', 380.))
MOM6_SIS2_EXE = os.environ.get('MOM6_SIS2_EXE', os.path.expanduser(
    '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build_sis2/mom6_sis2_isca.x'))
OM4 = os.environ.get('OM4_ZARR', 'https://nyu1.osn.mghpcc.org/m2lines-pubs/Samudra/v2026-09/om4_twodeg/OM4.zarr')
# rivers: 'ddm30' routes land runoff to the mouths of the real basins (isca.coupling.rivers);
# RIVERS=downhill uses the idealised default (steepest descent on the T42 topography)
RIVERS = None if os.environ.get('RIVERS', 'ddm30') == 'downhill' else os.environ.get('RIVERS', 'ddm30')

sc_dir = os.path.join(GFDL_BASE, 'exp/test_cases/simple_clouds')
land_file = os.path.join(sc_dir, 'input/era_land_t42_filtered.nc')
with netCDF4.Dataset(land_file) as nc:
    land = np.asarray(nc['land_mask'][:]) > 0.
    zsurf = np.asarray(nc['zsurf'][:])

cb = SocratesCodeBase.from_directory(GFDL_BASE)
ocean = MOM6SIS2EarthOcean(MOM6_SIS2_EXE, OM4, num_cores=NCORES_OCEAN)
exp = CoupledExperiment(EXP_NAME, codebase=cb, ocean=ocean, coupling_period=3600,
                        land_mask=land, zsurf=zsurf, river_network=RIVERS)
exp.inputfiles = [land_file, os.path.join(GFDL_BASE, 'input/rrtm_input_files/ozone_1990.nc')]

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
fields = [('dynamics', f) for f in ('ps', 'bk', 'pk', 'zsurf', 'sphum', 'ucomp', 'vcomp', 'omega', 'temp')]
fields += [('atmosphere', 'precipitation'), ('atmosphere', 'cpl_runoff'), ('atmosphere', 'bucket_depth')]
fields += [('mixed_layer', f) for f in ('t_surf', 'flux_t', 'flux_lhe', 'coupled_ocean_frac',
                                         'cpl_net_heat', 'cpl_evap', 'cpl_precip', 'albedo')]
fields += [('socrates', f) for f in (
    'soc_olr', 'soc_toa_sw', 'soc_toa_sw_up', 'soc_toa_sw_down', 'soc_olr_clr', 'soc_toa_sw_clr',
    'soc_toa_sw_up_clr', 'soc_surf_flux_lw', 'soc_surf_flux_sw', 'soc_surf_flux_lw_down',
    'soc_surf_flux_sw_down', 'soc_surf_flux_lw_clr', 'soc_surf_flux_sw_clr',
    'soc_surf_flux_lw_down_clr', 'soc_surf_flux_sw_down_clr')]
fields += [('cloud_simple', f) for f in ('cf', 'reff_rad', 'frac_liq', 'qcl_rad')]
fields += [('cloud_cover', f) for f in ('tot_cld_amt', 'high_cld_amt', 'mid_cld_amt', 'low_cld_amt')]
for mod, field in fields:
    diag.add_field(mod, field, time_avg=True)
exp.diag_table = diag

exp.namelist = Namelist({
    'main_nml': {'days': 30, 'hours': 0, 'minutes': 0, 'seconds': 0, 'dt_atmos': 720,
                 'current_date': [1, 1, 1, 0, 0, 0], 'calendar': 'thirty_day'},
    'socrates_rad_nml': {
        'stellar_constant': 1361.,
        'lw_spectral_filename': os.path.join(GFDL_BASE, 'src/atmos_param/socrates/src/trunk/data/spectra/ga3_1/sp_lw_ga3_1'),
        'sw_spectral_filename': os.path.join(GFDL_BASE, 'src/atmos_param/socrates/src/trunk/data/spectra/ga3_1/sp_sw_ga3_0'),
        'do_read_ozone': True, 'ozone_file_name': 'ozone_1990', 'ozone_field_name': 'ozone_1990',
        'co2_ppmv': CO2_PPMV, 'inc_ch4': True, 'inc_n2o': True,
        'dt_rad': 4320, 'store_intermediate_rad': True, 'chunk_size': 16,
        'use_pressure_interp_for_half_levels': False, 'tidally_locked': False},
    'idealized_moist_phys_nml': {
        'do_damping': True, 'turb': True, 'mixed_layer_bc': True, 'do_virtual': False,
        'do_simple': True, 'two_stream_gray': False, 'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', 'do_cloud_simple': True,
        'land_option': 'input', 'land_file_name': 'INPUT/era_land_t42_filtered.nc',
        'land_roughness_prefactor': 10.0,
        'roughness_mom': 2.e-04, 'roughness_heat': 2.e-04, 'roughness_moist': 2.e-04,
        # land bucket: its overflow is the river runoff routed to the river mouths
        'bucket': True, 'init_bucket_depth_land': 0.15, 'max_bucket_depth_land': 0.15,
        'surface_tiles': True},
    'cloud_simple_nml': {'do_qcl_with_temp': True, 'do_cloud_cover_diags': True,
                         'do_add_stratocumulus': True, 'reff_liq': 14, 'reff_ice': 25, 'qcl_val': 0.18},
    'large_scale_cloud_nml': {'cf_diag_formula_name': 'linear', 'do_adjust_cld_by_omega': False,
                              'do_freezedry': True, 'qv_polar_val': 0.006, 'freezedry_power': 2.5,
                              'do_fitted_rhcrit': False, 'linear_a_surf': 42, 'linear_a_top': 13,
                              'linear_power': 11},
    'marine_strat_cloud_nml': {'sc_diag_method': 'Park_ELF', 'intermediate_outputs_diags': False,
                               'dthdp_min_threshold': -0.08, 'park_a': 1.3, 'park_b': -0.1},
    'cloud_cover_diag_nml': {'overlap_assumption': 'maximum-random', 'mid_cld_bottom': 7.0e4,
                             'high_cld_bottom': 4.0e4, 'cf_min': 0},
    'vert_turb_driver_nml': {'do_mellor_yamada': False, 'do_diffusivity': True, 'do_simple': True,
                             'constant_gust': 0.0, 'use_tau': False},
    'diffusivity_nml': {'do_entrain': True, 'do_simple': True},
    'surface_flux_nml': {'use_virtual_temp': False, 'do_simple': True, 'old_dtaudv': True,
                         'land_humidity_prefactor': 1, 'land_evap_prefactor': 0.6},
    'atmosphere_nml': {'idealized_moist_model': True},
    # the ocean model replaces the slab and q-flux; SIS2 supplies the ice (no Isca ice albedo
    # or land mask from ice); the slab only remains on land and outside MOM6's domain
    'mixed_layer_nml': {
        'tconst': 285., 'prescribe_initial_dist': True, 'evaporation': True, 'depth': 20.0,
        'land_option': 'input', 'land_h_capacity_prefactor': 0.1,
        'albedo_value': 0.12, 'land_albedo_prefactor': 1.3,
        'do_qflux': False, 'load_qflux': False,
        'update_albedo_from_ice': False, 'update_land_mask_from_ice': False},
    'qe_moist_convection_nml': {'rhbm': 0.7, 'Tmin': 160., 'Tmax': 350.},
    'lscale_cond_nml': {'do_simple': True, 'do_evap': True},
    # tcmin/tcmax: with do_simple=False (saturation over ice below freezing) the default table
    # stops at -160 C, and the 0.1 hPa level reaches 113 K in the southern polar night
    'sat_vapor_pres_nml': {'do_simple': True, 'construct_table_wrt_liq_and_ice': True,
                           'show_all_bad_values': True, 'tcmin': -200, 'tcmax': 150},
    'damping_driver_nml': {'do_rayleigh': True, 'trayfric': -0.5, 'sponge_pbottom': 150.,
                           'do_conserve_energy': True},
    'diag_manager_nml': {'mix_snapshot_average_fields': False},
    'fms_nml': {'domains_stack_size': 600000},
    'fms_io_nml': {'threading_write': 'single', 'fileset_write': 'single'},
    'spectral_dynamics_nml': {
        'damping_order': 4, 'water_correction_limit': 200.e2, 'reference_sea_level_press': 1.0e5,
        'num_levels': 25, 'valid_range_t': [100., 800.], 'initial_sphum': [2.e-6],
        'vert_coord_option': 'uneven_sigma', 'surf_res': 0.03, 'scale_heights': 11.0,
        'exponent': 7.0, 'robert_coeff': 0.03, 'ocean_topog_smoothing': 0.8},
    'spectral_init_cond_nml': {'topog_file_name': 'era_land_t42_filtered.nc', 'topography_option': 'input'},
})

for item in filter(None, os.environ.get('TUNE', '').split(',')):
    key, value = item.split('=', 1)
    nml, name = key.split(':')
    exp.namelist.setdefault(nml.strip(), {})[name.strip()] = ast.literal_eval(value.strip())

if __name__ == '__main__':
    if not os.environ.get('SKIP_COMPILE'):
        cb.compile()
    first = int(os.environ.get('FIRST_MONTH', 1))
    for i in range(first, NMONTHS + 1):
        exp.run(i, use_restart=(i > 1), num_cores=NCORES, overwrite_data=True,
                mpirun_opts=os.environ.get('MPIRUN_OPTS', '--bind-to none'))

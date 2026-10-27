"""Uncoupled check of the thermodynamic sea ice (Neil Lewis, PR #266) after merging it
into the coupled-ocean branch: his grey-radiation slab test case, run for NMONTHS months
with the IscaCodeBase build the coupled runs use.

Usage: NMONTHS=2 python thermo_ice_uncoupled.py
"""
import os

from isca import IscaCodeBase, DiagTable, Experiment, GFDL_BASE

src = open(os.path.join(GFDL_BASE, 'exp/test_cases/thermodynamic_sea_ice/grey_thermo_ice_test_case.py')).read()
# reuse only the namelist from Neil's script (it compiles and runs on import otherwise)
ns = {}
exec(src[src.index('exp.namelist = namelist = Namelist('):src.index('# Lets do a run!')]
     .replace('exp.namelist = namelist = ', 'namelist = '), {'Namelist': __import__('isca').Namelist}, ns)
nml = ns['namelist']
nml['mixed_layer_nml'].pop('do_var_thermo_ice_albedo', None)   # not an option of the merged code
nml['main_nml']['days'] = 30
nml['qe_moist_convection_nml'].pop('precision', None)   # removed from Isca since the PR was written

cb = IscaCodeBase.from_directory(GFDL_BASE)
exp = Experiment(os.environ.get('EXP_NAME', 'thermo_ice_uncoupled'), codebase=cb)
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
for mod, f in [('mixed_layer', 't_surf'), ('mixed_layer', 'h_thermo_ice'), ('mixed_layer', 't_ml'),
               ('mixed_layer', 'albedo'), ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
               ('two_stream', 'olr'), ('two_stream', 'swdn_toa'), ('atmosphere', 'precipitation')]:
    diag.add_field(mod, f, time_avg=True)
exp.diag_table = diag
exp.namelist = nml

if __name__ == '__main__':
    for i in range(1, int(os.environ.get('NMONTHS', 2)) + 1):
        exp.run(i, use_restart=(i > 1), num_cores=8, overwrite_data=True, mpirun_opts='--bind-to none')

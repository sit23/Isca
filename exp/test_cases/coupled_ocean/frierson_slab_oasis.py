"""Phase-1 test of Isca's OASIS3-MCT ocean coupling.

Runs the Frierson aquaplanet twice from the same initial state:
  * 'reference': Isca's own slab ocean (mixed_layer.F90), no coupling;
  * 'coupled'  : the same slab physics in a separate program
                 (isca/coupling/slab_ocean.py) coupled through OASIS3-MCT
                 every atmospheric timestep.
The two climates should agree to within internal variability. Differences
can only come from the explicit (coupled) vs implicit (built-in) surface
temperature step and the one-step SST lag, both small for these settings.

Needs an Isca build with OASIS (GFDL_ENV=oasis_conda, OASIS_ROOT set) and
pyOASIS on the PYTHONPATH.
"""
import os
import sys

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.coupling import CoupledExperiment, SlabOcean

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp', 'test_cases', 'frierson'))

NCORES = 8
DT_ATMOS = 720
NMONTHS = int(os.environ.get('NMONTHS', 1))
DEPTH = 2.5

cb = IscaCodeBase.from_directory(GFDL_BASE)


def make_diag():
    diag = DiagTable()
    diag.add_file('atmos_monthly', 30, 'days', time_units='days')
    diag.add_field('dynamics', 'ps', time_avg=True)
    diag.add_field('dynamics', 'bk')
    diag.add_field('dynamics', 'pk')
    diag.add_field('atmosphere', 'precipitation', time_avg=True)
    diag.add_field('mixed_layer', 't_surf', time_avg=True)
    diag.add_field('mixed_layer', 'flux_t', time_avg=True)
    diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
    diag.add_field('dynamics', 'sphum', time_avg=True)
    diag.add_field('dynamics', 'ucomp', time_avg=True)
    diag.add_field('dynamics', 'vcomp', time_avg=True)
    diag.add_field('dynamics', 'temp', time_avg=True)
    return diag


def make_namelist(days):
    # The Frierson test case namelist (exp/test_cases/frierson), with the run length set here.
    from frierson_test_case import exp as frierson_exp
    nml = Namelist(frierson_exp.namelist.copy())
    nml['main_nml']['days'] = days
    nml['main_nml']['dt_atmos'] = DT_ATMOS
    nml['mixed_layer_nml']['depth'] = DEPTH
    return nml


def run(exp, nmonths):
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2, nmonths + 1):
        exp.run(i, num_cores=NCORES, overwrite_data=True)


if __name__ == '__main__':
    cb.compile()

    which = sys.argv[1:] or ['reference', 'coupled']

    if 'reference' in which:
        ref = Experiment('coupled_ocean_phase1_reference', codebase=cb)
        ref.diag_table = make_diag()
        ref.namelist = make_namelist(30)
        run(ref, NMONTHS)

    if 'coupled' in which:
        ocean = SlabOcean(depth=DEPTH, tconst=285., delta_t=40.)
        cpl = CoupledExperiment('coupled_ocean_phase1_slab', codebase=cb, ocean=ocean,
                                coupling_period=DT_ATMOS)
        cpl.diag_table = make_diag()
        cpl.namelist = make_namelist(30)
        run(cpl, NMONTHS)

"""Isca (Frierson grey-radiation aquaplanet, T42) coupled to a MOM6 aquaplanet ocean.

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
from isca.coupling.mom6 import MOM6Ocean

sys.path.insert(0, os.path.join(GFDL_BASE, 'exp', 'test_cases', 'frierson'))
from frierson_test_case import exp as frierson_exp

NCORES = 8
NCORES_OCEAN = int(os.environ.get('NCORES_OCEAN', 8))
NMONTHS = int(os.environ.get('NMONTHS', 1))
DAYS = int(os.environ.get('DAYS', 30))
DT_ATMOS = 720
COUPLING_PERIOD = 3600
MOM6_EXE = os.environ.get('MOM6_ISCA_EXE', os.path.expanduser(
    '~/project_specific_codes/isca_ocean_coupling/mom6_isca/build/mom6_isca.x'))
EXP_NAME = os.environ.get('EXP_NAME', 'coupled_ocean_mom6_aquaplanet')

cb = IscaCodeBase.from_directory(GFDL_BASE)

ocean = MOM6Ocean(MOM6_EXE, num_cores=NCORES_OCEAN, max_lat=80.)
exp = CoupledExperiment(EXP_NAME, codebase=cb, ocean=ocean, coupling_period=COUPLING_PERIOD)

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
for mod, field in [('dynamics', 'ps'), ('dynamics', 'bk'), ('dynamics', 'pk'),
                   ('atmosphere', 'precipitation'), ('mixed_layer', 't_surf'),
                   ('mixed_layer', 'flux_t'), ('mixed_layer', 'flux_lhe'),
                   ('mixed_layer', 'coupled_ocean_frac'), ('two_stream', 'olr'),
                   ('two_stream', 'swdn_toa'), ('two_stream', 'swdn_sfc'), ('two_stream', 'lwdn_sfc'),
                   ('two_stream', 'lwup_sfc'), ('dynamics', 'sphum'), ('dynamics', 'ucomp'),
                   ('dynamics', 'vcomp'), ('dynamics', 'temp'), ('dynamics', 'height')]:
    diag.add_field(mod, field, time_avg=True)
exp.diag_table = diag

exp.namelist = frierson_exp.namelist.copy()
exp.namelist['main_nml']['days'] = DAYS
exp.namelist['main_nml']['dt_atmos'] = DT_ATMOS
# slab depth used only for the polar caps outside the MOM6 domain
exp.namelist['mixed_layer_nml']['depth'] = 20.

if __name__ == '__main__':
    cb.compile()
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2, NMONTHS + 1):
        exp.run(i, num_cores=NCORES, overwrite_data=True)

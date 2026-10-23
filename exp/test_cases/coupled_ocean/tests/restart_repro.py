"""Restart reproducibility: one run of 2N days vs two runs of N days, for
uncoupled Isca, Isca + python slab, and Isca + MOM6. Prints which restart
variables differ. Usage: python restart_repro.py {uncoupled|slab|mom6} N"""
import os, sys, tarfile, io
import numpy as np, netCDF4
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from isca import Experiment, DiagTable, GFDL_DATA
from isca.coupling import CoupledExperiment, SlabOcean
from isca.coupling.mom6 import MOM6Ocean
import frierson_slab_oasis as f

kind, N = sys.argv[1], int(sys.argv[2])

def make(name, days):
    if kind == 'uncoupled':
        e = Experiment(name, codebase=f.cb)
    elif kind == 'slab':
        e = CoupledExperiment(name, codebase=f.cb, ocean=SlabOcean(depth=2.5), coupling_period=f.DT_ATMOS)
    else:
        e = CoupledExperiment(name, codebase=f.cb, coupling_period=3600, ocean=MOM6Ocean(
            os.path.expanduser('~/project_specific_codes/isca_ocean_coupling/mom6_isca/build/mom6_isca.x'),
            num_cores=4))
    e.diag_table = DiagTable(); e.diag_table.add_file('atmos_daily', 1, 'days', time_units='days')
    e.diag_table.add_field('mixed_layer', 't_surf', time_avg=True)
    e.namelist = f.make_namelist(days)
    return e

a = make('rr_%s_A' % kind, 2 * N); a.run(1, use_restart=False, num_cores=8, overwrite_data=True, mpirun_opts='--bind-to none')
b = make('rr_%s_B' % kind, N)
b.run(1, use_restart=False, num_cores=8, overwrite_data=True, mpirun_opts='--bind-to none')
b.run(2, num_cores=8, overwrite_data=True, mpirun_opts='--bind-to none')

def load(path):
    out = {}
    with tarfile.open(path) as t:
        for m in t.getmembers():
            if m.isfile() and m.name.endswith('.nc'):
                ds = netCDF4.Dataset('x', memory=t.extractfile(m).read())
                for v in ds.variables:
                    x = np.asarray(ds[v][:])
                    if x.dtype.kind == 'f':
                        out[m.name + ':' + v] = x
    return out
ra = load(os.path.join(GFDL_DATA, 'rr_%s_A' % kind, 'restarts', 'res0001.tar.gz'))
rb = load(os.path.join(GFDL_DATA, 'rr_%s_B' % kind, 'restarts', 'res0002.tar.gz'))
bad = [(k, float(np.nanmax(abs(ra[k] - rb[k])))) for k in ra if k in rb and not np.array_equal(ra[k], rb[k], equal_nan=True)]
print('RESULT %s: %d of %d restart variables differ' % (kind, len(bad), len(ra)))
for k, d in sorted(bad, key=lambda x: -x[1])[:12]:
    print('   %-50s %.3e' % (k, d))

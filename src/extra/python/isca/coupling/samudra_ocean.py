#!/usr/bin/env python3
"""The Samudra ML ocean emulator (M2LiNES) as an OASIS3-MCT component coupled to Isca.

Samudra steps the ocean in 5-day means on a 1-degree grid. A checkpoint with
2 input steps and 1 output step, forced by (tauuo, tauvo, hfds), predicts the
state of interval n from the states of intervals n-2, n-1 and the forcing over
those same intervals. Coupled to Isca:

  * at date d_n (start of interval n) Samudra receives Isca's fluxes averaged
    by OASIS over interval n-1 (LOCTRANS AVERAGE, LAG = one atmospheric step),
  * predicts interval n and immediately sends its top-level temperature as the
    SST Isca uses during interval n (LAG = 0: the coupling is sequential, Isca
    waits a few seconds every 5 model days; no lagged SST).

The initial ocean state and the forcing of the two intervals before the start
come from the OM4 data Samudra was trained on (--init-date). The emulator's
state (two normalised states and the forcing of the older interval) is saved
in --state-file at the end of each run; the forcing of the last interval of a
run reaches the next run through Isca's OASIS restart file.

Caveats (see PLAN.md): Earth geometry and OM4's present-day climate only; no
sea ice (OM4's hfds under ice is the ice-ocean flux, Isca sends the air-sea
flux); no freshwater forcing, so Isca's runoff/P-E do not affect salinity.

Run with the Python environment that has Samudra and torch, e.g.
    mpirun -np 16 isca.x : -np 1 python samudra_ocean.py --config eval.yaml ...
"""

import argparse
import os

import numpy as np
from mpi4py import MPI
import pyoasis
from pyoasis import OASIS

HLV = 2.500e6
FLUX_FIELDS = ['OCN_TAUX', 'OCN_TAUY', 'OCN_SWNET', 'OCN_LWNET', 'OCN_SENS', 'OCN_EVAP',
               'OCN_PRECIP', 'OCN_PSURF', 'OCN_RUNOFF']
GRID_NAME = 'smdr'


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--config', required=True, help='Samudra eval-style YAML (data + model)')
    p.add_argument('--ckpt', required=True)
    p.add_argument('--runtime', type=int, required=True, help='run length [s]')
    p.add_argument('--coupling-period', type=int, default=5 * 86400)
    p.add_argument('--state-file', default='samudra_state.pt')
    p.add_argument('--threads', type=int, default=8)
    args = p.parse_args()

    comp = pyoasis.Component('samudra', True, MPI.COMM_WORLD)

    import torch
    from samudra.config import EvalConfig
    from samudra.datasets import InferenceDataset
    from samudra.utils.data import BatchPreprocessor
    torch.set_num_threads(args.threads)

    # ---- model, normalisation and initial data, as Samudra's own Eval builds them
    cfg = EvalConfig.from_yaml(args.config)
    bundle = cfg.data.build(cfg.experiment.resolved_data_root, include_inference=True)
    layout = bundle.data_layout
    source = bundle.inference_source
    prog_names = list(layout.prognostic_var_names)
    bnd_names = list(layout.boundary_var_names)
    if bnd_names != ['tauuo', 'tauvo', 'hfds']:
        raise SystemExit('samudra_ocean.py needs a checkpoint forced by tauuo, tauvo, hfds; got %s'
                         % bnd_names)
    n_in, n_out = cfg.data.input_steps, cfg.data.output_steps
    if (n_in, n_out) != (2, 1):
        raise SystemExit('samudra_ocean.py supports 2 input / 1 output steps only')
    pre = BatchPreprocessor(source, prognostic_var_names=layout.prognostic_var_names,
                            boundary_var_names=layout.boundary_var_names)
    model = cfg.model.build(prog_channels=n_in * len(prog_names),
                            boundary_channels=n_in * len(bnd_names),
                            out_channels=n_out * len(prog_names), input_steps=n_in,
                            grid_sizes=[source.grid_size])
    ckpt = torch.load(args.ckpt, map_location='cpu', weights_only=False)
    model.load_state_dict({k.removeprefix('module.'): v for k, v in ckpt['model'].items()})
    model.eval()
    dataset = InferenceDataset(source=source, prognostic_var_names=layout.prognostic_var_names,
                               boundary_var_names=layout.boundary_var_names, input_steps=n_in,
                               output_steps=n_out, normalize_before_mask=cfg.data.normalize_before_mask,
                               masked_fill_value=cfg.data.masked_fill_value, long_rollout=True)
    i_sst = prog_names.index('thetao_0')

    # ---- the 1-degree grid and surface wet mask, for OASIS
    coords = source.to_xarray_dataset()

    def as2d(name_y, name_x):
        # coordinates may be stored 2-D (y, x) or as 1-D axes
        a, b = (np.asarray(coords[n].values, dtype=np.float64) for n in (name_y, name_x))
        if a.ndim == 1:
            b, a = np.meshgrid(b, a)
        return a, b

    lat_b, lon_b = as2d('lat_b', 'lon_b')                          # (ny+1, nx+1)
    lat, lon = as2d('lat', 'lon')                                  # (ny, nx)
    area = np.asarray(coords['areacello'].values, dtype=np.float64)
    wet = np.asarray(coords['mask_0'].values, dtype=bool) if 'mask_0' in coords else (
        np.isfinite(area) & (np.nan_to_num(area) > 0))
    ny, nx = lat.shape
    npts = nx * ny

    def F(a2d):            # (ny, nx) C-order -> OASIS (nx, ny) Fortran-order flat
        return np.asarray(a2d).T.ravel(order='F').copy() if a2d.ndim == 2 else a2d

    def unF(flat):         # OASIS flat -> (ny, nx)
        return np.asarray(flat).reshape((nx, ny), order='F').T

    partition = pyoasis.SerialPartition(npts)
    grid = pyoasis.Grid(GRID_NAME, nx, ny, lon.T.copy(), lat.T.copy(), partition)
    clo = np.stack([lon_b[:-1, :-1], lon_b[:-1, 1:], lon_b[1:, 1:], lon_b[1:, :-1]], -1)
    cla = np.stack([lat_b[:-1, :-1], lat_b[:-1, 1:], lat_b[1:, 1:], lat_b[1:, :-1]], -1)
    grid.set_corners(pyoasis.asarray(clo.transpose(1, 0, 2).copy()),
                     pyoasis.asarray(cla.transpose(1, 0, 2).copy()))
    grid.set_mask(np.where(wet, 0, 1).astype(np.int32).T.copy())
    grid.set_area(np.nan_to_num(area).T.copy())
    grid.write()

    with open('namcouple') as f:
        words = set(f.read().split())
    flux_vars = {n: pyoasis.Var(n, partition, OASIS.IN) for n in FLUX_FIELDS if n in words}
    sst_var = pyoasis.Var('SST_OCN', partition, OASIS.OUT)
    frac_var = pyoasis.Var('FRAC_OCN', partition, OASIS.OUT)
    comp.enddef()

    # ---- state: prognostic history (1, 2*77, ny, nx) and the older interval's forcing (raw)
    if os.path.exists(args.state_file):
        st = torch.load(args.state_file, weights_only=False)
        prog, older_forcing = st['prog'], st['older_forcing']
        print('samudra_ocean: restarting from %s' % args.state_file, flush=True)
        cold = False
    else:
        prog = dataset.initial_prognostic
        # raw (unnormalised) forcing of the two initial intervals from OM4
        raw0 = pre.unnormalize_tensor_boundary(
            dataset.get_boundary(0).reshape(1, n_in, len(bnd_names), ny, nx), fill_value=0.)
        older_forcing, newer_forcing_om4 = raw0[:, 0], raw0[:, 1]
        cold = True
        print('samudra_ocean: cold start from OM4 at %s' % str(dataset.get_target_time(0, 1)[0]),
              flush=True)

    fields = {n: pyoasis.asarray(np.zeros(npts)) for n in flux_vars}
    dt = args.coupling_period
    with torch.no_grad():
        for date in range(0, args.runtime, dt):
            for n, var in flux_vars.items():
                var.get(date, fields[n])
            if date == 0 and cold:
                newer = newer_forcing_om4                 # Isca has not run yet
            else:
                hfds = (fields['OCN_SWNET'] + fields['OCN_LWNET'] - fields['OCN_SENS']
                        - HLV * fields['OCN_EVAP'])
                newer = torch.from_numpy(np.stack([unF(fields['OCN_TAUX']), unF(fields['OCN_TAUY']),
                                                   unF(hfds)])[None].astype(np.float32))
            raw = torch.stack([older_forcing, newer], dim=1)          # (1, 2, 3, ny, nx)
            boundary = pre.prepare_boundary(raw, torch.device('cpu'))
            out = model.forward_once(prog, boundary, dataset.ctx)
            pred = model._assemble_prediction(prog, out)
            prog = model._advance_prognostic_history(prog, pred)
            older_forcing = newer

            phys = pre.unnormalize_tensor_prognostic(pred[:, None], fill_value=float('nan'))[0, 0]
            sst = phys[i_sst].numpy() + 273.15
            sst = np.where(wet & np.isfinite(sst), sst, 0.)
            sst_var.put(date, pyoasis.asarray(F(sst)))
            frac_var.put(date, pyoasis.asarray(F(wet.astype(np.float64))))
            if date % (6 * dt) == 0:
                print('samudra_ocean: day %d, mean SST %.3f C' %
                      (date // 86400, float(np.mean(sst[wet]) - 273.15)), flush=True)

    torch.save({'prog': prog, 'older_forcing': older_forcing}, args.state_file)
    print('samudra_ocean: finished, wrote %s' % args.state_file, flush=True)
    del comp


if __name__ == '__main__':
    main()

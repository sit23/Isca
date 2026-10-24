"""Is the two-step (10-day) SST oscillation of coupled Isca+Samudra intrinsic to the emulator?

Restarts Samudra offline from samudra_state.pt files saved by a coupled run, holds the
forcing fixed at the last forcing Isca gave it (both input slots), and rolls out NSTEPS
5-day steps with the same temperature bounds as samudra_ocean.py. Reports the global-mean
SST and the amplitude of its odd-even component, |SST_n - (SST_n-1 + SST_n+1)/2|. Growth
under fixed forcing means an emulator instability; decay means the atmosphere drives it.

With RA_NU set, a Robert-Asselin filter of strength RA_NU is applied to the older state of
the history after every step (x_n += nu * (x_n-1 - 2 x_n + x_n+1)), to damp the 2-step mode.

Usage (Samudra environment): [RA_NU=0.1] python samudra_mode_test.py eval.yaml ckpt.pt NSTEPS state.pt [state.pt ...]
"""
import os
import sys

import numpy as np
import torch
from samudra.config import EvalConfig
from samudra.datasets import InferenceDataset
from samudra.utils.data import BatchPreprocessor

torch.set_num_threads(8)
cfg = EvalConfig.from_yaml(sys.argv[1])
nsteps = int(sys.argv[3])
ra_nu = float(os.environ.get('RA_NU', 0.))
bundle = cfg.data.build(cfg.experiment.resolved_data_root, include_inference=True)
layout, source = bundle.data_layout, bundle.inference_source
prog_names = list(layout.prognostic_var_names)
pre = BatchPreprocessor(source, prognostic_var_names=layout.prognostic_var_names,
                        boundary_var_names=layout.boundary_var_names)
model = cfg.model.build(prog_channels=2 * len(prog_names), boundary_channels=6,
                        out_channels=len(prog_names), input_steps=2, grid_sizes=[source.grid_size])
ckpt = torch.load(sys.argv[2], map_location='cpu', weights_only=False)
model.load_state_dict({k.removeprefix('module.'): v for k, v in ckpt['model'].items()})
model.eval()
ds = InferenceDataset(source=source, prognostic_var_names=layout.prognostic_var_names,
                      boundary_var_names=layout.boundary_var_names, input_steps=2, output_steps=1,
                      normalize_before_mask=cfg.data.normalize_before_mask,
                      masked_fill_value=cfg.data.masked_fill_value, long_rollout=True)
i_sst = prog_names.index('thetao_0')
i_temp = [k for k, n in enumerate(prog_names) if n.startswith('thetao_')]
t_mean = torch.as_tensor(pre._prognostic_mean_np, dtype=torch.float32)[i_temp]
t_std = torch.as_tensor(pre._prognostic_std_np, dtype=torch.float32)[i_temp]
t_lo = ((-2. - t_mean) / t_std).reshape(-1, *([1] * (3 - t_mean.dim())))
t_hi = ((35. - t_mean) / t_std).reshape(-1, *([1] * (3 - t_mean.dim())))
t_wet = torch.as_tensor(np.asarray(pre.prognostic_mask))[i_temp]
wet = np.asarray(pre.prognostic_mask)[i_sst]
area = np.nan_to_num(np.asarray(source.to_xarray_dataset()['areacello'].values)) * wet


def gm(prog_slot):
    sst = pre.unnormalize_tensor_prognostic(prog_slot[:, None], fill_value=float('nan'))[0, 0, i_sst].numpy()
    return float(np.nansum(sst * area) / area.sum()), sst


for fn in sys.argv[4:]:
    st = torch.load(fn, weights_only=False)
    prog, forcing = st['prog'], st['older_forcing']
    boundary = pre.prepare_boundary(torch.stack([forcing, forcing], 1), torch.device('cpu'))
    n = len(prog_names)
    series = [gm(prog[:, :n])[0], gm(prog[:, n:])[0]]
    maps = [gm(prog[:, n:])[1]]
    with torch.no_grad():
        for k in range(nsteps):
            pred = model._assemble_prediction(prog, model.forward_once(prog, boundary, ds.ctx))
            temp = pred[0, i_temp]
            pred[0, i_temp] = torch.where(t_wet, torch.minimum(torch.maximum(temp, t_lo), t_hi), temp)
            if ra_nu:
                # filter interval n (the newer input) with n-1 and the prediction n+1
                prog[:, n:] += ra_nu * (prog[:, :n] - 2 * prog[:, n:] + pred)
            prog = model._advance_prognostic_history(prog, pred)
            g, m = gm(pred)
            series.append(g); maps.append(m)
    s = np.array(series)
    zig = np.abs(s[1:-1] - 0.5 * (s[:-2] + s[2:]))
    mp = np.array(maps)
    zmap = np.abs(mp[1:-1] - 0.5 * (mp[:-2] + mp[2:]))
    rms = np.sqrt(np.nanmean(zmap**2, axis=(1, 2)))
    print(fn)
    print('  global SST', ' '.join('%.2f' % v for v in s[::4]))
    print('  odd-even global amp', ' '.join('%.2f' % v for v in zig[::4]))
    print('  odd-even local rms ', ' '.join('%.2f' % v for v in rms[::4]))

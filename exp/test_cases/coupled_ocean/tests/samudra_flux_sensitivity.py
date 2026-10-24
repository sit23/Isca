"""How strongly does Samudra's SST respond to the heat flux it is forced with?

Starts from the OM4 state (as a cold-start coupled run does), perturbs hfds of the
newer forcing interval by +/-DF W/m2 over the whole ocean, and reports the SST
response after one 5-day step as an effective mixed-layer depth
h = dt * dF / (rho * c_p * dSST). Explicit 5-day coupling is unstable where
dt * lambda / (rho c_p h) > 2, lambda being the atmosphere's dF/dSST.

Usage (Samudra environment): python samudra_flux_sensitivity.py eval.yaml ckpt.pt [DF]
"""
import sys

import numpy as np
import torch
from samudra.config import EvalConfig
from samudra.datasets import InferenceDataset
from samudra.utils.data import BatchPreprocessor

RHO_CP = 1025. * 3989.
DT = 5 * 86400.
torch.set_num_threads(8)
cfg = EvalConfig.from_yaml(sys.argv[1])
dF = float(sys.argv[3]) if len(sys.argv) > 3 else 50.
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
ny, nx = source.grid_size
prog0 = ds.initial_prognostic
raw0 = pre.unnormalize_tensor_boundary(ds.get_boundary(0).reshape(1, 2, 3, ny, nx), fill_value=0.)
i_sst = prog_names.index('thetao_0')
wet = np.asarray(pre.prognostic_mask)[i_sst]
lat = np.asarray(source.to_xarray_dataset()['lat'].values)
lat2 = lat if lat.ndim == 2 else np.repeat(lat[:, None], nx, 1)


def step(prog, older, newer):
    with torch.no_grad():
        b = pre.prepare_boundary(torch.stack([older, newer], 1), torch.device('cpu'))
        pred = model._assemble_prediction(prog, model.forward_once(prog, b, ds.ctx))
    sst = pre.unnormalize_tensor_prognostic(pred[:, None], fill_value=float('nan'))[0, 0, i_sst].numpy()
    return pred, sst


older, newer = raw0[:, 0], raw0[:, 1]
res = {}
for sgn in (-1., 0., 1.):
    nw = newer.clone()
    nw[:, 2] += sgn * dF * torch.from_numpy(wet.astype(np.float32))
    pred1, sst1 = step(prog0, older, nw)
    # second step: the perturbed interval is now the "older" forcing, the new one unperturbed
    prog1 = model._advance_prognostic_history(prog0, pred1)
    _, sst2 = step(prog1, nw, newer)
    res[sgn] = (sst1, sst2)

for k, name in ((0, 'step 1 (perturbed interval is the newer forcing)'),
                (1, 'step 2 (perturbed interval is the older forcing)')):
    dT = (res[1.][k] - res[-1.][k]) / 2.
    print(name)
    for lo, hi in ((-70, -40), (-40, -15), (-15, 15), (15, 40), (40, 70)):
        m = wet & (lat2 >= lo) & (lat2 < hi) & np.isfinite(dT)
        d = np.median(dT[m])
        print('  lat %4d..%3d: median dSST %+.3f K per %g W/m2 -> effective depth %6.1f m '
              '(p10 %.3f, p90 %.3f)' % (lo, hi, d, dF, DT * dF / (RHO_CP * d) if d > 0 else np.inf,
                                        np.percentile(dT[m], 10), np.percentile(dT[m], 90)))

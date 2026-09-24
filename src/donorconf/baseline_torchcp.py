"""Published-package baseline: the conformal engine of the conformalized single-cell annotator.

López-De-Castro et al. (Bioinformatics 2025, btaf521) calibrate with torchCP 1.0.2 `SplitPredictor`,
`ClassWisePredictor` and `ClusteredPredictor` on a held-out set of cells. This module runs exactly those torchCP
classes on the head's probabilities so the comparison isolates the CALIBRATION step. It does NOT include the
annotator's out-of-distribution detector or its neural classifier; say so wherever results are reported.

torchCP imports `torchsort` at package import only for a training loss that calibration never uses; when torchsort
is not installed a stub module is inserted so the import succeeds.
"""
from __future__ import annotations

import sys
import types

import numpy as np


def _import_torchcp():
    if "torchsort" not in sys.modules:
        try:
            import torchsort  # noqa: F401
        except Exception:
            m = types.ModuleType("torchsort")
            m.soft_rank = m.soft_sort = None
            sys.modules["torchsort"] = m
    import torch
    from torchcp.classification.predictor import ClassWisePredictor, ClusteredPredictor, SplitPredictor
    from torchcp.classification.score import APS, THR
    return torch, SplitPredictor, ClassWisePredictor, ClusteredPredictor, THR, APS


VARIANTS = {"tcp_standard_thr": ("standard", "THR"), "tcp_standard_aps": ("standard", "APS"),
            "tcp_classwise_thr": ("classwise", "THR"), "tcp_cluster_thr": ("cluster", "THR")}


def torchcp_sets(P_cal, y_cal, P_tgt, alpha: float, variant: str, seed: int = 0) -> np.ndarray:
    """Boolean (n_target, n_classes) prediction sets from a torchCP predictor calibrated on (P_cal, y_cal)."""
    torch, Split, Class, Clust, THR, APS = _import_torchcp()
    kind, score = VARIANTS[variant]
    torch.manual_seed(seed)
    logp = lambda P: torch.log(torch.as_tensor(np.clip(P, 1e-12, 1.0), dtype=torch.float32))
    ds = torch.utils.data.TensorDataset(logp(P_cal), torch.as_tensor(y_cal, dtype=torch.long))
    loader = torch.utils.data.DataLoader(ds, batch_size=2048)
    sf = THR() if score == "THR" else APS()

    class _Pass(torch.nn.Module):          # torchCP reads model.parameters() to find the device; carry one dummy
        def __init__(self):
            super().__init__()
            self.dummy = torch.nn.Parameter(torch.zeros(1), requires_grad=False)

        def forward(self, x):
            return x

    ident = _Pass()
    if kind == "standard":
        pred = Split(sf, ident)
    elif kind == "classwise":
        pred = Class(sf, ident)
    else:
        pred = Clust(sf, ident, num_clusters=max(2, P_cal.shape[1] // 2))
    pred.calibrate(loader, alpha)
    out = pred.predict(logp(P_tgt))
    if hasattr(out, "cpu"):                 # torchCP returns a 0/1 integer tensor of shape (n, k)
        return out.cpu().numpy().astype(bool)
    n, k = P_tgt.shape
    S = np.zeros((n, k), dtype=bool)
    for i, s in enumerate(out):
        S[i, list(np.asarray(s, dtype=int))] = True
    return S

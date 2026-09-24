"""Turn a labelled embedding table into the .npz inputs donorconf expects.

Needs REAL data you have obtained and reviewed. Nothing here invents labels.

Input: an .h5ad (requires anndata) or a .npz with arrays X (cells x dims), labels
(cell-type strings), donors, dataset. Reference donors train the classification
head; calibration donors calibrate; target donors are held out.

    python scripts/prepare_inputs.py --input atlas.h5ad --embedding X_scGPT \
        --label-key cell_type --donor-key donor_id --dataset-key dataset \
        --reference DS_A --calibration DS_B --target DS_C --out data/run1

Donor ids are namespaced with the dataset so ids from different studies never
collide. Cell types not present in the reference dataset are dropped from
calibration/target and counted in the report, never silently relabelled.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from donorconf.simulate import Cohort, fit_head, head_probs  # noqa: E402


def load(a):
    if a.input.endswith(".h5ad"):
        try:
            import anndata
        except ImportError:
            raise SystemExit("pip install anndata to read .h5ad")
        ad = anndata.read_h5ad(a.input)
        X = np.asarray(ad.obsm[a.embedding])
        return X, ad.obs[a.label_key].astype(str).to_numpy(), \
            ad.obs[a.donor_key].astype(str).to_numpy(), ad.obs[a.dataset_key].astype(str).to_numpy()
    z = np.load(a.input, allow_pickle=False)
    return z["X"], z["labels"].astype(str), z["donors"].astype(str), z["dataset"].astype(str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--embedding", default="X")
    ap.add_argument("--label-key", default="cell_type")
    ap.add_argument("--donor-key", default="donor_id")
    ap.add_argument("--dataset-key", default="dataset")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--calibration", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    X, lab, don, ds = load(a)
    don = np.char.add(np.char.add(ds, "::"), don)
    ref_m = ds == a.reference
    classes = np.unique(lab[ref_m])
    cmap = {c: i for i, c in enumerate(classes)}
    report = {"classes": classes.tolist(), "dropped_unseen_type_cells": {}}

    def cohort(name):
        m = ds == name
        keep = m & np.isin(lab, classes)
        report["dropped_unseen_type_cells"][name] = int(m.sum() - keep.sum())
        return Cohort(X[keep].astype(float), np.array([cmap[c] for c in lab[keep]]), don[keep])

    ref, cal, tgt = cohort(a.reference), cohort(a.calibration), cohort(a.target)
    for nm, c in (("reference", ref), ("calibration", cal), ("target", tgt)):
        report[f"n_donors_{nm}"] = c.n_donors
        report[f"n_cells_{nm}"] = int(c.y.size)
    head = fit_head(ref)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    K = len(classes)
    np.savez_compressed(out / "calibration.npz", probs=head_probs(head, cal.X, K), labels=cal.y, donors=cal.donor)
    np.savez_compressed(out / "target.npz", probs=head_probs(head, tgt.X, K), labels=tgt.y, donors=tgt.donor)
    (out / "prepare_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

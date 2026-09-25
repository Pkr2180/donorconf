"""Merge Geneformer/scGPT embeddings (fetched from the Modal volume) into a dataset's npz.

Usage (from donorconf/):
    python -m modal volume get donorconf-data <tag>.npz results/real/data/<tag>.npz
    python -m modal volume get donorconf-data <tag>_geneformer.npy /tmp/<tag>_geneformer.npy
    python -m modal volume get donorconf-data <tag>_scgpt.npy /tmp/<tag>_scgpt.npy
    python scripts/merge_fm_embeddings.py <tag>

Adds `emb_geneformer` and `emb_scgpt` keys to results/real/data/<tag>.npz (which already has
emb_tf-sapiens, emb_tf-exemplar-human, emb_scvi, labels, donors, disease, tissue, sex from
modal_census.py::pull_raw). Refuses to merge if a new embedding's row count doesn't match the
existing arrays, or if it contains any non-finite values (silent corruption would otherwise
propagate straight into the benchmark).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main(tag: str):
    npz_path = ROOT / "results" / "real" / "data" / f"{tag}.npz"
    gf_path = Path(f"/tmp/{tag}_geneformer.npy")
    sc_path = Path(f"/tmp/{tag}_scgpt.npy")

    d = dict(np.load(npz_path, allow_pickle=False))
    n_cells = d["labels"].shape[0]

    added = []
    for name, path in (("emb_geneformer", gf_path), ("emb_scgpt", sc_path)):
        if not path.exists():
            print(f"SKIP {name}: {path} not found")
            continue
        arr = np.load(path)
        if arr.shape[0] != n_cells:
            raise ValueError(f"{name}: {arr.shape[0]} rows but {tag}.npz has {n_cells} cells")
        n_bad = (~np.isfinite(arr)).any(axis=1).sum()
        if n_bad:
            raise ValueError(f"{name}: {n_bad}/{n_cells} rows have non-finite values -- investigate before merging")
        # a degenerate all-zero row is finite (passes the check above) but is not a real embedding --
        # seen once on lung_I's Geneformer output (27% of cells, cause unresolved as of 2026-09-25).
        n_zero = (np.abs(arr).sum(axis=1) == 0).sum()
        if n_zero:
            raise ValueError(f"{name}: {n_zero}/{n_cells} rows are all-zero (degenerate, not real embeddings) -- investigate before merging")
        d[name] = arr.astype(np.float32)
        added.append(f"{name} (dim={arr.shape[1]})")

    if not added:
        print(f"Nothing to merge for {tag} -- no geneformer/scgpt .npy found at the expected /tmp paths")
        return

    np.savez_compressed(npz_path, **d)
    print(f"{tag}.npz updated: added {', '.join(added)}. Keys now: {sorted(d.keys())}")


if __name__ == "__main__":
    main(sys.argv[1])

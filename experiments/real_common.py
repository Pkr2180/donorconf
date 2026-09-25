"""Shared helpers for the real-data (CELLxGENE Census) experiments."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from _common import ROOT  # noqa: F401  (also puts src/ on sys.path)
from donorconf import CellLevelConformal, ClasswiseConformal, DonorLevelCRC, DonorWeightedQuantile
from donorconf.head import EmbeddingHead
from donorconf.metrics import coverage_summary
from donorconf.shift import donor_summaries, shift_screen

DATA = ROOT / "results" / "real" / "data"
EMBEDDINGS = ("emb_tf-sapiens", "emb_tf-exemplar-human", "emb_scvi", "emb_geneformer", "emb_scgpt")


def load(tag: str) -> dict:
    z = np.load(DATA / f"{tag}.npz", allow_pickle=False)
    return {k: z[k] for k in z.files}


def choose_classes(labels, donors, min_cells: int, min_donors: int):
    """Label vocabulary fixed by frequency; independent of any donor split."""
    types = np.unique(labels)
    keep = []
    for t in types:
        m = labels == t
        if m.sum() >= min_cells and np.unique(donors[m]).size >= min_donors:
            keep.append(t)
    return np.array(keep)


def encode(labels, classes):
    idx = {c: i for i, c in enumerate(classes)}
    y = np.array([idx.get(l, -1) for l in labels])
    return y


def evaluate_split(X, y, donors, ref_d, cal_d, tgt_d, alpha, n_classes, seed, pca=64, torchcp=False):
    """Fit head on reference donors, calibrate on calibration donors, score target donors."""
    m_ref, m_cal, m_tgt = (np.isin(donors, s) for s in (ref_d, cal_d, tgt_d))
    head = EmbeddingHead(n_components=pca, seed=seed).fit(X[m_ref], y[m_ref], n_classes)
    F_cal, F_tgt = head.transform(X[m_cal]), head.transform(X[m_tgt])
    P_cal, P_tgt = head.predict_proba(None, F_cal), head.predict_proba(None, F_tgt)
    y_cal, y_tgt, d_cal, d_tgt = y[m_cal], y[m_tgt], donors[m_cal], donors[m_tgt]

    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(y_cal))
    a, b = perm[: len(perm) // 2], perm[len(perm) // 2:]
    leak = CellLevelConformal(alpha).fit(P_cal[a], y_cal[a])
    reported = float(leak.predict_sets(P_cal[b])[np.arange(len(b)), y_cal[b]].mean())

    rows = []
    for name, model in (("cell_pooled", CellLevelConformal(alpha)),
                        ("classwise", ClasswiseConformal(alpha)),
                        ("donor_weighted", DonorWeightedQuantile(alpha)),
                        ("donor_crc", DonorLevelCRC(alpha))):
        model.fit(P_cal, y_cal, d_cal)
        cs = coverage_summary(model.predict_sets(P_tgt), y_tgt, d_tgt, alpha)
        rows.append({"method": name, **cs,
                     "reported_random_split_coverage": reported if name == "cell_pooled" else np.nan})
    if torchcp:
        from donorconf.baseline_torchcp import VARIANTS, torchcp_sets
        for name in VARIANTS:
            try:
                sets = torchcp_sets(P_cal, y_cal, P_tgt, alpha, name, seed=seed)
            except Exception as e:          # e.g. clustered calibration with too few cells for a class
                print("torchcp", name, "failed:", repr(e)[:80], flush=True)
                continue
            rows.append({"method": name, **coverage_summary(sets, y_tgt, d_tgt, alpha),
                         "reported_random_split_coverage": np.nan})
    centroid = F_cal.mean(axis=0)
    _, sc, _ = donor_summaries(P_cal, d_cal, F_cal, centroid)
    _, st, _ = donor_summaries(P_tgt, d_tgt, F_tgt, centroid)
    scr = shift_screen(sc, st, level=0.05, n_perm=499, max_exact=0, seed=seed)
    head_acc = float((P_tgt.argmax(axis=1) == y_tgt).mean())
    return rows, {"screen_p": scr["p_value"], "screen_flagged": scr["flagged"], "head_accuracy_target": head_acc,
                  "head_accuracy_calibration": float((P_cal.argmax(axis=1) == y_cal).mean()),
                  "conf_shift": float(P_tgt.max(axis=1).mean() - P_cal.max(axis=1).mean()),
                  "entropy_shift": float(-(P_tgt * np.log(np.clip(P_tgt, 1e-12, 1))).sum(axis=1).mean()
                                         + (P_cal * np.log(np.clip(P_cal, 1e-12, 1))).sum(axis=1).mean())}


def save_meta(path: Path, obj: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))


def prepare_pair(sets, scheme, min_cells, min_donors, emb):
    """Concatenate several Census datasets on a shared coarse label vocabulary, loading ONE embedding.

    Returns (X, y, donors, dataset_of_cell, classes, info). Cells whose label maps to '' or falls outside the
    shared vocabulary are dropped and counted in info['per_set'].
    """
    from donorconf.lineage import map_labels
    zs = {t: np.load(DATA / f"{t}.npz", allow_pickle=False) for t in sets}
    lab = {t: map_labels(zs[t]["labels"], scheme) for t in sets}
    don = {t: zs[t]["donors"] for t in sets}
    vocab = None
    for t in sets:
        ok = set()
        for c in np.unique(lab[t]):
            if c == "":
                continue
            m = lab[t] == c
            if m.sum() >= min_cells and np.unique(don[t][m]).size >= min_donors:
                ok.add(c)
        vocab = ok if vocab is None else vocab & ok
    classes = np.array(sorted(vocab))
    if classes.size < 3:
        raise ValueError(f"only {classes.size} shared types for {sets}")
    Xs, ys, ds, dsets, info = [], [], [], [], {"classes": classes.tolist(), "per_set": {}}
    for t in sets:
        y = encode(lab[t], classes)
        k = y >= 0
        info["per_set"][t] = {"cells": int(len(y)), "kept": int(k.sum()), "donors": int(np.unique(don[t][k]).size)}
        Xs.append(zs[t][emb][k]); ys.append(y[k]); ds.append(don[t][k]); dsets.append(np.full(k.sum(), t))
    return np.vstack(Xs), np.concatenate(ys), np.concatenate(ds), np.concatenate(dsets), classes, info

"""Audit lung_F..I against the fixed inclusion rule in ANALYSIS_LOCK.md Addendum C2, BEFORE any coverage is computed.

Rule: >= 15 donors with >= 20 mapped cells each, and >= 60% of all cells mapped (non-excluded) after lung
restriction. Appends per-label rows to config/lineage_audit_lung.csv (same format as lung_A..E) and prints the
per-dataset verdict. No hand-typed numbers: everything below is read from the npz files.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np

from donorconf.lineage import lung_coarse

R = Path(__file__).resolve().parents[1]
DATA = R / "results" / "real" / "data"
AUDIT = R / "config" / "lineage_audit_lung.csv"
TAGS = ["lung_F", "lung_G", "lung_H", "lung_I"]


def audit_one(tag: str):
    d = np.load(DATA / f"{tag}.npz", allow_pickle=True)
    labels, donors = d["labels"], d["donors"]
    n = len(labels)
    coarse = np.array([lung_coarse(str(s)) or "<excluded>" for s in labels])
    mapped = coarse != "<excluded>"
    frac_mapped = mapped.mean()

    per_donor_mapped = Counter(donors[mapped])
    donors_ok = sum(1 for c in per_donor_mapped.values() if c >= 20)

    qualifies = donors_ok >= 15 and frac_mapped >= 0.60
    rows = []
    label_counts = Counter(zip(labels.tolist(), coarse.tolist()))
    for (lab, co), cnt in sorted(label_counts.items()):
        rows.append({"dataset": tag, "source_label": lab, "coarse": co, "n_cells": cnt})
    return qualifies, frac_mapped, donors_ok, len(set(donors)), n, rows


def main():
    existing = AUDIT.read_text(encoding="utf-8").splitlines()
    header = existing[0]
    body = existing[1:]
    verdicts = []
    total_rows = 0
    for tag in TAGS:
        qualifies, frac_mapped, donors_ok, n_donors, n_cells, rows = audit_one(tag)
        verdicts.append((tag, qualifies, frac_mapped, donors_ok, n_donors, n_cells))
        body += [",".join(f'"{r[c]}"' if isinstance(r[c], str) and ("," in r[c]) else str(r[c])
                           for c in ("dataset", "source_label", "coarse", "n_cells")) for r in rows]
        total_rows += len(rows)
    AUDIT.write_text("\n".join([header] + body) + "\n", encoding="utf-8")

    print(f"{'tag':8s} {'qualifies':10s} {'frac_mapped':12s} {'donors>=20cells':16s} {'n_donors':9s} {'n_cells':8s}")
    for tag, q, fm, dok, nd, nc in verdicts:
        print(f"{tag:8s} {str(q):10s} {fm:12.3f} {dok:16d} {nd:9d} {nc:8d}")
    qualifying = [t for t, q, *_ in verdicts if q]
    print("\nQualifying lung_F..I datasets:", qualifying if qualifying else "NONE")
    print(f"appended {total_rows} label rows to {AUDIT}")


if __name__ == "__main__":
    main()

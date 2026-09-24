"""Command line: calibrate / predict / screen / report on .npz inputs.

Input arrays (float32/64 probs, int labels, donor ids of any dtype):
  calibration .npz : probs (cells x classes), labels, donors
  target      .npz : probs, donors            (labels only for `report`)
Probabilities come from any classification head on foundation-model embeddings;
this tool does not extract embeddings.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .conformal import DonorLevelCRC
from .metrics import coverage_summary, donor_bootstrap_coverage_ci
from .shift import donor_summaries, shift_screen


def _load(path, need):
    z = np.load(path, allow_pickle=False)
    missing = [k for k in need if k not in z.files]
    if missing:
        raise SystemExit(f"{path} is missing arrays: {missing}")
    return {k: z[k] for k in z.files}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="donorconf", description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("calibrate", help="fit donor-level CRC threshold")
    c.add_argument("--input", required=True)
    c.add_argument("--alpha", type=float, default=0.1)
    c.add_argument("--score", default="lac", choices=["lac", "aps"])
    c.add_argument("--out", required=True)

    p = sub.add_parser("predict", help="prediction sets for new cells")
    p.add_argument("--model", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)

    s = sub.add_parser("screen", help="label-free donor-level shift screen")
    s.add_argument("--cal", required=True)
    s.add_argument("--target", required=True)
    s.add_argument("--level", type=float, default=0.05)
    s.add_argument("--out", required=True)

    r = sub.add_parser("report", help="coverage on LABELED target donors")
    r.add_argument("--model", required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--out", required=True)

    a = ap.parse_args(argv)

    if a.cmd == "calibrate":
        d = _load(a.input, ["probs", "labels", "donors"])
        m = DonorLevelCRC(a.alpha, a.score).fit(d["probs"], d["labels"], d["donors"])
        json.dump(m.to_dict(), open(a.out, "w"), indent=2)
        if m.vacuous_:
            print(f"WARNING: alpha={a.alpha} < 1/(n+1)={m.min_alpha_:.3f}; sets contain all classes",
                  file=sys.stderr)
    elif a.cmd == "predict":
        m = DonorLevelCRC.from_dict(json.load(open(a.model)))
        d = _load(a.input, ["probs"])
        np.savez_compressed(a.out, sets=m.predict_sets(d["probs"]))
    elif a.cmd == "screen":
        cal = _load(a.cal, ["probs", "donors"])
        tgt = _load(a.target, ["probs", "donors"])
        _, sc, names = donor_summaries(cal["probs"], cal["donors"])
        _, st, _ = donor_summaries(tgt["probs"], tgt["donors"])
        res = shift_screen(sc, st, level=a.level)
        res["summary_columns"] = names
        json.dump(res, open(a.out, "w"), indent=2)
    elif a.cmd == "report":
        m = DonorLevelCRC.from_dict(json.load(open(a.model)))
        d = _load(a.input, ["probs", "labels", "donors"])
        sets = m.predict_sets(d["probs"])
        rep = coverage_summary(sets, d["labels"], d["donors"], m.alpha)
        mean, lo, hi = donor_bootstrap_coverage_ci(sets, d["labels"], d["donors"])
        rep.update({"donor_mean_coverage_ci95": [lo, hi], "target_coverage": 1 - m.alpha,
                    "note": "CI is over labeled target donors; it does not repair failed calibration."})
        json.dump(rep, open(a.out, "w"), indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

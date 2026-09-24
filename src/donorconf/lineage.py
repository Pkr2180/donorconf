"""Transparent coarse cell-type mapping so different studies share a label vocabulary.

Census `cell_type` names differ in granularity between studies ("CD4-positive, alpha-beta T
cell" vs "T cell"). Rules below are plain substring tests applied in order; labels that match
nothing map to None and are excluded (and counted by the caller). The full source-label ->
coarse-label table is exported by `audit_table` so every decision can be inspected.

These mappings are a methodological choice, not ground truth: they discard fine distinctions.
"""
from __future__ import annotations

import re

import numpy as np


def _has(s: str, *needles: str) -> bool:
    return any(n in s for n in needles)


_T = re.compile(r"\bt cell")      # word boundary: "mast cell" must NOT match
_B = re.compile(r"\bb cell")


def blood_fine(label: str):
    s = label.lower()
    if _has(s, "plasmacytoid"):
        return None
    if "natural killer cell" in s:
        return "NK"
    if _T.search(s) and "cd4-positive" in s:
        return "CD4_T"
    if "regulatory t cell" in s and "double negative" not in s:
        return "CD4_T"
    if _T.search(s) and "cd8-positive" in s:
        return "CD8_T"
    if _B.search(s):
        return "B"
    if _has(s, "non-classical monocyte", "cd14-low, cd16-positive monocyte"):
        return "mono_nonclassical"
    if _has(s, "cd14-positive monocyte", "classical monocyte"):
        return "mono_classical"
    if _has(s, "conventional dendritic cell", "cd1c-positive myeloid dendritic cell"):
        return "cDC"
    return None


def oral_coarse(label: str):
    s = label.lower()
    if _has(s, "plasmacytoid"):
        return None
    if _has(s, "nk t cell", "natural killer t"):
        return None
    if _has(s, "plasma cell", "plasmablast"):
        return "plasma"
    if "natural killer cell" in s:
        return "NK"
    if _B.search(s):
        return "B"
    if _T.search(s) or _has(s, "t-helper", "helper t"):
        return "T"
    if _has(s, "myeloid cell", "monocyte", "macrophage", "dendritic cell", "langerhans"):
        return "myeloid"
    if "mast cell" in s:
        return "mast"
    if _has(s, "epithelial cell", "keratinocyte", "squamous cell"):
        return "epithelial"
    if "fibroblast" in s:
        return "fibroblast"
    if "endothelial cell" in s:
        return "endothelial"
    if _has(s, "pericyte", "smooth muscle", "mural cell"):
        return "mural"
    return None


def lung_coarse(label: str):
    s = label.lower()
    if _has(s, "plasmacytoid", "nk t cell", "natural killer t"):
        return None
    if _has(s, "type ii pneumocyte", "type 2 pneumocyte", "alveolar type 2"):
        return "AT2"
    if _has(s, "type i pneumocyte", "type 1 pneumocyte", "alveolar type 1"):
        return "AT1"
    if _has(s, "macrophage"):
        return "macrophage"
    if _has(s, "monocyte"):
        return "monocyte"
    if "natural killer cell" in s:
        return "NK"
    if _B.search(s):
        return "B"
    if _T.search(s):
        return "T"
    if _has(s, "ciliated"):
        return "ciliated"
    if _has(s, "club cell", "goblet cell", "secretory cell"):
        return "secretory"
    if "basal cell" in s:
        return "basal"
    if "fibroblast" in s:
        return "fibroblast"
    if "endothelial cell" in s:
        return "endothelial"
    if _has(s, "pericyte", "smooth muscle", "mural cell"):
        return "mural"
    return None


SCHEMES = {"none": lambda s: s, "blood_fine": blood_fine, "oral_coarse": oral_coarse, "lung_coarse": lung_coarse}


def map_labels(labels, scheme: str) -> np.ndarray:
    """Array of coarse labels; unmapped cells get the string '' (excluded by callers)."""
    f = SCHEMES[scheme]
    cache = {}
    out = []
    for l in labels:
        if l not in cache:
            m = f(str(l))
            cache[l] = "" if m is None else m
        out.append(cache[l])
    return np.array(out, dtype=str)


def audit_table(labels, scheme: str):
    """(source_label, coarse_label or '<excluded>', n_cells) rows for the audit CSV."""
    f = SCHEMES[scheme]
    u, c = np.unique(np.asarray(labels).astype(str), return_counts=True)
    return [(str(a), f(str(a)) or "<excluded>", int(b)) for a, b in zip(u, c)]

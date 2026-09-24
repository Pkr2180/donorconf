"""Query NCBI GEO (E-utilities) for accession metadata so the manifest is verified, not assumed.

Reports title, organism, sample count and summary for each accession. The sample
count is NOT a donor count: donors must be read from the sample metadata, which
this script does not attempt to infer. Prints JSON; writes nothing else.

    python scripts/verify_geo_accessions.py GSE178360 GSE127465 --out geo_check.json
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _get(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r)


def lookup(acc: str) -> dict:
    q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]", "retmode": "json"})
    ids = _get(BASE + "esearch.fcgi?" + q)["esearchresult"]["idlist"]
    if not ids:
        return {"accession": acc, "found": False}
    q = urllib.parse.urlencode({"db": "gds", "id": ids[0], "retmode": "json"})
    rec = _get(BASE + "esummary.fcgi?" + q)["result"][ids[0]]
    return {"accession": acc, "found": True, "title": rec.get("title"), "taxon": rec.get("taxon"),
            "n_samples": rec.get("n_samples"), "gdstype": rec.get("gdstype"),
            "summary": (rec.get("summary") or "")[:400]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="+")
    ap.add_argument("--out")
    a = ap.parse_args()
    out = []
    for acc in a.accessions:
        try:
            out.append(lookup(acc))
        except Exception as e:  # network or API failure must be visible, never silently skipped
            out.append({"accession": acc, "found": None, "error": repr(e)})
        time.sleep(0.4)
    text = json.dumps(out, indent=2)
    print(text)
    if a.out:
        open(a.out, "w").write(text)


if __name__ == "__main__":
    main()

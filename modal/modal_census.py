"""Modal jobs that retrieve REAL data from CELLxGENE Census (workspace: pradeepaiperio).

Run from donorconf/:
    set MODAL_PROFILE=pradeepaiperio
    python -m modal run modal/modal_census.py::probe
    python -m modal run modal/modal_census.py::survey --tissue-like "gingiva|oral|tongue|mouth|saliv|tonsil"

Outputs land on the Modal volume `donorconf-data` and are pulled to results/real/.
Nothing here fabricates data: every row comes from the Census, and the census
version is recorded with each output.
"""
from __future__ import annotations

import json
import modal

app = modal.App("donorconf-census")
vol = modal.Volume.from_name("donorconf-data", create_if_missing=True)
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("cellxgene-census", "numpy", "pandas", "scipy", "pyarrow"))
VOL = "/data"


@app.function(image=image, timeout=1800, memory=8192)
def probe() -> dict:
    """Which Census versions and hosted embeddings exist right now."""
    import cellxgene_census as cc
    import cellxgene_census.experimental  # noqa: F401  (not imported by default)
    out = {"cellxgene_census": cc.__version__}
    for ver in ("stable", "latest"):
        try:
            d = cc.get_census_version_description(ver)
            out[f"{ver}_resolves_to"] = d.get("release_build") or d.get("release_date") or str(d)[:200]
        except Exception as e:
            out[f"{ver}_resolve_error"] = repr(e)
    try:
        out["versions"] = {k: v.get("release_date") if isinstance(v, dict) else str(v)
                           for k, v in cc.get_census_version_directory().items()}
    except Exception as e:
        out["versions_error"] = repr(e)
    for ver in ("stable", "latest"):
        try:
            emb = cc.experimental.get_all_available_embeddings(ver)
            out[f"embeddings_{ver}"] = [
                {k: e.get(k) for k in ("embedding_name", "experiment_name", "n_embeddings", "n_dimensions", "id")}
                if isinstance(e, dict) else str(e) for e in emb]
        except Exception as e:
            out[f"embeddings_{ver}_error"] = repr(e)
    return out


@app.local_entrypoint()
def probe_main():
    print(json.dumps(probe.remote(), indent=2, default=str))


ORAL_REGEX = "gingiv|oral cavity|tongue|saliv|tonsil|dental|periodont|buccal|lip$|palate|oropharyn"
ORAL_TISSUES = ["gingiva", "buccal mucosa", "periodontium", "dental pulp", "tongue", "anterior part of tongue",
                "posterior part of tongue", "mucosa of dorsum of tongue", "tonsil", "oropharynx", "oral cavity",
                "minor salivary gland", "salivary gland epithelium", "saliva", "hard palate", "soft palate",
                "mucosa of lip"]


@app.function(image=image, timeout=3600, memory=32768, volumes={VOL: vol})
def survey(census_version: str = "stable", tissue_generals: list | None = None,
           oral_regex: str = ORAL_REGEX) -> dict:
    """Real donor/dataset structure in the Census for oral-related tissues and large multi-donor tissues."""
    import re
    import cellxgene_census as cc
    import pandas as pd

    tissue_generals = tissue_generals or ["blood", "lung"]
    out = {"census_version_requested": census_version}
    with cc.open_soma(census_version=census_version) as census:
        out["census_info"] = {k: str(v)[:120] for k, v in census["census_info"]["summary"].read().concat().to_pandas()
                              .set_index("label")["value"].items()}
        counts = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
        hs = counts[(counts.organism == "homo_sapiens") & (counts.category == "tissue") & (counts.total_cell_count > 0)]
        oral_tissues = sorted(t for t in hs.label.astype(str).unique() if t in ORAL_TISSUES)
        out["oral_like_tissue_labels"] = oral_tissues
        out["oral_like_tissue_cell_counts"] = (hs[hs.label.isin(oral_tissues)]
                                                .groupby("label").total_cell_count.max().to_dict())
        obs = census["census_data"]["homo_sapiens"].obs
        cols = ["dataset_id", "donor_id", "tissue", "tissue_general", "cell_type", "disease", "assay"]
        frames = []
        if oral_tissues:
            frames.append(obs.read(value_filter=f"is_primary_data == True and tissue in {oral_tissues!r}",
                                   column_names=cols).concat().to_pandas().assign(group="oral_like"))
        for tg in tissue_generals:
            frames.append(obs.read(value_filter=f"is_primary_data == True and tissue_general == '{tg}'",
                                   column_names=cols).concat().to_pandas().assign(group=tg))
        df = pd.concat(frames, ignore_index=True)
    for c in cols + ["group"]:
        df[c] = df[c].astype(str)      # Census returns categoricals; plain strings avoid groupby blow-up
    df["donor_ns"] = df.dataset_id + "::" + df.donor_id
    tab = (df.groupby(["group", "dataset_id"], observed=True)
             .agg(n_cells=("donor_ns", "size"), n_donors=("donor_ns", "nunique"),
                  n_cell_types=("cell_type", "nunique"),
                  tissues=("tissue", lambda x: "; ".join(sorted(set(x))[:4])),
                  diseases=("disease", lambda x: "; ".join(sorted(set(x))[:4])),
                  assays=("assay", lambda x: "; ".join(sorted(set(x))[:3])))
             .reset_index().sort_values(["group", "n_donors"], ascending=[True, False]))
    tab.to_csv(f"{VOL}/survey_{census_version}.csv", index=False)
    vol.commit()
    out["n_rows"] = int(len(tab))
    out["oral_like_datasets"] = tab[tab.group == "oral_like"].to_dict("records")
    out["top_datasets_by_donors"] = {g: tab[tab.group == g].head(15).to_dict("records") for g in tissue_generals}
    return out


@app.local_entrypoint()
def survey_main(census_version: str = "stable", out: str = "../results/real/survey.json"):
    import os
    res = survey.remote(census_version)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w"), indent=2, default=str)
    print(json.dumps({k: v for k, v in res.items() if k not in ("top_datasets_by_donors",)}, indent=2, default=str)[:6000])


@app.function(image=image, timeout=1800, memory=8192)
def inspect_counts(census_version: str = "stable", regex: str = ORAL_REGEX) -> dict:
    import cellxgene_census as cc
    import pandas as pd
    with cc.open_soma(census_version=census_version) as census:
        c = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
    c["cnt"] = pd.to_numeric(c["total_cell_count"], errors="coerce")
    c = c[(c.organism.astype(str) == "homo_sapiens") & c.category.astype(str).isin(["tissue", "tissue_general"])
          & (c.cnt > 0)]
    c["label"] = c.label.astype(str)
    hits = c[c.label.str.contains(regex, case=False, regex=True)]
    return {"census_version": census_version,
            "n_human_tissue_labels": int(c.label.nunique()),
            "hit_counts": (hits.groupby(["category", "label"], observed=True).cnt.max().reset_index()
                           .sort_values("cnt", ascending=False).astype({"cnt": int}).to_dict("records"))}


@app.local_entrypoint()
def inspect_main(out: str = "../results/real/oral_tissue_counts.json"):
    import os
    res = inspect_counts.remote()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w"), indent=2)
    print("human tissue labels with cells:", res["n_human_tissue_labels"])
    for r in res["hit_counts"]:
        print(f"{r['category']:14s} {r['label'][:48]:48s} {r['cnt']:>10,d}")


@app.function(image=image, timeout=1800, memory=8192)
def embedding_api(census_version: str = "2025-11-08") -> dict:
    import cellxgene_census.experimental as cx
    out = {}
    allv = cx.get_all_available_embeddings(census_version)
    out["first_entry_raw"] = allv[0] if allv else None
    out["names"] = [(e.get("embedding_name"), e.get("experiment_name")) for e in allv]
    for name in ("tf-sapiens", "tf-exemplar-human", "scvi"):
        try:
            m = cx.get_embedding_metadata_by_name(name, "homo_sapiens", census_version=census_version)
            out[name] = {k: (v if not hasattr(v, "shape") else str(v.shape)) for k, v in dict(m).items()}
        except Exception as e:
            out[name] = "ERR " + repr(e)
    return out


@app.local_entrypoint()
def embedding_api_main():
    print(json.dumps(embedding_api.remote(), indent=2, default=str)[:5000])


EMBEDDINGS = ["tf-sapiens", "tf-exemplar-human", "scvi"]
OBS_COLS = ["soma_joinid", "dataset_id", "donor_id", "cell_type", "disease", "tissue", "assay", "sex"]


@app.function(image=image, timeout=7200, memory=32768, volumes={VOL: vol})
def pull(dataset_id: str, tag: str, max_cells_per_donor: int = 300, max_donors: int = 0,
         seed: int = 0, census_version: str = "2025-11-08", tissue_general: str = "") -> dict:
    """Subsample cells per donor from one Census dataset and attach hosted embeddings.

    Every array comes from the Census. Subsampling is random with a recorded seed.
    Cells with any missing embedding row are dropped and counted.
    """
    import numpy as np
    import pandas as pd
    import cellxgene_census as cc

    rng = np.random.default_rng(seed)
    with cc.open_soma(census_version=census_version) as census:
        obs = census["census_data"]["homo_sapiens"].obs
        vf = f"is_primary_data == True and dataset_id == '{dataset_id}'"
        if tissue_general:                     # multi-organ atlases: keep one organ only (added after a lung audit)
            vf += f" and tissue_general == '{tissue_general}'"
        df = obs.read(value_filter=vf,
                      column_names=OBS_COLS).concat().to_pandas()
        n_total, n_donors_total = len(df), df.donor_id.astype(str).nunique()
        for c in OBS_COLS[1:]:
            df[c] = df[c].astype(str)
        donors = np.array(sorted(df.donor_id.unique()))
        if max_donors and donors.size > max_donors:
            donors = rng.choice(donors, size=max_donors, replace=False)
        parts = []
        for d in donors:
            g = df[df.donor_id == d]
            parts.append(g if len(g) <= max_cells_per_donor else g.sample(max_cells_per_donor, random_state=int(rng.integers(1e9))))
        sub = pd.concat(parts).sort_values("soma_joinid").reset_index(drop=True)
        ad = cc.get_anndata(census, organism="Homo sapiens", obs_coords=sub.soma_joinid.to_numpy(),
                            var_value_filter="feature_id == 'ENSG00000000003'",
                            obs_embeddings=EMBEDDINGS, column_names={"obs": ["soma_joinid"]})
    emb = {}
    ok = np.ones(ad.n_obs, bool)
    for e in EMBEDDINGS:
        M = np.asarray(ad.obsm[e], dtype=np.float32)
        emb[e] = M
        ok &= np.isfinite(M).all(axis=1)
    order = pd.Series(np.arange(ad.n_obs), index=ad.obs["soma_joinid"].to_numpy())
    idx = order.loc[sub.soma_joinid.to_numpy()].to_numpy()
    sub = sub.assign(_ok=ok[idx])
    keep = sub._ok.to_numpy()
    out = {f"emb_{e}": emb[e][idx][keep] for e in EMBEDDINGS}
    out.update(labels=sub.cell_type.to_numpy()[keep].astype(str),
               donors=(dataset_id[:8] + "::" + sub.donor_id).to_numpy()[keep].astype(str),
               disease=sub.disease.to_numpy()[keep].astype(str),
               tissue=sub.tissue.to_numpy()[keep].astype(str),
               sex=sub.sex.to_numpy()[keep].astype(str))
    path = f"{VOL}/{tag}.npz"
    np.savez_compressed(path, **out)
    vol.commit()
    meta = {"tag": tag, "dataset_id": dataset_id, "census_version": census_version, "seed": seed,
            "n_cells_in_dataset": int(n_total), "n_donors_in_dataset": int(n_donors_total),
            "n_donors_kept": int(np.unique(out["donors"]).size), "n_cells_kept": int(keep.sum()),
            "n_cells_dropped_missing_embedding": int((~keep).sum()),
            "max_cells_per_donor": max_cells_per_donor,
            "embedding_dims": {e: int(emb[e].shape[1]) for e in EMBEDDINGS},
            "n_cell_types": int(np.unique(out["labels"]).size), "tissue_general_filter": tissue_general}
    (ad_path := f"{VOL}/{tag}.meta.json")
    open(ad_path, "w").write(json.dumps(meta, indent=2))
    vol.commit()
    return meta


@app.local_entrypoint()
def pull_main(dataset_id: str, tag: str, max_cells_per_donor: int = 300, max_donors: int = 0, seed: int = 0,
              tissue_general: str = ""):
    import os
    meta = pull.remote(dataset_id, tag, max_cells_per_donor, max_donors, seed, tissue_general=tissue_general)
    os.makedirs("../results/real/meta", exist_ok=True)
    json.dump(meta, open(f"../results/real/meta/{tag}.json", "w"), indent=2)
    print(json.dumps(meta, indent=2))

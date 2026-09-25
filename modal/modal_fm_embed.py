"""Modal GPU jobs: Geneformer and scGPT embeddings for donorconf real datasets (workspace: pradeepaiperio).

Neither model is Census-hosted (unlike tf-sapiens/tf-exemplar-human/scvi in modal_census.py), so this
reads the raw-counts bundle `{tag}_raw.npz` written by `modal_census.py::pull_raw` (same cell subsample
as the existing embeddings, guaranteed by reusing pull()'s exact rng/order logic) and runs each
foundation model's own pretrained checkpoint to produce per-cell embeddings, aligned by `cell_id`
(soma_joinid) back to the same row order as `{tag}.npz`.

Design choices made after checking each project's current docs (2026-09-24), to avoid burning GPU time
on dependency dead ends:
  - Geneformer V2: TranscriptomeTokenizer accepts .h5ad directly (file_format="h5ad") -- no loom needed.
    Requires var["ensembl_id"] and obs["n_counts"]. Checkpoint: ctheodoris/Geneformer, subfolder
    "Geneformer-V2-104M" (general-purpose, non-cancer).
  - scGPT: flash-attention is an OPTIONAL dependency (per bowang-lab docs) -- explicitly pass
    use_fast_transformer=False to run on plain PyTorch attention and skip the flash-attn build
    entirely (its pinned CUDA/version requirements are a known source of broken installs).
    Checkpoint: mirrored on Hugging Face at wanglab/scGPT-human (args.json, best_model.pt, vocab.json)
    -- avoids scripting a Google Drive download from a headless container.

Run from donorconf/:
    set MODAL_PROFILE=pradeepaiperio
    python -m modal run modal/modal_fm_embed.py::embed_geneformer_main --tag lung_C
    python -m modal run modal/modal_fm_embed.py::embed_scgpt_main --tag lung_C
"""
from __future__ import annotations

import json
import modal

app = modal.App("donorconf-fm-embed")
vol = modal.Volume.from_name("donorconf-data", create_if_missing=True)
VOL = "/data"
GPU = "A10G"

geneformer_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "git-lfs")
    # geneformer's own setup.py pins no transformers version, so an unpinned install pulls the latest
    # release -- which has dropped the top-level `SpecialTokensMixin` re-export that geneformer's
    # collator code imports. Pin to a version still compatible.
    .pip_install("numpy", "pandas", "scipy", "anndata", "scanpy",
                 "torch", "transformers<4.50", "datasets", "accelerate", "huggingface_hub")
    # `pip install git+https://...` makes pip do its own partial clone (--filter=blob:none), which
    # Hugging Face's git server rejects for this LFS-heavy repo ("expected 'packfile'" / promisor-remote
    # error). Clone it ourselves instead and `pip install .` the package code directly.
    # NOTE: do NOT set GIT_LFS_SKIP_SMUDGE here -- the package's own gene-median/token dictionary
    # pickles (loaded by TranscriptomeTokenizer) are themselves stored via git-lfs *inside* this repo,
    # not fetched separately. Skipping smudge left LFS pointer text in place of the real pickles
    # (UnpicklingError: invalid load key, 'v' -- the 'v' of "version https://git-lfs..."). A normal
    # `git lfs install` + clone (not partial) smudges them correctly; `git lfs pull` is a safety net
    # in case shallow-clone checkout missed any.
    .run_commands(
        "git lfs install --skip-repo",
        "git clone --depth 1 https://huggingface.co/ctheodoris/Geneformer /tmp/Geneformer",
        "cd /tmp/Geneformer && git lfs pull && pip install .",
    )
)

scgpt_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("numpy", "pandas", "scipy", "anndata", "scanpy",
                 # scgpt's tokenizer imports torchtext, whose compiled .so is ABI-locked to the torch
                 # version it was built against. torchtext 0.18.0 (its final release -- the project is
                 # abandoned) was built for torch 2.3.0; an unpinned "torch" install pulled 2.14.0,
                 # causing "OSError: Could not load this library: .../libtorchtext.so" at import time.
                 # Pin both explicitly to a matched pair.
                 "torch==2.3.0", "torchtext==0.18.0",
                 "huggingface_hub", "ipython")  # scgpt.utils.util imports IPython unconditionally
    .pip_install("scgpt")
)


def _load_raw_bundle(tag: str):
    import numpy as np
    import scipy.sparse as sp

    d = np.load(f"{VOL}/{tag}_raw.npz", allow_pickle=False)
    X = sp.csr_matrix((d["X_data"], d["X_indices"], d["X_indptr"]), shape=tuple(d["X_shape"]))
    return X, d["var_feature_id"], d["var_feature_name"], d["cell_id"], d["labels"], d["donors"]


@app.function(image=geneformer_image, gpu=GPU, timeout=7200, memory=32768, volumes={VOL: vol})
def embed_geneformer(tag: str, model_dir_name: str = "Geneformer-V2-104M", forward_batch_size: int = 24) -> dict:
    import os
    import shutil
    import numpy as np
    import pandas as pd
    import anndata as ad
    from huggingface_hub import snapshot_download
    from geneformer import TranscriptomeTokenizer, EmbExtractor

    X, feat_id, feat_name, cell_id, labels, donors = _load_raw_bundle(tag)
    n_cells = X.shape[0]

    work = f"/tmp/gf_{tag}"
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(f"{work}/in", exist_ok=True)
    os.makedirs(f"{work}/tok", exist_ok=True)
    os.makedirs(f"{work}/out", exist_ok=True)

    a = ad.AnnData(X=X.tocsr())
    a.var["ensembl_id"] = [str(x).split(".")[0] for x in feat_id]  # strip Ensembl version suffix if present
    a.var_names = a.var["ensembl_id"].to_numpy()
    a.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
    a.obs["cell_id"] = [str(x) for x in cell_id]
    a.obs_names = a.obs["cell_id"].to_numpy()
    a.write_h5ad(f"{work}/in/{tag}.h5ad")

    model_dir = snapshot_download("ctheodoris/Geneformer", allow_patterns=f"{model_dir_name}/*")
    model_path = f"{model_dir}/{model_dir_name}"

    tk = TranscriptomeTokenizer({"cell_id": "cell_id"}, nproc=4, model_version="V2")
    tk.tokenize_data(f"{work}/in", f"{work}/tok", tag, file_format="h5ad")

    # emb_label is required for EmbExtractor to carry the custom "cell_id" attribute (set via
    # TranscriptomeTokenizer's custom_attr_name_dict above) through into extract_embs' output columns.
    # Without it, the returned DataFrame has no cell_id column at all and the realignment-by-cell_id
    # below silently matches zero rows (caught by this script's own n_matched/n_dropped check, which
    # is exactly why that check exists rather than trusting row order blindly).
    ee = EmbExtractor(model_type="Pretrained", emb_mode="cls", max_ncells=None,
                       forward_batch_size=forward_batch_size, nproc=4, model_version="V2",
                       emb_label=["cell_id"])
    emb_df = ee.extract_embs(model_path, f"{work}/tok/{tag}.dataset", f"{work}/out", f"{tag}_emb")

    if not isinstance(emb_df, pd.DataFrame):
        emb_df = pd.read_csv(f"{work}/out/{tag}_emb.csv")
    dim_cols = [c for c in emb_df.columns if str(c).isdigit() or (isinstance(c, int))]
    if not dim_cols:
        dim_cols = [c for c in emb_df.columns if c not in ("cell_id",)]
    emb_df = emb_df.set_index("cell_id") if "cell_id" in emb_df.columns else emb_df
    emb_df.index = emb_df.index.astype(str)

    out_dim = len(dim_cols)
    result = np.full((n_cells, out_dim), np.nan, dtype=np.float32)
    order = {c: i for i, c in enumerate(a.obs["cell_id"].to_numpy())}
    n_matched = 0
    for cid, row in emb_df.iterrows():
        i = order.get(str(cid))
        if i is not None:
            result[i] = row[dim_cols].to_numpy(dtype=np.float32)
            n_matched += 1

    np.save(f"{VOL}/{tag}_geneformer.npy", result)
    vol.commit()
    meta = {"tag": tag, "model": model_dir_name, "n_cells": int(n_cells), "n_matched": int(n_matched),
            "emb_dim": int(out_dim), "n_dropped_by_tokenizer": int(n_cells - n_matched)}
    open(f"{VOL}/{tag}_geneformer.meta.json", "w").write(json.dumps(meta, indent=2))
    vol.commit()
    return meta


@app.local_entrypoint()
def embed_geneformer_main(tag: str, forward_batch_size: int = 24):
    meta = embed_geneformer.remote(tag, forward_batch_size=forward_batch_size)
    print(json.dumps(meta, indent=2))


@app.function(image=scgpt_image, gpu=GPU, timeout=7200, memory=32768, volumes={VOL: vol})
def embed_scgpt(tag: str, batch_size: int = 64) -> dict:
    import numpy as np
    import anndata as ad
    from huggingface_hub import snapshot_download
    from scgpt.tasks import embed_data

    X, feat_id, feat_name, cell_id, labels, donors = _load_raw_bundle(tag)
    n_cells = X.shape[0]

    a = ad.AnnData(X=X.tocsr())
    a.var["feature_name"] = [str(x) for x in feat_name]
    a.var_names = a.var["feature_name"].to_numpy()
    a.obs["cell_id"] = [str(x) for x in cell_id]
    a.obs_names = a.obs["cell_id"].to_numpy()
    a.var_names_make_unique()

    model_dir = snapshot_download("wanglab/scGPT-human")

    out = embed_data(a, model_dir, gene_col="feature_name", max_length=1200, batch_size=batch_size,
                      obs_to_save=["cell_id"], device="cuda", use_fast_transformer=False,
                      return_new_adata=True)
    # with return_new_adata=True, embed_data returns the embeddings as .X (not obsm["X_scGPT"] --
    # that key is only populated when return_new_adata=False, modifying the input adata in place).
    raw_emb = np.asarray(out.X, dtype=np.float32)

    # embed_data uses a SequentialSampler (no shuffling), so row order should already match `a`, but
    # realign explicitly by cell_id (same defensive pattern as embed_geneformer) rather than assume it.
    out_cell_id = out.obs["cell_id"].astype(str).to_numpy()
    expected_cell_id = a.obs["cell_id"].to_numpy()
    if out_cell_id.shape[0] != n_cells or not np.array_equal(out_cell_id, expected_cell_id):
        order = {c: i for i, c in enumerate(expected_cell_id)}
        emb = np.full((n_cells, raw_emb.shape[1]), np.nan, dtype=np.float32)
        for i, cid in enumerate(out_cell_id):
            j = order.get(cid)
            if j is not None:
                emb[j] = raw_emb[i]
    else:
        emb = raw_emb

    np.save(f"{VOL}/{tag}_scgpt.npy", emb)
    vol.commit()
    meta = {"tag": tag, "model": "wanglab/scGPT-human", "n_cells": int(n_cells), "emb_dim": int(emb.shape[1])}
    open(f"{VOL}/{tag}_scgpt.meta.json", "w").write(json.dumps(meta, indent=2))
    vol.commit()
    return meta


@app.local_entrypoint()
def embed_scgpt_main(tag: str, batch_size: int = 64):
    meta = embed_scgpt.remote(tag, batch_size=batch_size)
    print(json.dumps(meta, indent=2))

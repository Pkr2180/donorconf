#!/bin/bash
# Scale the validated lung_C Geneformer/scGPT pipeline to the remaining 10 datasets used in the
# C benchmark (blood_A..E, lung_E/F/G/H/I -- lung_C already done). Run from donorconf/ root:
#   bash scripts/run_fm_pipeline.sh
#
# Per dataset: pull_raw (writes {tag}_raw.npz to the Modal volume; {tag}.npz is also (re)written
# there but the LOCAL results/real/data/{tag}.npz already exists from the original pull() run and
# is guaranteed cell-identical -- same seed/rng path, verified exactly on lung_C -- so it is not
# re-fetched) -> embed_geneformer -> embed_scgpt -> fetch both .npy -> merge_fm_embeddings.py.
#
# Known recurring transient failure: a local DNS blip ("getaddrinfo failed") that resolves itself
# within seconds. Each step retries up to 3 times before giving up on that dataset. A dataset that
# fails does not block the others -- failures are logged to STATUS_FILE and reported at the end.
#
# RESUMABLE: safe to Ctrl-C or lose (laptop sleep/shutdown) and re-run verbatim. A tag already
# merged (results/real/data/{tag}.npz has both emb_geneformer and emb_scgpt) is skipped, so nothing
# already paid for in GPU time is repeated. The status file is appended to, not truncated.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

export MODAL_PROFILE=pradeepaiperio
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8

STATUS_FILE="results/real/meta/fm_pipeline_status.tsv"
[ -f "$STATUS_FILE" ] || echo -e "tag\tstep\tstatus\ttimestamp" > "$STATUS_FILE"

log_status() {
    echo -e "$1\t$2\t$3\t$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$STATUS_FILE"
}

already_done() {
    python3 -c "
import numpy as np, sys
try:
    d = np.load('results/real/data/$1.npz', allow_pickle=False)
    sys.exit(0 if ('emb_geneformer' in d and 'emb_scgpt' in d) else 1)
except Exception:
    sys.exit(1)
"
}

# tag -> "dataset_id max_donors tissue_general" (tissue_general may be empty string)
declare -A DS=(
    [blood_A]="3faad104-2ab8-4434-816d-474d8d2641db 80 "
    [blood_B]="c838aec3-03ef-4398-b882-0e3912abfff0 80 "
    [blood_C]="218acb0f-9f2f-4f76-b90b-15a4b7c7f629 80 "
    [blood_D]="ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3 80 "
    [blood_E]="c7775e88-49bf-4ba2-a03b-93f00447c958 80 "
    [lung_E]="d8da613f-e681-4c69-b463-e94f5e66847f 27 "
    [lung_F]="9f222629-9e39-47d0-b83f-e08d610c7479 40 lung"
    [lung_G]="f14bc322-1322-4184-8d16-409557525ea5 40 lung"
    [lung_H]="d68a8b48-abf0-4bd5-8834-155d34ec448a 21 lung"
    [lung_I]="1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1 40 lung"
)
ORDER="blood_A blood_B blood_C blood_D blood_E lung_E lung_F lung_G lung_H lung_I"

run_with_retry() {
    # $1 = log file, rest = command
    local logf="$1"; shift
    local tries=0 max_tries=3
    while [ $tries -lt $max_tries ]; do
        tries=$((tries + 1))
        if "$@" > "$logf" 2>&1; then
            echo "EXIT=0" >> "$logf"
            return 0
        fi
        echo "attempt $tries/$max_tries failed, see $logf" >&2
        sleep 8
    done
    echo "EXIT=1" >> "$logf"
    return 1
}

for tag in $ORDER; do
    if already_done "$tag"; then
        log_status "$tag" "DONE" "SKIP_ALREADY_MERGED"
        continue
    fi

    read -r dsid maxd tg <<< "${DS[$tag]}"
    echo "=== $tag (dataset_id=$dsid max_donors=$maxd tissue_general='$tg') ==="

    # 1. pull_raw (writes {tag}_raw.npz + {tag}.npz to the volume; runs from modal/ so the
    #    entrypoint's relative "../results/real/meta/{tag}.json" write lands in the right place)
    if (cd modal && run_with_retry "../results/real/meta/${tag}_pullraw.log" \
            python -m modal run modal_census.py::pull_raw_main \
            --dataset-id "$dsid" --tag "$tag" --max-cells-per-donor 150 --max-donors "$maxd" \
            --tissue-general "$tg"); then
        log_status "$tag" "pull_raw" "OK"
    else
        log_status "$tag" "pull_raw" "FAIL"
        continue
    fi

    # 2. embed_geneformer
    if (cd modal && run_with_retry "../results/real/meta/${tag}_geneformer.log" \
            python -m modal run modal_fm_embed.py::embed_geneformer_main --tag "$tag"); then
        # confirm n_matched is real, not a silent 0-match (the exact bug caught on lung_C attempt 7)
        n_matched=$(grep -o '"n_matched": [0-9]*' "results/real/meta/${tag}_geneformer.log" | tail -1 | grep -o '[0-9]*$')
        if [ -z "$n_matched" ] || [ "$n_matched" -eq 0 ]; then
            log_status "$tag" "embed_geneformer" "FAIL_ZERO_MATCH"
            continue
        fi
        log_status "$tag" "embed_geneformer" "OK n_matched=$n_matched"
    else
        log_status "$tag" "embed_geneformer" "FAIL"
        continue
    fi

    # 3. embed_scgpt
    if (cd modal && run_with_retry "../results/real/meta/${tag}_scgpt.log" \
            python -m modal run modal_fm_embed.py::embed_scgpt_main --tag "$tag"); then
        log_status "$tag" "embed_scgpt" "OK"
    else
        log_status "$tag" "embed_scgpt" "FAIL"
        continue
    fi

    # 4. fetch both .npy from the volume, place at D:\tmp so native Python's Path("/tmp/...")
    #    (used by merge_fm_embeddings.py) resolves to the same files bash downloads to /tmp
    mkdir -p /d/tmp
    ok=1
    (cd modal && run_with_retry "../results/real/meta/${tag}_fetch_gf.log" \
        python -m modal volume get donorconf-data "${tag}_geneformer.npy" "/tmp/${tag}_geneformer.npy" --force) || ok=0
    (cd modal && run_with_retry "../results/real/meta/${tag}_fetch_sc.log" \
        python -m modal volume get donorconf-data "${tag}_scgpt.npy" "/tmp/${tag}_scgpt.npy" --force) || ok=0
    if [ "$ok" -ne 1 ]; then
        log_status "$tag" "fetch" "FAIL"
        continue
    fi
    cp "/tmp/${tag}_geneformer.npy" "/tmp/${tag}_scgpt.npy" /d/tmp/
    log_status "$tag" "fetch" "OK"

    # 5. merge -- this script's own guards (row-count match, all-finite) are the last line of
    #    defense against any silent corruption in the steps above
    if python scripts/merge_fm_embeddings.py "$tag" > "results/real/meta/${tag}_merge.log" 2>&1; then
        log_status "$tag" "merge" "OK"
    else
        log_status "$tag" "merge" "FAIL"
        cat "results/real/meta/${tag}_merge.log" >&2
        continue
    fi

    log_status "$tag" "DONE" "OK"
done

echo "=== pipeline finished, status summary: ==="
cat "$STATUS_FILE"

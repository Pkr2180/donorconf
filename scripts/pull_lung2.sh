#!/bin/bash
# second lung batch, restricted to tissue_general == lung; retried on transient network errors
cd "$(dirname "$0")/../modal" || exit 1
export MODAL_PROFILE=pradeepaiperio PYTHONUTF8=1
for pair in "9f222629-9e39-47d0-b83f-e08d610c7479 lung_F" "f14bc322-1322-4184-8d16-409557525ea5 lung_G" \
            "d68a8b48-abf0-4bd5-8834-155d34ec448a lung_H" "1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1 lung_I"; do
  set -- $pair
  for i in 1 2 3 4; do
    python -m modal run modal_census.py::pull_main --dataset-id "$1" --tag "$2" --max-cells-per-donor 150 --max-donors 40 \
      --tissue-general lung > "../results/real/meta/$2.log" 2>&1
    [ -f "../results/real/meta/$2.json" ] && break
    sleep 20
  done
done
echo done > ../results/real/meta/lung2.done

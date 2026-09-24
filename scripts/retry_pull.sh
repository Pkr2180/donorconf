#!/bin/bash
# usage: retry_pull.sh <dataset_id> <tag>   (retries a Modal Census pull up to 4 times; transient DNS errors happen)
cd "$(dirname "$0")/../modal" || exit 1
export MODAL_PROFILE=pradeepaiperio PYTHONUTF8=1
for i in 1 2 3 4; do
  python -m modal run modal_census.py::pull_main --dataset-id "$1" --tag "$2" --max-cells-per-donor 150 --max-donors 40 \
    > "../results/real/meta/$2.log" 2>&1
  [ -f "../results/real/meta/$2.json" ] && break
  sleep 20
done

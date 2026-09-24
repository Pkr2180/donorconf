#!/bin/bash
# Pull five human lung datasets from CELLxGENE Census 2025-11-08 via Modal (workspace pradeepaiperio), one at a time.
cd "$(dirname "$0")/../modal" || exit 1
export MODAL_PROFILE=pradeepaiperio PYTHONUTF8=1
mkdir -p ../results/real/meta
for pair in "350237e0-9f48-4cbd-9140-3b44495549f3 lung_A" "6725ee8e-ef5b-4e68-8901-61bd14a1fe73 lung_B" \
            "eb499fd8-7000-419f-8854-926b9b61b11f lung_C" "1b350d0a-4535-4879-beb6-1142f3f94947 lung_D" \
            "d8da613f-e681-4c69-b463-e94f5e66847f lung_E"; do
  set -- $pair
  python -m modal run modal_census.py::pull_main --dataset-id "$1" --tag "$2" --max-cells-per-donor 150 --max-donors 40 \
    > "../results/real/meta/$2.log" 2>&1
done
echo done > ../results/real/meta/lung.done

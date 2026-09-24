# Data

The raw `.npz` extracts pulled from CELLxGENE Census (2025-11-08 release) are not committed to this
repository — they are tens of MB each (roughly 900 MB in total across all datasets) and are not suited to
a git repository.

To regenerate them, run `modal/modal_census.py` against Modal (see the project README for the exact
`pull_main` invocations used for each dataset tag). Alternatively, the archived data files are available
from the Zenodo record for this project (DOI to be added here once minted).

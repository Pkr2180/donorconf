"""donorconf: donor-level conformal calibration for transferred single-cell annotations.

Pre-release research code. No result produced by this package has been
validated on real data; see PROTOCOL.md for the tests that must pass first.
"""
from .adaptive import TargetAdaptiveCRC, WeightedDonorCRC
from .conformal import (CellLevelConformal, ClasswiseConformal, DonorLevelCRC, DonorWeightedQuantile,
                        conformal_quantile)
from .metrics import coverage_summary, donor_bootstrap_coverage_ci, per_donor_coverage
from .ppi import (cell_pooled_diff, confusion_corrected_diff, donor_proportions,
                  naive_plugin_diff, ppi_group_diff, ppi_mean)
from .scores import score_matrix, true_label_scores
from .shift import donor_permutation_test, donor_summaries, shift_screen

__version__ = "0.0.1"
__all__ = [
    "CellLevelConformal", "ClasswiseConformal", "DonorLevelCRC", "DonorWeightedQuantile",
    "TargetAdaptiveCRC", "WeightedDonorCRC", "conformal_quantile",
    "coverage_summary", "donor_bootstrap_coverage_ci", "per_donor_coverage",
    "cell_pooled_diff", "confusion_corrected_diff", "donor_proportions",
    "naive_plugin_diff", "ppi_group_diff", "ppi_mean",
    "score_matrix", "true_label_scores",
    "donor_permutation_test", "donor_summaries", "shift_screen",
]

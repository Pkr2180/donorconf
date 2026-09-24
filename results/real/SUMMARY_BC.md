# Framings B and C, real data (auto-generated)

**Exploratory.** Written after R1-R4; see `ANALYSIS_LOCK.md` Addendum B. Embedding and pretraining caveats from `SUMMARY_REAL.md` apply. Nothing here is a confirmatory result.

## B: target-adaptive calibration with k labeled target donors

### Blood, leave-dataset-out (A, B, C, permuted roles)
Target 0.90; scored on target donors NOT used for adaptation; pooled over embeddings.

| k | method | reps | coverage | cov_sd | set_size | p_rep_under5 | mean_omega |
|---|---|---|---|---|---|---|---|
| 3 | adaptive_crc | 30 | 0.924 | 0.029 | 1.384 | 0.011 | 0.875 |
| 3 | adaptive_plugin | 30 | 0.896 | 0.026 | 1.254 | 0.067 | nan |
| 3 | cell_pooled | 30 | 0.865 | 0.094 | 1.222 | 0.211 | nan |
| 3 | donor_crc | 30 | 0.914 | 0.055 | 1.370 | 0.067 | nan |
| 3 | donor_weighted | 30 | 0.869 | 0.093 | 1.230 | 0.200 | nan |
| 3 | target_only_k | 30 | 0.887 | 0.043 | 1.241 | 0.144 | nan |
| 6 | adaptive_crc | 30 | 0.918 | 0.021 | 1.345 | 0.000 | 0.794 |
| 6 | adaptive_plugin | 30 | 0.893 | 0.019 | 1.250 | 0.000 | nan |
| 6 | cell_pooled | 30 | 0.864 | 0.093 | 1.226 | 0.211 | nan |
| 6 | donor_crc | 30 | 0.913 | 0.055 | 1.374 | 0.067 | nan |
| 6 | donor_weighted | 30 | 0.868 | 0.091 | 1.235 | 0.200 | nan |
| 6 | target_only_k | 30 | 0.887 | 0.028 | 1.239 | 0.122 | nan |

### Oral: gingiva -> oropharynx -> oral SCC
Target 0.90; scored on target donors NOT used for adaptation; pooled over embeddings.

| k | method | reps | coverage | cov_sd | set_size | p_rep_under5 | mean_omega |
|---|---|---|---|---|---|---|---|
| 3 | adaptive_crc | 30 | 0.925 | 0.031 | 1.028 | 0.022 | 0.672 |
| 3 | adaptive_plugin | 30 | 0.896 | 0.037 | 0.948 | 0.122 | nan |
| 3 | cell_pooled | 30 | 0.865 | 0.034 | 0.900 | 0.289 | nan |
| 3 | donor_crc | 30 | 0.974 | 0.008 | 1.636 | 0.000 | nan |
| 3 | donor_weighted | 30 | 0.863 | 0.034 | 0.898 | 0.311 | nan |
| 3 | target_only_k | 30 | 0.901 | 0.043 | 0.963 | 0.111 | nan |
| 6 | adaptive_crc | 30 | 0.915 | 0.033 | 0.987 | 0.056 | 0.628 |
| 6 | adaptive_plugin | 30 | 0.896 | 0.038 | 0.948 | 0.111 | nan |
| 6 | cell_pooled | 30 | 0.864 | 0.041 | 0.899 | 0.289 | nan |
| 6 | donor_crc | 30 | 0.974 | 0.009 | 1.636 | 0.000 | nan |
| 6 | donor_weighted | 30 | 0.862 | 0.041 | 0.897 | 0.311 | nan |
| 6 | target_only_k | 30 | 0.900 | 0.039 | 0.954 | 0.089 | nan |

### Oral: gingiva -> multi-site oral -> oral SCC
Target 0.90; scored on target donors NOT used for adaptation; pooled over embeddings.

| k | method | reps | coverage | cov_sd | set_size | p_rep_under5 | mean_omega |
|---|---|---|---|---|---|---|---|
| 3 | adaptive_crc | 30 | 0.980 | 0.012 | 1.118 | 0.000 | 0.250 |
| 3 | adaptive_plugin | 30 | 0.933 | 0.025 | 0.949 | 0.000 | nan |
| 3 | cell_pooled | 30 | 0.976 | 0.015 | 1.095 | 0.000 | nan |
| 3 | donor_crc | 30 | 0.997 | 0.005 | 1.947 | 0.000 | nan |
| 3 | donor_weighted | 30 | 0.976 | 0.014 | 1.081 | 0.000 | nan |
| 3 | target_only_k | 30 | 0.890 | 0.046 | 0.900 | 0.178 | nan |
| 6 | adaptive_crc | 30 | 0.979 | 0.012 | 1.102 | 0.000 | 0.250 |
| 6 | adaptive_plugin | 30 | 0.926 | 0.023 | 0.939 | 0.011 | nan |
| 6 | cell_pooled | 30 | 0.977 | 0.014 | 1.095 | 0.000 | nan |
| 6 | donor_crc | 30 | 0.997 | 0.004 | 1.943 | 0.000 | nan |
| 6 | donor_weighted | 30 | 0.977 | 0.013 | 1.080 | 0.000 | nan |
| 6 | target_only_k | 30 | 0.901 | 0.030 | 0.910 | 0.056 | nan |

## C: blood benchmark (reported vs realised coverage under dataset shift)
Ordered dataset triples: 24; rows: 576. Target coverage 0.90. Exploratory.

### Coverage and set size by method (mean over triples, replications, embeddings)
| method | coverage | cov_sd | set_size | p_rep_under5 | frac_poor_donors |
|---|---|---|---|---|---|
| cell_pooled | 0.876 | 0.088 | 1.006 | 0.188 | 0.168 |
| classwise | 0.832 | 0.123 | 1.016 | 0.299 | 0.258 |
| donor_crc | 0.915 | 0.058 | 1.117 | 0.097 | 0.090 |
| donor_weighted | 0.877 | 0.088 | 1.008 | 0.174 | 0.163 |

### By embedding
| emb | method | coverage | set_size | p_rep_under5 |
|---|---|---|---|---|
| TF-Exemplar | cell_pooled | 0.886 | 1.099 | 0.188 |
| TF-Exemplar | classwise | 0.808 | 1.076 | 0.438 |
| TF-Exemplar | donor_crc | 0.918 | 1.244 | 0.083 |
| TF-Exemplar | donor_weighted | 0.886 | 1.102 | 0.188 |
| TF-Sapiens | cell_pooled | 0.849 | 0.958 | 0.333 |
| TF-Sapiens | classwise | 0.812 | 0.987 | 0.354 |
| TF-Sapiens | donor_crc | 0.901 | 1.080 | 0.208 |
| TF-Sapiens | donor_weighted | 0.851 | 0.961 | 0.292 |
| scVI | cell_pooled | 0.892 | 0.960 | 0.042 |
| scVI | classwise | 0.876 | 0.986 | 0.104 |
| scVI | donor_crc | 0.926 | 1.026 | 0.000 |
| scVI | donor_weighted | 0.893 | 0.962 | 0.042 |

### Reported minus realised coverage (cell_pooled)
mean 0.027, sd 0.087; share with |gap| > 0.03: 0.479; share over-reported (reported > realised + 0.03): 0.292; share under-reported (reported < realised - 0.03): 0.188

### Can a label-free signal predict shortfall? (kill rule: |Spearman| >= 0.4)
Correlations are across replications that reuse the same five datasets, so p-values are optimistic.
| label_free_signal | spearman_vs_shortfall | p_value |
|---|---|---|
| conf_shift | -0.566 | 0.000 |
| entropy_shift | 0.539 | 0.000 |
| screen_p | -0.104 | 0.213 |
| head_accuracy_calibration | 0.633 | 0.000 |
| screen flagged: AUROC for >5pt failure | 0.577 | nan |

Same correlations with replications averaged first (the fairer unit; the same five datasets still recur):
| unit | label_free_signal | spearman_vs_shortfall | p_value | n |
|---|---|---|---|---|
| per embedding x triple (n=72) | conf_shift | -0.545 | 0.000 | 72 |
| per embedding x triple (n=72) | entropy_shift | 0.521 | 0.000 | 72 |
| per embedding x triple (n=72) | screen_p | -0.099 | 0.406 | 72 |
| per embedding x triple (n=72) | head_accuracy_calibration | 0.701 | 0.000 | 72 |
| per triple, embeddings averaged (n=24) | conf_shift | -0.732 | 0.000 | 24 |
| per triple, embeddings averaged (n=24) | entropy_shift | 0.730 | 0.000 | 24 |
| per triple, embeddings averaged (n=24) | screen_p | -0.297 | 0.159 | 24 |
| per triple, embeddings averaged (n=24) | head_accuracy_calibration | 0.614 | 0.001 | 24 |

Figures: `figures/fig1_reported_vs_realised.png`, `fig2_shortfall_vs_confidence_shift.png`, `fig3_coverage_vs_size.png`.

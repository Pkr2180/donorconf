# Published-package baselines and second tissue (auto-generated)

**Exploratory** (`ANALYSIS_LOCK.md` Addendum C). tcp_* rows are the torchCP 1.0.2 predictors used by the conformalized single-cell annotator (Bioinformatics 2025, btaf521) on this study's head and splits; the annotator's OOD detector and neural classifier are not included. Embeddings are Census-hosted and probably saw these cells.

## Blood (5 datasets, 24 triples) with published-package baselines
Rows 1152; triples 24; target coverage 0.90.

### Methods (mean over triples, replications, embeddings)
| method | reps | coverage | cov_sd | set_size | p_fail5 | frac_poor_donors |
|---|---|---|---|---|---|---|
| cell_pooled | 144 | 0.876 | 0.088 | 1.006 | 0.188 | 0.168 |
| classwise | 144 | 0.832 | 0.123 | 1.016 | 0.299 | 0.258 |
| donor_crc | 144 | 0.915 | 0.058 | 1.117 | 0.097 | 0.090 |
| donor_weighted | 144 | 0.877 | 0.088 | 1.008 | 0.174 | 0.163 |
| tcp_classwise_thr | 144 | 0.835 | 0.119 | 1.027 | 0.285 | 0.253 |
| tcp_cluster_thr | 144 | 0.843 | 0.110 | 1.139 | 0.312 | 0.238 |
| tcp_standard_aps | 144 | 0.893 | 0.044 | 1.130 | 0.139 | 0.079 |
| tcp_standard_thr | 144 | 0.876 | 0.088 | 1.006 | 0.188 | 0.168 |

### By embedding (coverage / set size)
| emb | method | coverage | set_size | p_fail5 |
|---|---|---|---|---|
| TF-Exemplar | cell_pooled | 0.886 | 1.099 | 0.188 |
| TF-Exemplar | classwise | 0.808 | 1.076 | 0.438 |
| TF-Exemplar | donor_crc | 0.918 | 1.244 | 0.083 |
| TF-Exemplar | donor_weighted | 0.886 | 1.102 | 0.188 |
| TF-Exemplar | tcp_classwise_thr | 0.812 | 1.096 | 0.417 |
| TF-Exemplar | tcp_cluster_thr | 0.819 | 1.201 | 0.438 |
| TF-Exemplar | tcp_standard_aps | 0.890 | 1.198 | 0.208 |
| TF-Exemplar | tcp_standard_thr | 0.886 | 1.099 | 0.188 |
| TF-Sapiens | cell_pooled | 0.849 | 0.958 | 0.333 |
| TF-Sapiens | classwise | 0.812 | 0.987 | 0.354 |
| TF-Sapiens | donor_crc | 0.901 | 1.080 | 0.208 |
| TF-Sapiens | donor_weighted | 0.851 | 0.961 | 0.292 |
| TF-Sapiens | tcp_classwise_thr | 0.818 | 0.999 | 0.333 |
| TF-Sapiens | tcp_cluster_thr | 0.827 | 1.109 | 0.375 |
| TF-Sapiens | tcp_standard_aps | 0.890 | 1.117 | 0.208 |
| TF-Sapiens | tcp_standard_thr | 0.849 | 0.958 | 0.333 |
| scVI | cell_pooled | 0.892 | 0.960 | 0.042 |
| scVI | classwise | 0.876 | 0.986 | 0.104 |
| scVI | donor_crc | 0.926 | 1.026 | 0.000 |
| scVI | donor_weighted | 0.893 | 0.962 | 0.042 |
| scVI | tcp_classwise_thr | 0.876 | 0.986 | 0.104 |
| scVI | tcp_cluster_thr | 0.883 | 1.106 | 0.125 |
| scVI | tcp_standard_aps | 0.898 | 1.076 | 0.000 |
| scVI | tcp_standard_thr | 0.892 | 0.960 | 0.042 |

### Reported vs realised (cell_pooled)
reported 0.903 vs realised 0.876; |gap| > 0.03 in 47.9% of replications (over-reported 29.2%, under-reported 18.8%)

### Label-free signals vs shortfall (per triple; kill rule |Spearman| >= 0.4)
| signal | spearman_per_triple | p_value | n |
|---|---|---|---|
| conf_shift | -0.732 | 0.000 | 24 |
| entropy_shift | 0.730 | 0.000 | 24 |
| screen_p | -0.297 | 0.159 | 24 |
| head_accuracy_calibration | 0.614 | 0.001 | 24 |

Cross-check: torchCP standard-THR vs this package's pooled-cell coverage, mean absolute difference 0.0000 over 144 runs.

## Lung (6 datasets: C, E, F, G, H, I), second tissue
Rows 1152; triples 24; target coverage 0.90.

### Methods (mean over triples, replications, embeddings)
| method | reps | coverage | cov_sd | set_size | p_fail5 | frac_poor_donors |
|---|---|---|---|---|---|---|
| cell_pooled | 144 | 0.822 | 0.172 | 1.079 | 0.382 | 0.286 |
| classwise | 144 | 0.784 | 0.175 | 1.012 | 0.542 | 0.374 |
| donor_crc | 144 | 0.898 | 0.116 | 1.266 | 0.215 | 0.138 |
| donor_weighted | 144 | 0.837 | 0.163 | 1.115 | 0.368 | 0.258 |
| tcp_classwise_thr | 144 | 0.788 | 0.167 | 1.044 | 0.542 | 0.370 |
| tcp_cluster_thr | 144 | 0.898 | 0.166 | 2.261 | 0.208 | 0.143 |
| tcp_standard_aps | 144 | 0.890 | 0.046 | 1.237 | 0.125 | 0.100 |
| tcp_standard_thr | 144 | 0.823 | 0.169 | 1.080 | 0.382 | 0.284 |

### By embedding (coverage / set size)
| emb | method | coverage | set_size | p_fail5 |
|---|---|---|---|---|
| TF-Exemplar | cell_pooled | 0.802 | 1.153 | 0.458 |
| TF-Exemplar | classwise | 0.769 | 1.048 | 0.542 |
| TF-Exemplar | donor_crc | 0.883 | 1.354 | 0.271 |
| TF-Exemplar | donor_weighted | 0.819 | 1.193 | 0.438 |
| TF-Exemplar | tcp_classwise_thr | 0.781 | 1.123 | 0.542 |
| TF-Exemplar | tcp_cluster_thr | 0.878 | 2.303 | 0.229 |
| TF-Exemplar | tcp_standard_aps | 0.891 | 1.324 | 0.083 |
| TF-Exemplar | tcp_standard_thr | 0.806 | 1.158 | 0.458 |
| TF-Sapiens | cell_pooled | 0.800 | 1.120 | 0.333 |
| TF-Sapiens | classwise | 0.766 | 1.047 | 0.521 |
| TF-Sapiens | donor_crc | 0.892 | 1.349 | 0.208 |
| TF-Sapiens | donor_weighted | 0.813 | 1.160 | 0.354 |
| TF-Sapiens | tcp_classwise_thr | 0.767 | 1.069 | 0.521 |
| TF-Sapiens | tcp_cluster_thr | 0.895 | 2.274 | 0.167 |
| TF-Sapiens | tcp_standard_aps | 0.889 | 1.297 | 0.146 |
| TF-Sapiens | tcp_standard_thr | 0.801 | 1.120 | 0.333 |
| scVI | cell_pooled | 0.862 | 0.963 | 0.354 |
| scVI | classwise | 0.816 | 0.941 | 0.562 |
| scVI | donor_crc | 0.918 | 1.096 | 0.167 |
| scVI | donor_weighted | 0.879 | 0.992 | 0.312 |
| scVI | tcp_classwise_thr | 0.816 | 0.941 | 0.562 |
| scVI | tcp_cluster_thr | 0.921 | 2.206 | 0.229 |
| scVI | tcp_standard_aps | 0.892 | 1.090 | 0.146 |
| scVI | tcp_standard_thr | 0.862 | 0.963 | 0.354 |

### Reported vs realised (cell_pooled)
reported 0.903 vs realised 0.822; |gap| > 0.03 in 71.5% of replications (over-reported 47.2%, under-reported 24.3%)

### Label-free signals vs shortfall (per triple; kill rule |Spearman| >= 0.4)
| signal | spearman_per_triple | p_value | n |
|---|---|---|---|
| conf_shift | -0.776 | 0.000 | 24 |
| entropy_shift | 0.788 | 0.000 | 24 |
| screen_p | 0.106 | 0.622 | 24 |
| head_accuracy_calibration | 0.354 | 0.090 | 24 |

Cross-check: torchCP standard-THR vs this package's pooled-cell coverage, mean absolute difference 0.0013 over 144 runs.

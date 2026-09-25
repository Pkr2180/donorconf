# Published-package baselines and second tissue (auto-generated)

**Exploratory** (`ANALYSIS_LOCK.md` Addendum C). tcp_* rows are the torchCP 1.0.2 predictors used by the conformalized single-cell annotator (Bioinformatics 2025, btaf521) on this study's head and splits; the annotator's OOD detector and neural classifier are not included. Embeddings are Census-hosted and probably saw these cells.

## Blood (5 datasets, 24 triples) with published-package baselines
Rows 1920; triples 24; target coverage 0.90.

### Methods (mean over triples, replications, embeddings)
| method | reps | coverage | cov_sd | set_size | p_fail5 | frac_poor_donors |
|---|---|---|---|---|---|---|
| cell_pooled | 240 | 0.870 | 0.094 | 1.097 | 0.246 | 0.198 |
| classwise | 240 | 0.821 | 0.120 | 1.050 | 0.388 | 0.303 |
| donor_crc | 240 | 0.906 | 0.075 | 1.217 | 0.158 | 0.121 |
| donor_weighted | 240 | 0.871 | 0.094 | 1.098 | 0.237 | 0.194 |
| tcp_classwise_thr | 240 | 0.823 | 0.118 | 1.063 | 0.375 | 0.299 |
| tcp_cluster_thr | 240 | 0.833 | 0.111 | 1.182 | 0.396 | 0.282 |
| tcp_standard_aps | 240 | 0.887 | 0.064 | 1.229 | 0.183 | 0.114 |
| tcp_standard_thr | 240 | 0.870 | 0.094 | 1.115 | 0.246 | 0.198 |

### By embedding (coverage / set size)
| emb | method | coverage | set_size | p_fail5 |
|---|---|---|---|---|
| Geneformer | cell_pooled | 0.840 | 1.315 | 0.375 |
| Geneformer | classwise | 0.793 | 1.062 | 0.521 |
| Geneformer | donor_crc | 0.870 | 1.444 | 0.354 |
| Geneformer | donor_weighted | 0.840 | 1.315 | 0.375 |
| Geneformer | tcp_classwise_thr | 0.794 | 1.091 | 0.521 |
| Geneformer | tcp_cluster_thr | 0.804 | 1.231 | 0.500 |
| Geneformer | tcp_standard_aps | 0.859 | 1.504 | 0.333 |
| Geneformer | tcp_standard_thr | 0.841 | 1.407 | 0.375 |
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
| scGPT | cell_pooled | 0.884 | 1.152 | 0.292 |
| scGPT | classwise | 0.815 | 1.137 | 0.521 |
| scGPT | donor_crc | 0.914 | 1.293 | 0.146 |
| scGPT | donor_weighted | 0.884 | 1.153 | 0.292 |
| scGPT | tcp_classwise_thr | 0.817 | 1.142 | 0.500 |
| scGPT | tcp_cluster_thr | 0.830 | 1.261 | 0.542 |
| scGPT | tcp_standard_aps | 0.896 | 1.253 | 0.167 |
| scGPT | tcp_standard_thr | 0.884 | 1.152 | 0.292 |
| scVI | cell_pooled | 0.892 | 0.960 | 0.042 |
| scVI | classwise | 0.876 | 0.986 | 0.104 |
| scVI | donor_crc | 0.926 | 1.026 | 0.000 |
| scVI | donor_weighted | 0.893 | 0.962 | 0.042 |
| scVI | tcp_classwise_thr | 0.876 | 0.986 | 0.104 |
| scVI | tcp_cluster_thr | 0.883 | 1.106 | 0.125 |
| scVI | tcp_standard_aps | 0.898 | 1.076 | 0.000 |
| scVI | tcp_standard_thr | 0.892 | 0.960 | 0.042 |

### Reported vs realised (cell_pooled)
reported 0.902 vs realised 0.870; |gap| > 0.03 in 56.2% of replications (over-reported 34.2%, under-reported 22.1%)

### Label-free signals vs shortfall (per triple; kill rule |Spearman| >= 0.4)
| signal | spearman_per_triple | p_value | n |
|---|---|---|---|
| conf_shift | -0.446 | 0.029 | 24 |
| entropy_shift | 0.451 | 0.027 | 24 |
| screen_p | -0.199 | 0.351 | 24 |
| head_accuracy_calibration | 0.522 | 0.009 | 24 |

Cross-check: torchCP standard-THR vs this package's pooled-cell coverage, mean absolute difference 0.0001 over 240 runs.

## Lung (6 datasets: C, E, F, G, H, I), second tissue
Rows 1920; triples 24; target coverage 0.90.

### Methods (mean over triples, replications, embeddings)
| method | reps | coverage | cov_sd | set_size | p_fail5 | frac_poor_donors |
|---|---|---|---|---|---|---|
| cell_pooled | 240 | 0.829 | 0.162 | 1.153 | 0.379 | 0.281 |
| classwise | 240 | 0.788 | 0.172 | 1.044 | 0.525 | 0.370 |
| donor_crc | 240 | 0.899 | 0.111 | 1.340 | 0.200 | 0.144 |
| donor_weighted | 240 | 0.842 | 0.153 | 1.184 | 0.354 | 0.257 |
| tcp_classwise_thr | 240 | 0.792 | 0.167 | 1.072 | 0.525 | 0.365 |
| tcp_cluster_thr | 240 | 0.906 | 0.147 | 2.279 | 0.200 | 0.136 |
| tcp_standard_aps | 240 | 0.888 | 0.059 | 1.295 | 0.133 | 0.107 |
| tcp_standard_thr | 240 | 0.832 | 0.161 | 1.162 | 0.375 | 0.278 |

### By embedding (coverage / set size)
| emb | method | coverage | set_size | p_fail5 |
|---|---|---|---|---|
| Geneformer | cell_pooled | 0.823 | 1.481 | 0.438 |
| Geneformer | classwise | 0.770 | 1.085 | 0.542 |
| Geneformer | donor_crc | 0.885 | 1.693 | 0.250 |
| Geneformer | donor_weighted | 0.839 | 1.529 | 0.375 |
| Geneformer | tcp_classwise_thr | 0.775 | 1.122 | 0.542 |
| Geneformer | tcp_cluster_thr | 0.904 | 2.325 | 0.271 |
| Geneformer | tcp_standard_aps | 0.876 | 1.618 | 0.208 |
| Geneformer | tcp_standard_thr | 0.830 | 1.520 | 0.417 |
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
| scGPT | cell_pooled | 0.859 | 1.048 | 0.312 |
| scGPT | classwise | 0.820 | 1.097 | 0.458 |
| scGPT | donor_crc | 0.918 | 1.209 | 0.104 |
| scGPT | donor_weighted | 0.863 | 1.047 | 0.292 |
| scGPT | tcp_classwise_thr | 0.820 | 1.104 | 0.458 |
| scGPT | tcp_cluster_thr | 0.933 | 2.284 | 0.104 |
| scGPT | tcp_standard_aps | 0.895 | 1.147 | 0.083 |
| scGPT | tcp_standard_thr | 0.859 | 1.048 | 0.312 |
| scVI | cell_pooled | 0.862 | 0.963 | 0.354 |
| scVI | classwise | 0.816 | 0.941 | 0.562 |
| scVI | donor_crc | 0.918 | 1.096 | 0.167 |
| scVI | donor_weighted | 0.879 | 0.992 | 0.312 |
| scVI | tcp_classwise_thr | 0.816 | 0.941 | 0.562 |
| scVI | tcp_cluster_thr | 0.921 | 2.206 | 0.229 |
| scVI | tcp_standard_aps | 0.892 | 1.090 | 0.146 |
| scVI | tcp_standard_thr | 0.862 | 0.963 | 0.354 |

### Reported vs realised (cell_pooled)
reported 0.903 vs realised 0.829; |gap| > 0.03 in 72.5% of replications (over-reported 45.0%, under-reported 27.5%)

### Label-free signals vs shortfall (per triple; kill rule |Spearman| >= 0.4)
| signal | spearman_per_triple | p_value | n |
|---|---|---|---|
| conf_shift | -0.887 | 0.000 | 24 |
| entropy_shift | 0.891 | 0.000 | 24 |
| screen_p | 0.022 | 0.918 | 24 |
| head_accuracy_calibration | 0.321 | 0.126 | 24 |

Cross-check: torchCP standard-THR vs this package's pooled-cell coverage, mean absolute difference 0.0021 over 240 runs.

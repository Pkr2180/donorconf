# Real-data results, CELLxGENE Census 2025-11-08 (auto-generated)

**Exploratory.** Analysis parameters were fixed in `ANALYSIS_LOCK.md` before results, but the protocol is not externally registered. Embeddings are hosted by the Census; TranscriptFormer and scVI were trained on Census cells that probably include these donors, so evaluation is not independent of pretraining. Labels are the Census `cell_type` (curated mapping of the authors' labels), not model output.

## R1 gingiva: random donor splits (12 ref / 12 cal / 10 target donors)
_data: {"n_classes": 17, "cells_total": 10200, "cells_dropped_outside_vocabulary": 615, "n_donors": 34}_
_classes: ['B cell', 'T cell', 'capillary endothelial cell', 'conventional dendritic cell', 'endothelial cell of artery', 'endothelial cell of lymphatic vessel', 'endothelial cell of venule', 'fibroblast', 'helper T cell', 'keratinocyte', 'macrophage', 'mast cell', 'natural killer cell', 'neutrophil', 'pericyte', 'plasma cell', 'vascular associated smooth muscle cell']_

Target coverage 0.90. coverage = donor-mean coverage on TARGET donors. p_rep_under5 = share of replications more than 5 points below target.

| embedding | method | reps | coverage | cov_sd | set_size | p_rep_under5 | frac_poor_donors | reported_random_split |
|---|---|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 100 | 0.897 | 0.018 | 1.037 | 0.010 | 0.034 | 0.901 |
| emb_scvi | donor_crc | 100 | 0.973 | 0.008 | 1.463 | 0.000 | 0.000 | nan |
| emb_scvi | donor_weighted | 100 | 0.897 | 0.018 | 1.038 | 0.010 | 0.034 | nan |
| emb_tf-exemplar-human | cell_pooled | 100 | 0.898 | 0.018 | 1.038 | 0.000 | 0.027 | 0.901 |
| emb_tf-exemplar-human | donor_crc | 100 | 0.974 | 0.009 | 1.931 | 0.000 | 0.000 | nan |
| emb_tf-exemplar-human | donor_weighted | 100 | 0.898 | 0.018 | 1.039 | 0.000 | 0.027 | nan |
| emb_tf-sapiens | cell_pooled | 100 | 0.897 | 0.019 | 1.034 | 0.010 | 0.024 | 0.900 |
| emb_tf-sapiens | donor_crc | 100 | 0.973 | 0.008 | 1.852 | 0.000 | 0.000 | nan |
| emb_tf-sapiens | donor_weighted | 100 | 0.897 | 0.019 | 1.035 | 0.010 | 0.024 | nan |

Reported (within-calibration random split) minus realised (target donors) coverage for cell_pooled, by embedding:
| embedding | reported_minus_realised |
|---|---|
| emb_scvi | 0.004 |
| emb_tf-exemplar-human | 0.003 |
| emb_tf-sapiens | 0.003 |

Shift screen vs coverage failure (fail = donor-mean coverage >5 points below target). A useful screen has p_fail_given_flagged well above p_fail_given_not_flagged:
| embedding | method_failing | screen_flag_rate | p_fail | p_fail_given_flagged | p_fail_given_not_flagged | head_accuracy_target |
|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 0.020 | 0.010 | 0.500 | 0.000 | 0.882 |
| emb_tf-exemplar-human | cell_pooled | 0.070 | 0.000 | 0.000 | 0.000 | 0.883 |
| emb_tf-sapiens | cell_pooled | 0.050 | 0.010 | 0.200 | 0.000 | 0.884 |
| emb_scvi | donor_crc | 0.020 | 0.000 | 0.000 | 0.000 | 0.882 |
| emb_tf-exemplar-human | donor_crc | 0.070 | 0.000 | 0.000 | 0.000 | 0.883 |
| emb_tf-sapiens | donor_crc | 0.050 | 0.000 | 0.000 | 0.000 | 0.884 |

## R2 blood: leave-dataset-out (30 donors per role, 7 coarse types)
_data: {"n_classes": 6, "sets": ["blood_A", "blood_B", "blood_C"], "per_set": {"blood_A": {"cells": 12000, "kept": 11579, "donors": 80}, "blood_B": {"cells": 12000, "kept": 10049, "donors": 80}, "blood_C": {"cells": 12000, "kept": 11663, "donors": 80}}}_
_classes: ['B', 'CD4_T', 'CD8_T', 'NK', 'mono_classical', 'mono_nonclassical']_

Target coverage 0.90. coverage = donor-mean coverage on TARGET donors. p_rep_under5 = share of replications more than 5 points below target.

| embedding | method | reps | coverage | cov_sd | set_size | p_rep_under5 | frac_poor_donors | reported_random_split |
|---|---|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 60 | 0.893 | 0.031 | 0.979 | 0.083 | 0.065 | 0.897 |
| emb_scvi | donor_crc | 60 | 0.929 | 0.020 | 1.059 | 0.000 | 0.020 | nan |
| emb_scvi | donor_weighted | 60 | 0.897 | 0.030 | 0.985 | 0.050 | 0.058 | nan |
| emb_tf-exemplar-human | cell_pooled | 60 | 0.870 | 0.089 | 1.439 | 0.200 | 0.191 | 0.898 |
| emb_tf-exemplar-human | donor_crc | 60 | 0.920 | 0.053 | 1.633 | 0.050 | 0.073 | nan |
| emb_tf-exemplar-human | donor_weighted | 60 | 0.876 | 0.087 | 1.452 | 0.183 | 0.178 | nan |
| emb_tf-sapiens | cell_pooled | 60 | 0.845 | 0.124 | 1.318 | 0.283 | 0.262 | 0.899 |
| emb_tf-sapiens | donor_crc | 60 | 0.907 | 0.067 | 1.493 | 0.117 | 0.114 | nan |
| emb_tf-sapiens | donor_weighted | 60 | 0.852 | 0.120 | 1.327 | 0.250 | 0.241 | nan |

Reported (within-calibration random split) minus realised (target donors) coverage for cell_pooled, by embedding:
| embedding | reported_minus_realised |
|---|---|
| emb_scvi | 0.004 |
| emb_tf-exemplar-human | 0.028 |
| emb_tf-sapiens | 0.055 |

Shift screen vs coverage failure (fail = donor-mean coverage >5 points below target). A useful screen has p_fail_given_flagged well above p_fail_given_not_flagged:
| embedding | method_failing | screen_flag_rate | p_fail | p_fail_given_flagged | p_fail_given_not_flagged | head_accuracy_target |
|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 0.917 | 0.083 | 0.091 | 0.000 | 0.908 |
| emb_tf-exemplar-human | cell_pooled | 0.850 | 0.200 | 0.216 | 0.111 | 0.829 |
| emb_tf-sapiens | cell_pooled | 0.817 | 0.283 | 0.306 | 0.182 | 0.843 |
| emb_scvi | donor_crc | 0.917 | 0.000 | 0.000 | 0.000 | 0.908 |
| emb_tf-exemplar-human | donor_crc | 0.850 | 0.050 | 0.059 | 0.000 | 0.829 |
| emb_tf-sapiens | donor_crc | 0.817 | 0.117 | 0.143 | 0.000 | 0.843 |

## R3 oral: gingiva -> oropharynx -> oral SCC
_data: {"n_classes": 6, "sets": ["gingiva", "oropharynx", "oral_scc"], "per_set": {"gingiva": {"cells": 10200, "kept": 7487, "donors": 34}, "oropharynx": {"cells": 3000, "kept": 2629, "donors": 10}, "oral_scc": {"cells": 3460, "kept": 3374, "donors": 12}}}_
_classes: ['B', 'T', 'endothelial', 'epithelial', 'fibroblast', 'myeloid']_

Target coverage 0.90. coverage = donor-mean coverage on TARGET donors. p_rep_under5 = share of replications more than 5 points below target.

| embedding | method | reps | coverage | cov_sd | set_size | p_rep_under5 | frac_poor_donors | reported_random_split |
|---|---|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 60 | 0.833 | 0.035 | 0.868 | 0.717 | 0.300 | 0.900 |
| emb_scvi | donor_crc | 60 | 0.969 | 0.007 | 1.423 | 0.000 | 0.000 | nan |
| emb_scvi | donor_weighted | 60 | 0.832 | 0.035 | 0.866 | 0.717 | 0.300 | nan |
| emb_tf-exemplar-human | cell_pooled | 60 | 0.882 | 0.022 | 0.918 | 0.083 | 0.098 | 0.902 |
| emb_tf-exemplar-human | donor_crc | 60 | 0.978 | 0.006 | 1.836 | 0.000 | 0.000 | nan |
| emb_tf-exemplar-human | donor_weighted | 60 | 0.880 | 0.022 | 0.916 | 0.083 | 0.100 | nan |
| emb_tf-sapiens | cell_pooled | 60 | 0.870 | 0.026 | 0.908 | 0.250 | 0.153 | 0.901 |
| emb_tf-sapiens | donor_crc | 60 | 0.973 | 0.007 | 1.674 | 0.000 | 0.000 | nan |
| emb_tf-sapiens | donor_weighted | 60 | 0.869 | 0.026 | 0.906 | 0.267 | 0.160 | nan |

Reported (within-calibration random split) minus realised (target donors) coverage for cell_pooled, by embedding:
| embedding | reported_minus_realised |
|---|---|
| emb_scvi | 0.067 |
| emb_tf-exemplar-human | 0.020 |
| emb_tf-sapiens | 0.031 |

Shift screen vs coverage failure (fail = donor-mean coverage >5 points below target). A useful screen has p_fail_given_flagged well above p_fail_given_not_flagged:
| embedding | method_failing | screen_flag_rate | p_fail | p_fail_given_flagged | p_fail_given_not_flagged | head_accuracy_target |
|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 0.417 | 0.717 | 1.000 | 0.514 | 0.914 |
| emb_tf-exemplar-human | cell_pooled | 0.333 | 0.083 | 0.150 | 0.050 | 0.932 |
| emb_tf-sapiens | cell_pooled | 0.483 | 0.250 | 0.345 | 0.161 | 0.928 |
| emb_scvi | donor_crc | 0.417 | 0.000 | 0.000 | 0.000 | 0.914 |
| emb_tf-exemplar-human | donor_crc | 0.333 | 0.000 | 0.000 | 0.000 | 0.932 |
| emb_tf-sapiens | donor_crc | 0.483 | 0.000 | 0.000 | 0.000 | 0.928 |

## R3b oral: gingiva -> multi-site normal oral -> oral SCC
_data: {"n_classes": 4, "sets": ["gingiva", "oral_normal_multisite", "oral_scc"], "per_set": {"gingiva": {"cells": 10200, "kept": 5620, "donors": 34}, "oral_normal_multisite": {"cells": 5731, "kept": 3675, "donors": 20}, "oral_scc": {"cells": 3460, "kept": 1796, "donors": 12}}}_
_classes: ['endothelial', 'epithelial', 'fibroblast', 'myeloid']_

Target coverage 0.90. coverage = donor-mean coverage on TARGET donors. p_rep_under5 = share of replications more than 5 points below target.

| embedding | method | reps | coverage | cov_sd | set_size | p_rep_under5 | frac_poor_donors | reported_random_split |
|---|---|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 60 | 0.971 | 0.011 | 1.006 | 0.000 | 0.000 | 0.901 |
| emb_scvi | donor_crc | 60 | 0.998 | 0.002 | 1.514 | 0.000 | 0.000 | nan |
| emb_scvi | donor_weighted | 60 | 0.972 | 0.010 | 1.008 | 0.000 | 0.000 | nan |
| emb_tf-exemplar-human | cell_pooled | 60 | 0.977 | 0.013 | 1.142 | 0.000 | 0.000 | 0.902 |
| emb_tf-exemplar-human | donor_crc | 60 | 0.996 | 0.004 | 2.117 | 0.000 | 0.000 | nan |
| emb_tf-exemplar-human | donor_weighted | 60 | 0.976 | 0.013 | 1.120 | 0.000 | 0.000 | nan |
| emb_tf-sapiens | cell_pooled | 60 | 0.978 | 0.014 | 1.084 | 0.000 | 0.000 | 0.903 |
| emb_tf-sapiens | donor_crc | 60 | 0.995 | 0.006 | 1.975 | 0.000 | 0.000 | nan |
| emb_tf-sapiens | donor_weighted | 60 | 0.977 | 0.014 | 1.069 | 0.000 | 0.000 | nan |

Reported (within-calibration random split) minus realised (target donors) coverage for cell_pooled, by embedding:
| embedding | reported_minus_realised |
|---|---|
| emb_scvi | -0.070 |
| emb_tf-exemplar-human | -0.075 |
| emb_tf-sapiens | -0.075 |

Shift screen vs coverage failure (fail = donor-mean coverage >5 points below target). A useful screen has p_fail_given_flagged well above p_fail_given_not_flagged:
| embedding | method_failing | screen_flag_rate | p_fail | p_fail_given_flagged | p_fail_given_not_flagged | head_accuracy_target |
|---|---|---|---|---|---|---|
| emb_scvi | cell_pooled | 1.000 | 0.000 | 0.000 | nan | 0.977 |
| emb_tf-exemplar-human | cell_pooled | 0.900 | 0.000 | 0.000 | 0.000 | 0.972 |
| emb_tf-sapiens | cell_pooled | 0.867 | 0.000 | 0.000 | 0.000 | 0.975 |
| emb_scvi | donor_crc | 1.000 | 0.000 | 0.000 | nan | 0.977 |
| emb_tf-exemplar-human | donor_crc | 0.900 | 0.000 | 0.000 | 0.000 | 0.972 |
| emb_tf-sapiens | donor_crc | 0.867 | 0.000 | 0.000 | 0.000 | 0.975 |

## R4 gingiva composition: normal vs periodontitis (nominal 0.95)
Reference = difference in mean per-donor proportions from all donors' Census labels. Cross-fitted predictions. confusion_proxy is NOT DCATS.

### Overall
| embedding | method | coverage | width |
|---|---|---|---|
| emb_scvi | cell_pooled | 0.882 | 0.019 |
| emb_scvi | confusion_proxy | 0.902 | 0.103 |
| emb_scvi | naive_plugin | 1.000 | 0.086 |
| emb_scvi | ppi_donor | 0.985 | 0.113 |
| emb_tf-sapiens | cell_pooled | 1.000 | 0.019 |
| emb_tf-sapiens | confusion_proxy | 0.919 | 0.100 |
| emb_tf-sapiens | naive_plugin | 1.000 | 0.085 |
| emb_tf-sapiens | ppi_donor | 0.980 | 0.111 |

### Five cell types with the largest reference difference: ['endothelial cell of venule', 'keratinocyte', 'plasma cell']
| embedding | method | coverage | width |
|---|---|---|---|
| emb_scvi | cell_pooled | 1.000 | 0.029 |
| emb_scvi | confusion_proxy | 0.950 | 0.246 |
| emb_scvi | naive_plugin | 1.000 | 0.214 |
| emb_scvi | ppi_donor | 0.958 | 0.273 |
| emb_tf-sapiens | cell_pooled | 1.000 | 0.029 |
| emb_tf-sapiens | confusion_proxy | 0.956 | 0.243 |
| emb_tf-sapiens | naive_plugin | 1.000 | 0.212 |
| emb_tf-sapiens | ppi_donor | 0.958 | 0.271 |

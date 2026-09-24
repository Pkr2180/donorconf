# Development-time simulation results (auto-generated)

**All numbers below come from simulated data with known ground truth.** They test whether the code meets its stated properties and show where methods break. They are not biological evidence and were seen before the real-data protocol was locked (see PROTOCOL.md).

## H1 coverage on new donors (target = 1 - alpha)
_generated 2026-09-21T05:41:01Z | reps=100 | alpha=0.1 | SIMULATION ONLY_

Columns: donor_mean_coverage = coverage averaged with equal donor weight; p_rep_undercovers_5pts = share of replications whose realised coverage fell >5 points below target; frac_donors_poorly_covered = share of donors covered >10 points below target.

### Overall (mean over all scenarios)
| method | donor_mean_coverage | mean_set_size | p_rep_undercovers_5pts | frac_donors_poorly_covered |
|---|---|---|---|---|
| cell_pooled | 0.894 | 2.107 | 0.147 | 0.119 |
| donor_crc | 0.968 | 3.896 | 0.002 | 0.029 |
| donor_weighted | 0.896 | 2.069 | 0.100 | 0.113 |

### By size_difficulty
| size_difficulty | cov / cell_pooled | cov / donor_crc | cov / donor_weighted | size / cell_pooled | size / donor_crc | size / donor_weighted | p_rep_under / cell_pooled | p_rep_under / donor_crc | p_rep_under / donor_weighted | frac_poor / cell_pooled | frac_poor / donor_crc | frac_poor / donor_weighted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| -0.500 | 0.862 | 0.967 | 0.895 | 1.757 | 3.972 | 2.085 | 0.329 | 0.002 | 0.127 | 0.186 | 0.037 | 0.131 |
| 0.000 | 0.898 | 0.970 | 0.898 | 1.862 | 3.666 | 1.857 | 0.057 | 0.000 | 0.043 | 0.080 | 0.013 | 0.080 |
| 0.500 | 0.922 | 0.967 | 0.896 | 2.701 | 4.051 | 2.266 | 0.056 | 0.003 | 0.131 | 0.090 | 0.035 | 0.127 |

### By n_cal_donors
| n_cal_donors | cov / cell_pooled | cov / donor_crc | cov / donor_weighted | size / cell_pooled | size / donor_crc | size / donor_weighted | p_rep_under / cell_pooled | p_rep_under / donor_crc | p_rep_under / donor_weighted | frac_poor / cell_pooled | frac_poor / donor_crc | frac_poor / donor_weighted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6.000 | 0.891 | 1.000 | 0.894 | 2.099 | 6.000 | 2.085 | 0.178 | 0.000 | 0.138 | 0.126 | 0.000 | 0.120 |
| 12.000 | 0.893 | 0.969 | 0.895 | 2.071 | 3.183 | 2.049 | 0.148 | 0.000 | 0.108 | 0.120 | 0.028 | 0.114 |
| 24.000 | 0.899 | 0.936 | 0.899 | 2.150 | 2.506 | 2.073 | 0.116 | 0.006 | 0.056 | 0.109 | 0.058 | 0.104 |

### By donor_sd
| donor_sd | cov / cell_pooled | cov / donor_crc | cov / donor_weighted | size / cell_pooled | size / donor_crc | size / donor_weighted | p_rep_under / cell_pooled | p_rep_under / donor_crc | p_rep_under / donor_weighted | frac_poor / cell_pooled | frac_poor / donor_crc | frac_poor / donor_weighted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.400 | 0.900 | 0.970 | 0.899 | 1.442 | 3.288 | 1.414 | 0.018 | 0.000 | 0.014 | 0.046 | 0.008 | 0.045 |
| 0.800 | 0.895 | 0.968 | 0.896 | 2.075 | 3.904 | 1.995 | 0.161 | 0.001 | 0.117 | 0.128 | 0.031 | 0.124 |
| 1.200 | 0.887 | 0.966 | 0.893 | 2.803 | 4.498 | 2.799 | 0.262 | 0.004 | 0.170 | 0.181 | 0.047 | 0.168 |

### What a within-dataset random cell split would report
reported coverage range across scenarios: 0.899 to 0.903; realised (new-donor) coverage of the same sets ranged 0.838 to 0.942.

- cell_pooled: worst scenario coverage 0.838; scenarios below target-0.05: 3 of 27
- donor_crc: worst scenario coverage 0.933; scenarios below target-0.05: 0 of 27
- donor_weighted: worst scenario coverage 0.886; scenarios below target-0.05: 0 of 27

## H2 label-free donor-level shift screen
_generated 2026-09-21T05:35:28Z | reps=150 | alpha=0.1 | SIMULATION ONLY_

p_fail = probability that realised donor-mean coverage fell >5 points below target. A useful screen has p_fail_given_flagged clearly above p_fail_given_not_flagged.

| shift | magnitude | target | screen_rejection_rate | coverage_no_repair | coverage_after_recal | p_fail_no_repair | p_fail_given_flagged | p_fail_given_not_flagged |
|---|---|---|---|---|---|---|---|---|
| none | 0.000 | 0.900 | 0.067 | 0.945 | 0.939 | 0.000 | 0.000 | 0.000 |
| dataset_offset | 0.500 | 0.900 | 0.087 | 0.938 | 0.933 | 0.000 | 0.000 | 0.000 |
| dataset_offset | 1.000 | 0.900 | 0.187 | 0.931 | 0.926 | 0.000 | 0.000 | 0.000 |
| dataset_offset | 2.000 | 0.900 | 0.680 | 0.895 | 0.896 | 0.067 | 0.088 | 0.021 |
| label_shift | 2.000 | 0.900 | 0.093 | 0.948 | 0.941 | 0.000 | 0.000 | 0.000 |
| label_shift | 5.000 | 0.900 | 0.187 | 0.951 | 0.946 | 0.000 | 0.000 | 0.000 |
| concept_shift | 0.300 | 0.900 | 0.253 | 0.919 | 0.916 | 0.000 | 0.000 | 0.000 |
| concept_shift | 0.600 | 0.900 | 0.280 | 0.871 | 0.876 | 0.260 | 0.262 | 0.259 |

## H3 composition-difference intervals (nominal 0.95)
_generated 2026-09-21T05:43:37Z | reps=150 | alpha=0.05 | SIMULATION ONLY_

confusion_proxy is a confusion-matrix correction standing in for the idea behind DCATS; DCATS itself has NOT been run.

### Overall
| method | ci_coverage | mean_ci_width |
|---|---|---|
| cell_pooled | 0.274 | 0.028 |
| confusion_proxy | 0.907 | 0.205 |
| naive_plugin | 0.851 | 0.117 |
| ppi_donor | 0.965 | 0.194 |

### By donors_per_group
| donors_per_group | cover / cell_pooled | cover / confusion_proxy | cover / naive_plugin | cover / ppi_donor | width / cell_pooled | width / confusion_proxy | width / naive_plugin | width / ppi_donor |
|---|---|---|---|---|---|---|---|---|
| 12.000 | 0.317 | 0.876 | 0.906 | 0.968 | 0.039 | 0.296 | 0.166 | 0.200 |
| 24.000 | 0.276 | 0.921 | 0.866 | 0.963 | 0.027 | 0.186 | 0.109 | 0.193 |
| 48.000 | 0.228 | 0.925 | 0.783 | 0.966 | 0.019 | 0.132 | 0.076 | 0.190 |

### By labeled_per_group
| labeled_per_group | cover / cell_pooled | cover / confusion_proxy | cover / naive_plugin | cover / ppi_donor | width / cell_pooled | width / confusion_proxy | width / naive_plugin | width / ppi_donor |
|---|---|---|---|---|---|---|---|---|
| 4.000 | 0.284 | 0.914 | 0.848 | 0.971 | 0.028 | 0.193 | 0.115 | 0.246 |
| 8.000 | 0.263 | 0.900 | 0.855 | 0.959 | 0.028 | 0.217 | 0.119 | 0.142 |

### By donor_sd
| donor_sd | cover / cell_pooled | cover / confusion_proxy | cover / naive_plugin | cover / ppi_donor | width / cell_pooled | width / confusion_proxy | width / naive_plugin | width / ppi_donor |
|---|---|---|---|---|---|---|---|---|
| 0.500 | 0.333 | 0.903 | 0.845 | 0.963 | 0.028 | 0.137 | 0.089 | 0.175 |
| 1.000 | 0.214 | 0.911 | 0.858 | 0.968 | 0.028 | 0.273 | 0.145 | 0.213 |

### By cell type (0 and 1 have a true difference; 2 is a true null)
| cell_type | cell_pooled | confusion_proxy | naive_plugin | ppi_donor |
|---|---|---|---|---|
| 0.000 | 0.229 | 0.904 | 0.806 | 0.967 |
| 1.000 | 0.251 | 0.901 | 0.801 | 0.963 |
| 2.000 | 0.341 | 0.917 | 0.948 | 0.966 |

## GEO record check (scripts/verify_geo_accessions.py)

- GSE178360: found=True samples=3 taxon=Homo sapiens | Human distal lung maps and lineage hierarchies reveal a bipotent progenitor [sin
- GSE127465: found=True samples=40 taxon=Homo sapiens; Mus musculus | Single cell transcriptomics of human and mouse lung cancers reveals conserved my
- GSE84133: found=True samples=6 taxon=Homo sapiens; Mus musculus | A single-cell transcriptomic map of the human and mouse pancreas reveals inter- 
- GSE180878: found=True samples=16 taxon=Homo sapiens | A human breast transcriptome atlas
- GSE165816: found=True samples=54 taxon=Homo sapiens | Single Cell Transcriptomic Landscape of Diabetic Foot Ulcers

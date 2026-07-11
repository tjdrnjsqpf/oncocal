# Revision / robustness re-analyses (2026-07)

Statistical robustness checks and figure rebuild produced during manuscript revision. These scripts
consume the main-pipeline outputs (`../../tables/`, `../../data/`) and reproduce the exact numbers
reported in the revised manuscript. They do **not** modify the core pipeline in `../../src/`.

| script | purpose | key result |
|---|---|---|
| `s1_diagnose.py`, `s1_unify.py` | Unify the definition of "AlphaMissense-missed" oncogenic positives (deduplicated unique variants, raw score < 0.5) | 1,071 / 3,922 (27.3%) missed; medians phyloP 2.51 vs 7.89, RSA 0.671 vs 0.185 |
| `s2_clustered.py` | Remove pseudo-replication: treat each gene as the independent unit (gene-collapsed Wilcoxon, cluster-robust GEE, gene-cluster bootstrap, mixed-effects logistic) | gene-clustered p = 4.8×10⁻²⁰ (phyloP) / 2.3×10⁻³⁵ (RSA); effect sizes robust across all three methods |
| `s3_threshold_free.py` | Re-test the unsupervised-vs-supervised gap without a fixed threshold (per-tool oncogene AUROC) | unsupervised oncogene AUROC 0.744 vs supervised 0.835 (permutation p = 1×10⁻⁴); not a calibration artifact, but GOF-specific excess is borderline |
| `s4_ascertainment.py` | Control for label ascertainment: does the conservation/exposure signal hold within single label sources and recurrence-matched strata? | effect holds in all strata (gene-clustered GEE p 4×10⁻⁶ – 4×10⁻⁴⁹) — not a labeling artifact |
| `s5_logo_ci.py` | Leave-one-out per-gene calibration (removes in-sample optimism) + bootstrap 95% CI on the AUROC gap | per-gene calibration +0.032 (OncoCal) / +0.016 (AlphaMissense) under within-gene leave-one-out; OncoCal AUROC 0.928 [0.923, 0.933] vs best tool 0.870 [0.862, 0.878], ΔAUROC 0.058 [0.052, 0.064] |
| `s6_fig2_rebuild.py`, `fig2_rebuild.R` | Rebuild Figure 2 on the unified 3,922-variant set and recompute the low-conservation/exposed corner | corner (phyloP < 4 & RSA > 0.4): 81% missed vs 15% elsewhere; `fig2A_scatter_3922.pdf`, `fig2B_corner_3922.pdf` |

Output tables (`s1_unified_mechanism.tsv`, `s2_clustered_stats.tsv`, `s3_threshold_free.tsv`,
`s3_tool_thresholdfree.tsv`, `s5_auroc_ci.tsv`, `s5_logo_bacc.tsv`, `s6_corner_stats.tsv`,
`fig2_master_3922.tsv`) accompany each script.

## Run
```bash
conda activate mutfunc
python s1_unify.py && python s2_clustered.py && python s3_threshold_free.py
python s4_ascertainment.py && python s5_logo_ci.py && python s6_fig2_rebuild.py
Rscript fig2_rebuild.R      # figures
```

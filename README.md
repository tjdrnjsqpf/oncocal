# OncoCal — cancer-calibrated missense oncogenicity benchmark

Code and openly licensed data for *"An openly licensed benchmark and per-gene calibration map for missense pathogenicity
predictors on activating cancer drivers"* (Sung-Gwon Lee).

Germline-trained missense predictors (AlphaMissense, REVEL, …) systematically under-score
activating gain-of-function cancer drivers (e.g., JAK2 V617F). This repository contains the
pipeline, the openly licensed benchmark, the per-gene calibration map, and the figure code.

## Repository layout
```
src/         analysis pipeline (data build, features, model, evaluation)
scripts/     data download / run scripts (scripts/revision/ = robustness re-analyses, 2026-07)
figures/     R (ggplot2) scripts to reproduce all figures (input = tables/, output = figures_R/)
tables/      result and summary tables (open data; see Data below; tables/r1/ = revision outputs)
data/        README only — restricted/large source data fetched from providers
environment.yml   conda environment, pinned to the versions that produced the results
environment-lock.txt  full resolved environment (linux-64, build strings)
```

## Data
- **Labels (released here):** CIViC (CC0), ClinVar (public), cancerhotspots (open), gnomAD (open).
- **COSMIC:** only **summary statistics** (per-variant sample counts) are included; raw COSMIC is
  **not redistributed** per its academic licence — obtain from https://www.cosmickb.org.
- **Restricted / large (fetch from source):** dbNSFP 5.3.1a (https://www.dbnsfp.org),
  AlphaFold DB (https://alphafold.ebi.ac.uk), ProteinGym (https://proteingym.org).
- OncoKB was **not** used.

## Reproduce
```bash
conda env create -f environment.yml && conda activate mutfunc
# 1. fetch restricted data into data/ (see data/README.md)
python -m src.data.build_dataset --recurrence-min 10        # -> data/processed/dataset.parquet
python -m src.data.make_group_split --level gene
python -m src.models.train_baseline                          # OncoCal, gene-group OOF
# figures
for f in figures/fig*.R; do Rscript "$f"; done   # run from the repo root
```

## Revision robustness analyses

`scripts/revision/` reproduces the statistical robustness checks and the rebuilt Figure 2 from the
first revision (gene-clustered statistics, threshold-free comparison, ascertainment controls,
leave-one-out calibration, bootstrap CIs). See `scripts/revision/README.md`.

`scripts/revision_r1/` holds the analyses added for the PLOS ONE revision, with their outputs in
`tables/r1/`:

| script | what it answers |
|---|---|
| `r1a_label_source.py` | oncogene-vs-TSG gap recomputed per label source (CIViC / cancerhotspots / COSMIC recurrence) |
| `r1b_direct_mechanism.py` | direct CIViC mechanism labels vs the gene-role proxy, and their concordance |
| `r1c_matched_negatives.py` | discrimination under progressively harder negative sets (rare-benign, gene-matched) |
| `r1d_role_census_core.py` | gene-role census, the core-role restriction, and the mechanism statistics under it |
| `r1e_pergene_core.py` | per-gene thresholds; quantifies the off-by-one in the original threshold function |
| `r1f_pergene_ci.py` | gene-clustered bootstrap CI on the per-gene calibration gain |
| `r1g_auroc_cluster_ci.py` | AUROC CIs under gene-clustered vs variant-level resampling |
| `r2a_rescue_repseed.py` | two-directional rescue accounting at matched operating points |
| `r2b_tool_gap_twostage.py` | per-tool gap CIs, both deduplication conventions |
| `r2c_family_contrast.py` | pre-specified unsupervised-vs-supervised family contrast; per-tool BH/Holm adjustment |
| `r2d_marginal_vs_gap_ci.py` | per-axis gene-clustered CIs for each predictor's oncogene and TSG AUROC (Fig 1A error bars, S1 Table) |
| `c2_within_gene_mechanism.py` | gene-stratified conditional logistic regression: the conservation and RSA contrasts estimated strictly within genes (S2 Table) |
| `r2a_rescue_am564.py` | rescue accounting at AlphaMissense's published 0.564 boundary, with a gene-clustered bootstrap |

`r2a_rescue_repseed.py` and `r2a_rescue_am564.py` run from the released tables alone:

```bash
python scripts/revision_r1/r2a_rescue_am564.py      # 0.564 cut, as reported
python scripts/revision_r1/r2a_rescue_repseed.py   # earlier 0.5 cut, superseded
```

`r2a_rescue_repseed.py` is pinned to the earlier 0.5 cut and prints the superseded
315 / 746 / net 227 figures; `r2a_rescue_am564.py` reproduces the 336 / 811 (41.4%),
net 250 figures quoted in the paper, with the gene-clustered interval.

The others additionally need `data/processed/dataset.parquet`, which is rebuilt from the
provider-restricted sources listed under **Data** (it is not redistributable).

### Threshold function correction

`youden_threshold()` as used for the original submission returned the score immediately *below*
the J-optimal cut — a threshold for `score > T` — while every call site applied it as `score >= T`.
On sparse per-gene data this shifts thresholds downward (median ~0.06, up to ~0.44 in simulation).
`scripts/revision_r1/r1e_pergene_core.py` carries both the corrected and the original function and
reports the difference (`tables/r1/r1e_threshold_fix_impact.tsv`); all per-gene numbers in the
revised manuscript use the corrected one.

## Citation
Lee S-G. An openly licensed benchmark and per-gene calibration map for missense pathogenicity
predictors on activating cancer drivers. Manuscript under review, 2026.
Archived at Zenodo (concept DOI): 10.5281/zenodo.20741850.

## License
Code: MIT (see LICENSE). Data tables: per source licence (CIViC CC0; ClinVar/cancerhotspots/gnomAD open;
COSMIC summary statistics only). dbNSFP, AlphaFold, ProteinGym remain under their providers' terms.

### Which tables the paper reports

Two tables in `tables/` are kept at their released state for provenance and are **not** the
ones the revised manuscript reports:

- `per_gene_calibration.tsv` is the 118-gene map from v1.1.0, computed before the
  `youden_threshold()` correction above and before the core-role restriction. The map the
  paper reports is `tables/r1/r1e_pergene_table_corrected.tsv` (104 genes; S6 Table), also
  copied to `tables/per_gene_calibration_corrected.tsv`.
- `external_validation_dms.tsv` previously carried `SRC_HUMAN_Ahler_CD_2019`, which
  duplicated `SRC_HUMAN_Ahler_2019` exactly and is not a ProteinGym v1.3 assay (v1.3 lists
  three SRC assays). It is removed as of v1.2.1, so the file now matches S8 Table:
  14 scored assays over 11 proteins.

### Environment

`environment.yml` pins the versions used for the released results (Python 3.11.15, numpy
2.4.6, pandas 3.0.3, scikit-learn 1.9.0, statsmodels 0.14.6, LightGBM 4.6.0).
`environment-lock.txt` is the full resolved environment including transitive dependencies
and build strings, and `figures/R-sessionInfo.txt` records the R session used for the
figures (R 4.4.2, ggplot2 4.0.2).

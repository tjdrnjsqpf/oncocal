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
environment.yml   conda environment (Python 3.11)
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

`r2a_rescue_repseed.py` runs from the released tables alone:

```bash
python scripts/revision_r1/r2a_rescue_repseed.py
```

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

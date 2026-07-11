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
figures/     R (ggplot2) scripts to reproduce all figures (input = tables/)
tables/      result and summary tables (open data; see Data below)
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
cd figures && for f in fig*.R; do Rscript "$f"; done
```

## Revision robustness analyses
`scripts/revision/` reproduces the statistical robustness checks and the rebuilt Figure 2 from the
revised manuscript (gene-clustered statistics, threshold-free comparison, ascertainment controls,
leave-one-out calibration, bootstrap CIs). See `scripts/revision/README.md`.

## Citation
Lee S-G. An openly licensed benchmark and per-gene calibration map for missense pathogenicity
predictors on activating cancer drivers. Manuscript under review, 2026.
Archived at Zenodo (concept DOI): 10.5281/zenodo.20741850.

## License
Code: MIT (see LICENSE). Data tables: per source licence (CIViC CC0; ClinVar/cancerhotspots/gnomAD open;
COSMIC summary statistics only). dbNSFP, AlphaFold, ProteinGym remain under their providers' terms.

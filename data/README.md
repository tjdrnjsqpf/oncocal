# Data

This repository ships **openly licensed labels and derived summary tables** in `../tables/`.
Large or restricted source data are **not** included; fetch them from the providers below and
place under `data/` to reproduce the pipeline.

| Source | Used for | License | Get from |
|--------|----------|---------|----------|
| dbNSFP 5.3.1a (GRCh38) | 49 predictor scores + conservation | academic, free | https://www.dbnsfp.org/download (academic 'a' build) |
| COSMIC v104 (GRCh38) | gene panel (CGC) + recurrence | academic; **no redistribution** | https://www.cosmickb.org (only summary stats shared here) |
| AlphaFold DB v6 | RSA / pLDDT | CC BY 4.0 | https://alphafold.ebi.ac.uk |
| ProteinGym v1.3 | external DMS validation | open | https://proteingym.org |
| CIViC (nightly) | oncogenic labels | CC0 | https://civicdb.org (also in ../tables) |
| ClinVar | benign labels | public | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ |
| cancerhotspots v2 | hotspot labels | open | https://www.cancerhotspots.org |
| gnomAD (via dbNSFP) | common-variant negatives | open | included in dbNSFP |

OncoKB was **not** used (licence prohibits model training).

Expected layout after download:
```
data/external/dbnsfp/dbNSFP5.3.1a_grch38.gz
data/external/cosmic/Cosmic_CancerGeneCensus_*.tsv  (+ recurrence files)
data/external/alphafold/  (per-protein PDBs)
data/external/proteingym/
```

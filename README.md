# Bilingualism and Cognition across Adolescence: A Longitudinal Analysis of the ABCD Study

This repository contains the R analysis pipeline for a longitudinal study examining bilingual language experience and cognitive development in children and adolescents using data from the [Adolescent Brain Cognitive Development (ABCD) Study](https://abcdstudy.org/), Release 6.0.

## Overview

The pipeline examines associations between bilingual language use and a broad range of cognitive outcomes — including NIH Toolbox measures, stop-signal inhibitory control, emotional Stroop interference, and emotional n-back working memory — across seven annual assessment waves (ages ~9–18). Analyses use linear mixed-effects models (LMMs) and generalized additive mixed models (GAMMs) with Bayesian multiple imputation via [BLIMP](https://www.appliedmissingdata.com/blimp).

Key features of the analysis include:

- A continuous **Bilingual Use** composite score (0 = exclusively non-English across friends and family; 8 = exclusively English) derived from ABCD acculturation items
- Longitudinal **growth curve modeling** with participant, family, and site as random effects
- A **"higher = better" sign convention** applied uniformly across all outcomes, so that positive bilingual_use coefficients consistently indicate a monolingual advantage
- **Generalized additive mixed models (GAMMs)** with penalized regression splines for age, used where AIC comparisons favor nonlinear trajectories
- **Equivalence testing (TOST)** for main effects of bilingual use
- **FDR correction** (Benjamini-Hochberg) across all 39 tests (13 outcomes × 3 predictors)
- Covariate control for sex, family education, household income, NIH Toolbox Oral Reading Recognition, and the top 10 genetic ancestry principal components

## Pipeline Structure

The script is organized into seven phases, with a required manual stop between Phases 4 and 5 to run BLIMP externally.

| Phase | Description |
|-------|-------------|
| 1 | Setup, data ingestion via `NBDCtools`, and variable selection |
| 2 | Variable construction: bilingual use score, monolingual recoding, balanced panel |
| 3 | Pre-processing: outlier cleaning (SST, N-back, Stroop), z-scoring, cost score computation |
| 4 | Export to BLIMP for multiple imputation (4 separate imputation runs) |
| — | **STOP: Run BLIMP externally before proceeding** |
| 5 | Import BLIMP output, post-imputation scoring, LMM modeling, CI tables |
| 6 | Equivalence testing and forest plot |
| 7 | GAMM fitting, AIC model comparison, table extraction, final trajectory plots |

## Data Access

This script requires access to **ABCD Study Release 6.0** data in Parquet format. ABCD data are not publicly available and must be obtained through the [NIMH Data Archive (NDA)](https://nda.nih.gov/). Researchers must have an approved NDA data use agreement.

Data are loaded using the `NBDCtools` and `NBDCtoolsData` packages, which provide utilities for working with ABCD NDA releases. See the [NBDCtools documentation](https://nbdc-datasharing.github.io/NBDCtools/) for setup instructions.

## Requirements

### R Packages

```r
# Data access and wrangling
NBDCtools
NBDCtoolsData
dplyr
tidyr
readr
purrr
tibble

# Multiple imputation
mitml

# Modeling
lme4
mgcv
MuMIn

# Visualization
ggplot2
viridis
patchwork

# Tables
xtable

# Parallelization
parallel
```

### External Software

- [**BLIMP**](https://www.appliedmissingdata.com/blimp) (v3 or later) — used for Bayesian multiple imputation of the longitudinal dataset. The script exports four CSV files at the end of Phase 4 and imports the corresponding stacked imputation output files at the start of Phase 5. BLIMP must be run manually between these phases.

## Setup

1. Clone this repository.
2. Set the two path variables near the top of the script:

```r
dir_abcd <- "/path/to/abcd/6_0/phenotype"  # directory containing ABCD parquet files
dir_proj <- "/path/to/project/data/"        # project output directory for CSVs, RData, and plots
```

3. Run Phases 1–4 to generate the BLIMP input CSVs.
4. Run BLIMP using your `.imp` syntax files on the four exported CSVs.
5. Place the stacked BLIMP output files in `dir_proj` and run Phases 5–7.

## Outputs

- **LaTeX tables** for LMM and GAMM results (printed to console; paste into your manuscript)
- **TIFF and PNG trajectory plots** for each outcome (original units, ages 9–18, by bilingual use level)
- **Equivalence test forest plot** (`bilingual_use_equiv.png`)
- **RT–accuracy coupling figures** (`sat_combined_by_age.tiff`, `sat_bilingual_by_age.tiff`)
- **Saved model objects** (`bilingual_model_objects.RData`, `gamm_fit_objects.RData`) for downstream use without refitting

## Citation

If you use this code, please cite the associated paper (citation to be added upon publication) and the ABCD Study:

> Volkow, N. D., et al. (2018). The conception of the ABCD study: From substance use to a broad NIH collaboration. *Developmental Cognitive Neuroscience*, 32, 4–7. https://doi.org/10.1016/j.dcn.2017.10.002

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Anthony Steven Dick
adick@fiu.edu
Florida International University  
[anthonystevendick](https://github.com/anthonystevendick)

################################################################################
# ABCD 6.0 DATA PIPELINE: BILINGUALISM & COGNITION
#
# Phase 1: Setup, Ingestion & Structuring
# Phase 2: Variable Construction (bilingual use score; monolingual fix applied)
# Phase 3: Pre-processing (RT/accuracy outlier cleaning, z-scoring, cost scores)
# Phase 4: Export to BLIMP for multiple imputation
# [STOP HERE MANUALLY TO RUN BLIMP]
# Phase 5: Import BLIMP output, post-imputation scoring, LMM modeling, tables
# Phase 6: Equivalence testing (TOST) for bilingual use main effects
# Phase 7: GAMM models (hybrid smooth/linear age), table extraction, final plots
################################################################################
# ==============================================================================
# PHASE 1: SETUP & INGESTION
# ==============================================================================
suppressPackageStartupMessages({
  library(NBDCtools)
  library(NBDCtoolsData)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(mitml)
  library(lme4)
  library(ggplot2)
  library(viridis)
})
# --- Paths (Update to match your local setup) ---
# dir_abcd: directory containing the ABCD 6.0 phenotype parquet files
# dir_proj: project output directory where CSVs, RData files, and plots are saved
dir_abcd <- "/path/to/abcd/6_0/phenotype"
dir_proj <- "/path/to/project/data/"
# --- Variable Definitions ---
pcs <- sprintf("gn_y_popstruct_pc__%02d", 1:10)
vars_main <- c(
  # Cognition (NIHTB + WISC)
  "nc_y_nihtb__lswmt__uncor_score",       # List Sorting Working Memory Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_nihtb__readr__uncor_score",       # Oral Reading Recognition Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_nihtb__picvcb__uncor_score",      # Picture Vocabulary Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_nihtb__comp__cryst__uncor_score", # NIH Toolbox: Crystallized composite - Uncorrected standard score
  "nc_y_nihtb__comp__fluid__uncor_score", # NIH Toolbox: Cognition fluid composite - Uncorrected standard score
  "nc_y_nihtb__crdst__uncor_score",       # Dimensional Change Card Sort Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_nihtb__flnkr__uncor_score",       # Flanker Inhibitory Control and Attention Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_nihtb__pttcp__uncor_score",       # Pattern Comparison Processing Speed Task (NIH Toolbox): Uncorrected Standard Score
  "nc_y_wisc__scaled_score",              # WISC-V Matrix Reasoning: total Scaled Score
  # Emotional Stroop (EST)
  "nc_y_est_acc",                         # EST: Accuracy rate (proportion correct) across all correct trials
  "nc_y_est_rt",                          # EST: Mean reaction time averaged across all correct trials
  "nc_y_est__congr_acc",                  # EST: Accuracy rate in congruent trials (across conditions)
  "nc_y_est__congr_rt",                   # EST: Mean reaction time across all correct congruent trials (across conditions)
  "nc_y_est__incongr_acc",               # EST: Accuracy rate in incongruent trials (across conditions)
  "nc_y_est__incongr_rt",               # EST: Mean reaction time across all correct incongruent trials (across conditions)
  # N-back (Behavioral tfMRI) - CORRECTED FOR 6.0
  "mr_y_tfmri__nback__beh__0b_acc",      # N-Back: Correct response rate to 0-back stimuli (runs 1 & 2)
  "mr_y_tfmri__nback__beh__0b_rt",       # N-Back: Avg reaction time for correct 0-back responses (runs 1 & 2)
  "mr_y_tfmri__nback__beh__2b_acc",      # N-Back: Correct response rate to 2-back stimuli (runs 1 & 2)
  "mr_y_tfmri__nback__beh__2b_rt",       # N-Back: Avg reaction time for correct 2-back responses (runs 1 & 2)
  "mr_y_tfmri__nback__beh__qc_indicator", # N-Back: Acceptable task performance indicator (replaces perf_ok, match_indicator,
  # tdiff_ign_indicator, imgincl_nback_include)
  # Bilingual / Acculturation (Youth)
  "fc_y_aclt_001",                        # How well do you speak English? [Youth]
  "fc_y_aclt_002",                        # Besides English, do you speak or understand another language or dialect? [Youth]
  "fc_y_aclt_002__01",                    # What other language or dialect do you speak or understand (besides English)? [Youth]
  "fc_y_aclt_002__02",                    # How well do you speak the other language? [Youth]
  "fc_y_aclt_002__03",                    # What language do you speak with most of your friends? [Youth]
  "fc_y_aclt_002__04",                    # What language do you speak with most of your family? [Youth]
  # Parent Reports
  "fc_p_aclt_001",                        # How well do you speak English? [Parent]
  "fc_p_aclt_002",                        # Besides English, do you speak or understand another language or dialect? [Parent]
  "fc_p_aclt_002__01",                    # What other language or dialect do you speak or understand (besides English)? [Parent]
  "fc_p_aclt_002__02",                    # How well do you speak the other language? [Parent]
  "fc_p_aclt_002__03",                    # What language do you speak with most of your friends? [Parent]
  "fc_p_aclt_002__04",                    # What language do you speak with most of your family? [Parent]
  "ab_p_demo__ntvlang_003",              # Composition of English vs. other languages spoken in house after child's birth
  "ab_p_demo__ntvlang_004",              # Does your child attend a dual-language or language immersion program at school?
  # SST (Behavioral tfMRI) - CORRECTED FOR 6.0
  "mr_y_tfmri__sst__beh__ssrt_intgr",    # SST: Stop Signal Reaction Time, integration estimation
  "mr_y_tfmri__sst__beh__ssrt_mean",     # SST: Stop Signal Reaction Time, mean estimation
  "mr_y_tfmri__sst__beh__violat_indicator", # SST: Race model violators where Stop Fail RT > Go RT
  "mr_y_tfmri__sst__beh__qc_indicator",  # SST: Acceptable task performance indicator (all trial types must meet criteria)
  "mr_y_tfmri__sst__beh__coderr_indicator" # SST: Task coding error
  # Removed: match_indicator, tdiff_ign_indicator — gone in 6.0
)

vars_covars <- c(
  "ab_p_demo_age",                        # Youth's age at data collection
  "ab_g_stc__cohort_sex",                 # Participant's sex
  "ab_g_dyn__cohort_edu__cgs",            # Highest education across caregivers
  "ab_g_dyn__cohort_income__hhold__3lvl", # Household income - 3 levels
  "ab_g_stc__cohort_ethn",               # Ethnicity (Hispanic or not Hispanic)
  pcs,                                    # Genetic ancestry principal components 1–10
  # Design variables
  "ab_g_stc__design_id__fam",            # Family ID (participants in same family share ID)
  "ab_g_dyn__design_site",               # Assessment site
  "ab_g_dyn__design_mr__serial",         # MRI scanner serial number (hashed)
  "ab_g_dyn__design_mr__manufact",       # MRI scanner manufacturer
  "ab_g_dyn__design_mr__software"        # MRI scanner software version
)
vars_all <- unique(c(vars_main, vars_covars))
# --- Check Files ---
dd6 <- NBDCtools::get_dd("abcd", release = "6.0")
tbls_needed <- dd6 %>%
  filter(name %in% vars_all) %>%
  distinct(table_name) %>%
  pull()
parquet_files <- list.files(dir_abcd, pattern = "\\.parquet$", full.names = FALSE)
missing_tbls <- setdiff(paste0(tbls_needed, ".parquet"), parquet_files)
if (length(missing_tbls)) warning("Missing parquet files:\n", paste(" -", missing_tbls, collapse = "\n"))
# --- Load Data ---
ds6_raw <- NBDCtools::create_dataset(
  dir_data          = dir_abcd,
  study             = "abcd",
  vars              = vars_main,
  vars_add          = vars_covars,
  categ_to_factor   = TRUE,
  add_labels        = TRUE,
  value_to_label    = FALSE,
  value_to_na       = TRUE,
  time_to_hms       = FALSE,
  remove_empty_rows = TRUE
)
# Filter for annual sessions
ds6_annual <- ds6_raw %>%
  NBDCtools::filter_events_abcd("annual") %>%
  NBDCtools::filter_empty_rows()
# ==============================================================================
# PHASE 2: VARIABLE CONSTRUCTION & BALANCING
# ==============================================================================
# 1. Enforce Session Factors
session_priority <- c(
  "ses-00A", "ses-01A", "ses-02A",
  "ses-03A", "ses-04A", "ses-05A", "ses-06A"
)
ds6_annual <- ds6_annual %>%
  mutate(session_id = factor(session_id, levels = session_priority))
# 2. Create Balanced Structure (Skeleton)
all_ids <- ds6_annual %>%
  distinct(participant_id) %>%
  arrange(participant_id)
all_sessions <- tibble(session_id = factor(session_priority, levels = session_priority))
skeleton <- tidyr::expand_grid(
  participant_id = all_ids$participant_id,
  session_id     = all_sessions$session_id
)
ds6_balanced <- skeleton %>%
  left_join(ds6_annual, by = c("participant_id", "session_id"))
# 3. Helper Function: Create "_complete" (Baseline/Best Available) vars
make_complete_by_id <- function(dat, var_name) {
  wide <- dat %>%
    select(participant_id, session_id, !!sym(var_name)) %>%
    distinct(participant_id, session_id, .keep_all = TRUE) %>%
    pivot_wider(names_from = session_id, values_from = !!sym(var_name))

  cname <- paste0(var_name, "_complete")

  wide %>%
    mutate(!!cname := coalesce(
      `ses-02A`, `ses-00A`, `ses-01A`,
      `ses-03A`, `ses-04A`, `ses-05A`, `ses-06A`
    )) %>%
    select(participant_id, !!sym(cname))
}
# 4. Process Gate Variables (Bilingual Logic)
vars_gate_base <- c(
  "fc_y_aclt_002", "fc_y_aclt_002__02", "fc_y_aclt_002__03",
  "fc_y_aclt_002__04", "ab_p_demo__ntvlang_003"
)
complete_by_id <- make_complete_by_id(ds6_balanced, vars_gate_base[1])
for (v in vars_gate_base[-1]) {
  complete_by_id <- complete_by_id %>%
    left_join(make_complete_by_id(ds6_balanced, v), by = "participant_id")
}
# *** UPSTREAM FIX: Monolinguals set to 5 ***
complete_by_id <- complete_by_id %>%
  mutate(across(ends_with("_complete"), ~ as.numeric(as.character(.)))) %>%
  mutate(
    fc_y_aclt_002__02_complete = ifelse(fc_y_aclt_002_complete == 0, 5, fc_y_aclt_002__02_complete),
    fc_y_aclt_002__03_complete = ifelse(fc_y_aclt_002_complete == 0, 5, fc_y_aclt_002__03_complete),
    fc_y_aclt_002__04_complete = ifelse(fc_y_aclt_002_complete == 0, 5, fc_y_aclt_002__04_complete)
  )
ds6_balanced <- ds6_balanced %>% left_join(complete_by_id, by = "participant_id")
# 5. Process Fixed Covariates (Create *_complete versions)
vars_ses02_fixed <- c(
  "ab_g_stc__cohort_sex", "ab_g_dyn__cohort_edu__cgs",
  "ab_g_dyn__cohort_income__hhold__3lvl", "ab_g_stc__cohort_ethn",
  "nc_y_nihtb__readr__uncor_score", "nc_y_nihtb__comp__cryst__uncor_score",
  "nc_y_nihtb__comp__fluid__uncor_score",
  paste0("gn_y_popstruct_pc__", sprintf("%02d", 1:10))
)
ses02_fixed_by_id <- make_complete_by_id(ds6_balanced, vars_ses02_fixed[1])
for (v in vars_ses02_fixed[-1]) {
  ses02_fixed_by_id <- ses02_fixed_by_id %>%
    left_join(make_complete_by_id(ds6_balanced, v), by = "participant_id")
}
ds6_balanced <- ds6_balanced %>% left_join(ses02_fixed_by_id, by = "participant_id")
# 5b. Family ID
fam_complete_by_id <- make_complete_by_id(ds6_balanced, "ab_g_stc__design_id__fam")
ds6_balanced <- ds6_balanced %>% left_join(fam_complete_by_id, by = "participant_id")
# 6. Fill Design Site (LOCF)
ds6_balanced <- ds6_balanced %>%
  group_by(participant_id) %>%
  arrange(session_id) %>%
  fill(ab_g_dyn__design_site, .direction = "downup") %>%
  ungroup()
site_complete <- make_complete_by_id(ds6_balanced, "ab_g_dyn__design_site")
ds6_balanced <- ds6_balanced %>%
  left_join(site_complete, by = "participant_id") %>%
  mutate(ab_g_dyn__design_site = coalesce(
    ab_g_dyn__design_site,
    ab_g_dyn__design_site_complete
  )) %>%
  select(-ab_g_dyn__design_site_complete)
# 7. Clean SST Variables - CORRECTED FOR 6.0
# match_indicator and tdiff_ign_indicator removed; qc_indicator subsumes them
ds6_balanced <- ds6_balanced %>%
  mutate(
    is_bad_sst = (as.character(mr_y_tfmri__sst__beh__qc_indicator) == "0") |
      (as.character(mr_y_tfmri__sst__beh__violat_indicator) == "1") |
      (as.character(mr_y_tfmri__sst__beh__coderr_indicator) == "1"),
    mr_y_tfmri__sst__beh__ssrt_intgr_cleaned = ifelse(is_bad_sst, NA,
      mr_y_tfmri__sst__beh__ssrt_intgr
    ),
    mr_y_tfmri__sst__beh__ssrt_mean_cleaned = ifelse(is_bad_sst, NA,
      mr_y_tfmri__sst__beh__ssrt_mean
    )
  )
# 8. N-Back Data Cleaning - CORRECTED FOR 6.0
# perf_ok, match_indicator, tdiff_ign_indicator, and imgincl_nback_include all
# removed; qc_indicator subsumes all of them
ds6_balanced <- ds6_balanced %>%
  mutate(
    is_bad_nback = (as.character(mr_y_tfmri__nback__beh__qc_indicator) == "0"),

    # Apply QC Indicator
    mr_y_tfmri__nback__beh__0b_acc_cleaned = ifelse(is_bad_nback, NA, mr_y_tfmri__nback__beh__0b_acc),
    mr_y_tfmri__nback__beh__0b_rt_cleaned = ifelse(is_bad_nback, NA, mr_y_tfmri__nback__beh__0b_rt),
    mr_y_tfmri__nback__beh__2b_acc_cleaned = ifelse(is_bad_nback, NA, mr_y_tfmri__nback__beh__2b_acc),
    mr_y_tfmri__nback__beh__2b_rt_cleaned = ifelse(is_bad_nback, NA, mr_y_tfmri__nback__beh__2b_rt),

    # --- NEW: Apply Physiological RT Cutoffs (150ms - 3000ms) ---
    mr_y_tfmri__nback__beh__0b_rt_cleaned = ifelse(
      mr_y_tfmri__nback__beh__0b_rt_cleaned < 150 | mr_y_tfmri__nback__beh__0b_rt_cleaned > 3000,
      NA, mr_y_tfmri__nback__beh__0b_rt_cleaned
    ),
    mr_y_tfmri__nback__beh__2b_rt_cleaned = ifelse(
      mr_y_tfmri__nback__beh__2b_rt_cleaned < 150 | mr_y_tfmri__nback__beh__2b_rt_cleaned > 3000,
      NA, mr_y_tfmri__nback__beh__2b_rt_cleaned
    )
  )

cat("N-back rows cleaned:", sum(ds6_balanced$is_bad_nback, na.rm = TRUE), "\n")
# 9. Drop problematic raw variables to save memory
drop_vars <- c(
  "nc_y_wisc__scaled_score", "ab_p_demo__ntvlang_003", "ab_p_demo__ntvlang_004",
  "ab_g_stc__cohort_sex", "ab_g_stc__cohort_ethn",
  paste0("gn_y_popstruct_pc__", sprintf("%02d", 1:10)),
  "ab_g_stc__design_id__fam"
)
ds6_balanced <- ds6_balanced %>% select(-any_of(drop_vars))
# ==============================================================================
# PHASE 3: PRE-PROCESSING (OUTLIER CLEANING, Z-SCORING, COST SCORES)
# ==============================================================================
# We maintain two clean dataframes:
#   ds6_clean_sst      — base pipeline; retains raw N-back cleaned columns but
#                        does NOT include z-scored N-back variables. Used for
#                        Runs 2, 3, and 4. Column structure matches those BLIMP
#                        input files exactly.
#   ds6_clean_combined — extends ds6_clean_sst with z-scored N-back variables
#                        and load cost difference scores. Used for Run 1 only.

char_vars <- ds6_balanced %>%
  select(where(is.character)) %>%
  names()

ds6_clean_sst <- ds6_balanced %>%
  mutate(
    session_id = as.integer(session_id),
    across(all_of(char_vars), ~ as.integer(factor(.x)))
  ) %>%
  select(-is_bad_sst, -is_bad_nback)

# ------------------------------------------------------------------------------
# CLEAN & CALCULATE INTERFERENCE COSTS
# ------------------------------------------------------------------------------
ds6_clean_sst <- ds6_clean_sst %>%
  mutate(
    nc_y_est_rt = ifelse(nc_y_est_rt > 2000 | nc_y_est_rt < 200, NA, nc_y_est_rt),
    nc_y_est_acc = ifelse(nc_y_est_acc > 1 | nc_y_est_acc < 0, NA, nc_y_est_acc),
    nc_y_est__congr_rt = ifelse(nc_y_est__congr_rt > 2000 | nc_y_est__congr_rt < 200,
      NA, nc_y_est__congr_rt
    ),
    nc_y_est__congr_acc = ifelse(nc_y_est__congr_acc > 1 | nc_y_est__congr_acc < 0,
      NA, nc_y_est__congr_acc
    ),
    nc_y_est__incongr_rt = ifelse(nc_y_est__incongr_rt > 2000 | nc_y_est__incongr_rt < 200,
      NA, nc_y_est__incongr_rt
    ),
    nc_y_est__incongr_acc = ifelse(nc_y_est__incongr_acc > 1 | nc_y_est__incongr_acc < 0,
      NA, nc_y_est__incongr_acc
    ),
    nc_y_est_interference_rt = nc_y_est__incongr_rt - nc_y_est__congr_rt,
    nc_y_est_interference_acc = nc_y_est__congr_acc - nc_y_est__incongr_acc
  )

# ------------------------------------------------------------------------------
# DIAGNOSTIC HISTOGRAMS (run before applying SSRT cutoffs)
# ------------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
# Note: ggplot2 is already loaded in Phase 1; patchwork is an additional dependency
# loaded here. Repeating library() calls is harmless but serves as a local reminder
# of dependencies for readers working through individual sections. +
  geom_histogram(bins = 100, fill = "steelblue", color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(150, 3000), color = "red", linetype = "dashed") +
  labs(title = "N-Back 0-back RT", x = "RT (ms)", y = "Count") +
  theme_minimal()
p2 <- ggplot(ds6_clean_sst, aes(x = mr_y_tfmri__nback__beh__2b_rt_cleaned)) +
  geom_histogram(bins = 100, fill = "steelblue", color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(150, 3000), color = "red", linetype = "dashed") +
  labs(title = "N-Back 2-back RT", x = "RT (ms)", y = "Count") +
  theme_minimal()
p3 <- ggplot(ds6_clean_sst, aes(x = mr_y_tfmri__sst__beh__ssrt_intgr_cleaned)) +
  geom_histogram(bins = 100, fill = "darkorange", color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(50, 700), color = "red", linetype = "dashed") +
  labs(title = "SSRT (Integration)", x = "SSRT (ms)", y = "Count") +
  theme_minimal()
p4 <- ggplot(ds6_clean_sst, aes(x = mr_y_tfmri__sst__beh__ssrt_mean_cleaned)) +
  geom_histogram(bins = 100, fill = "darkorange", color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(50, 700), color = "red", linetype = "dashed") +
  labs(title = "SSRT (Mean)", x = "SSRT (ms)", y = "Count") +
  theme_minimal()
print((p1 | p2) / (p3 | p4))

# ------------------------------------------------------------------------------
# APPLY SSRT CUTOFFS (after visual inspection confirms bounds are appropriate)
# ------------------------------------------------------------------------------
ds6_pre_ssrt_cutoff <- ds6_clean_sst

ds6_clean_sst <- ds6_clean_sst %>%
  mutate(
    mr_y_tfmri__sst__beh__ssrt_intgr_cleaned = ifelse(
      mr_y_tfmri__sst__beh__ssrt_intgr_cleaned < 50 |
        mr_y_tfmri__sst__beh__ssrt_intgr_cleaned > 700,
      NA, mr_y_tfmri__sst__beh__ssrt_intgr_cleaned
    ),
    mr_y_tfmri__sst__beh__ssrt_mean_cleaned = ifelse(
      mr_y_tfmri__sst__beh__ssrt_mean_cleaned < 50 |
        mr_y_tfmri__sst__beh__ssrt_mean_cleaned > 700,
      NA, mr_y_tfmri__sst__beh__ssrt_mean_cleaned
    )
  )

cat("SSRT intgr excluded:", sum(is.na(ds6_clean_sst$mr_y_tfmri__sst__beh__ssrt_intgr_cleaned)) -
  sum(is.na(ds6_pre_ssrt_cutoff$mr_y_tfmri__sst__beh__ssrt_intgr_cleaned)), "\n")
cat("SSRT mean excluded:", sum(is.na(ds6_clean_sst$mr_y_tfmri__sst__beh__ssrt_mean_cleaned)) -
  sum(is.na(ds6_pre_ssrt_cutoff$mr_y_tfmri__sst__beh__ssrt_mean_cleaned)), "\n")

# ------------------------------------------------------------------------------
# FINALIZE ds6_plotting_ref
# Reassign from ds6_clean_sst so SSRT cutoffs are included, then flip RT
# variables so higher = better (faster), matching score_run1/score_run2
# sign conventions used in modeling.
#
# Variables flipped (higher raw = slower = worse -> flip so higher = faster):
#   mr_y_tfmri__sst__beh__ssrt_intgr      — SST stop-signal RT
#   nc_y_est_interference_rt               — Stroop interference RT
#   mr_y_tfmri__nback__beh__0b_rt_cleaned  — N-back 0-back RT
#   mr_y_tfmri__nback__beh__2b_rt_cleaned  — N-back 2-back RT
#
# Variables NOT flipped (higher already = better):
#   All NIH Toolbox scores
#   nc_y_est_interference_acc
#   mr_y_tfmri__nback__beh__0b_acc_cleaned
#   mr_y_tfmri__nback__beh__2b_acc_cleaned
#
# Load cost variables have no original units (z-scored difference scores)
# and are handled separately in the plotting script.
# ------------------------------------------------------------------------------
ds6_plotting_ref <- ds6_clean_sst %>%
  mutate(
    mr_y_tfmri__sst__beh__ssrt_intgr = mr_y_tfmri__sst__beh__ssrt_intgr * -1,
    nc_y_est_interference_rt = nc_y_est_interference_rt * -1,
    mr_y_tfmri__nback__beh__0b_rt_cleaned = mr_y_tfmri__nback__beh__0b_rt_cleaned * -1,
    mr_y_tfmri__nback__beh__2b_rt_cleaned = mr_y_tfmri__nback__beh__2b_rt_cleaned * -1,
    bilingual_use = (fc_y_aclt_002__03_complete + fc_y_aclt_002__04_complete) - 2
  )

cat("ds6_plotting_ref finalized with flipped RT variables.\n")
cat("Rows:", nrow(ds6_plotting_ref), "\n")

# ------------------------------------------------------------------------------
# Z-SCORING
# ------------------------------------------------------------------------------
vars_to_scale_sst <- c(
  "ab_p_demo_age",
  "nc_y_nihtb__lswmt__uncor_score",
  "nc_y_nihtb__readr__uncor_score",
  "nc_y_nihtb__picvcb__uncor_score",
  "nc_y_nihtb__comp__cryst__uncor_score",
  "nc_y_nihtb__comp__fluid__uncor_score",
  "nc_y_nihtb__crdst__uncor_score",
  "nc_y_nihtb__flnkr__uncor_score",
  "nc_y_nihtb__pttcp__uncor_score",
  "nc_y_est_acc",
  "nc_y_est_rt",
  "nc_y_est_interference_rt",
  "nc_y_est_interference_acc",
  "nc_y_nihtb__readr__uncor_score_complete",
  "nc_y_nihtb__comp__cryst__uncor_score_complete",
  "nc_y_nihtb__comp__fluid__uncor_score_complete",
  "mr_y_tfmri__sst__beh__ssrt_intgr_cleaned",
  "mr_y_tfmri__sst__beh__ssrt_mean_cleaned",
  paste0("gn_y_popstruct_pc__", sprintf("%02d", 1:10), "_complete")
)
ds6_clean_sst[vars_to_scale_sst] <- scale(ds6_clean_sst[vars_to_scale_sst])
ds6_clean_sst <- ds6_clean_sst %>%
  rename_with(~ paste0(., "_z"), all_of(vars_to_scale_sst))

cat(
  "RT Interference Range (Z): ",
  range(ds6_clean_sst$nc_y_est_interference_rt_z, na.rm = TRUE), "\n"
)
cat(
  "SSRT Intgr Range (Z): ",
  range(ds6_clean_sst$mr_y_tfmri__sst__beh__ssrt_intgr_cleaned_z, na.rm = TRUE), "\n"
)

# ------------------------------------------------------------------------------
# BUILD ds6_clean_combined: add N-back z-scored variables for Run 1
# N-back load cost (raw) is also retained here for descriptive reporting.
# Sign convention (positive = worse performance under load):
#   Accuracy cost: 0-back minus 2-back (higher = bigger accuracy drop)
#   RT cost:       2-back minus 0-back (higher = bigger slowing)
# ------------------------------------------------------------------------------
vars_to_scale_nback <- c(
  "mr_y_tfmri__nback__beh__0b_acc_cleaned",
  "mr_y_tfmri__nback__beh__0b_rt_cleaned",
  "mr_y_tfmri__nback__beh__2b_acc_cleaned",
  "mr_y_tfmri__nback__beh__2b_rt_cleaned"
)
ds6_clean_combined <- ds6_clean_sst
ds6_clean_combined[vars_to_scale_nback] <- scale(ds6_clean_combined[vars_to_scale_nback])
ds6_clean_combined <- ds6_clean_combined %>%
  rename_with(~ paste0(., "_z"), all_of(vars_to_scale_nback)) %>%
  mutate(
    mr_y_tfmri__nback__beh__load_acc =
      mr_y_tfmri__nback__beh__0b_acc_cleaned_z - mr_y_tfmri__nback__beh__2b_acc_cleaned_z,
    mr_y_tfmri__nback__beh__load_rt =
      mr_y_tfmri__nback__beh__2b_rt_cleaned_z - mr_y_tfmri__nback__beh__0b_rt_cleaned_z
  )

cat(
  "N-Back 0b Acc Range (Z): ",
  range(ds6_clean_combined$mr_y_tfmri__nback__beh__0b_acc_cleaned_z, na.rm = TRUE), "\n"
)
cat(
  "N-Back 2b RT Range (Z): ",
  range(ds6_clean_combined$mr_y_tfmri__nback__beh__2b_rt_cleaned_z, na.rm = TRUE), "\n"
)

# ==============================================================================
# PHASE 4: EXPORT TO BLIMP
# ==============================================================================
setwd(dir_proj)

# --- Run 1: Main (1,3,5,7) — uses ds6_clean_combined (includes N-back) ---
ds6_balanced_blimp_1357 <- ds6_clean_combined %>%
  filter(session_id %in% c(1, 3, 5, 7)) %>%
  select(-mr_y_tfmri__nback__beh__load_acc, -mr_y_tfmri__nback__beh__load_rt) %>%
  arrange(participant_id, session_id)
# write.csv(ds6_balanced_blimp_1357, "abcd6_0_bilingual_balanced_impute_1357.csv", row.names=FALSE, na="NA")

# --- Runs 2, 3, 4 use ds6_clean_sst (original pipeline, no N-back) ---
# This exactly matches the column structure of the BLIMP files for these runs.

# --- Run 2: Stroop (2,4,6) ---
# Column order matches the original BLIMP export exactly.
# Note: no N-back variables, no ntvlang_003, no design_mr vars after fam_complete.
stroop_cols <- c(
  "participant_id", "session_id",
  "nc_y_nihtb__lswmt__uncor_score_z", "nc_y_nihtb__readr__uncor_score_z",
  "nc_y_nihtb__picvcb__uncor_score_z", "nc_y_nihtb__comp__cryst__uncor_score_z",
  "nc_y_nihtb__comp__fluid__uncor_score_z", "nc_y_nihtb__crdst__uncor_score_z",
  "nc_y_nihtb__flnkr__uncor_score_z", "nc_y_nihtb__pttcp__uncor_score_z",
  "nc_y_est_acc_z", "nc_y_est_rt_z",
  "nc_y_est__congr_acc", "nc_y_est__congr_rt",
  "nc_y_est__incongr_acc", "nc_y_est__incongr_rt",
  "fc_y_aclt_001", "fc_y_aclt_002", "fc_y_aclt_002__01",
  "fc_y_aclt_002__02", "fc_y_aclt_002__03", "fc_y_aclt_002__04",
  "fc_p_aclt_001", "fc_p_aclt_002", "fc_p_aclt_002__01",
  "fc_p_aclt_002__02", "fc_p_aclt_002__03", "fc_p_aclt_002__04",
  "mr_y_tfmri__sst__beh__ssrt_intgr", "mr_y_tfmri__sst__beh__ssrt_mean",
  "mr_y_tfmri__sst__beh__violat_indicator", "mr_y_tfmri__sst__beh__qc_indicator",
  "mr_y_tfmri__sst__beh__coderr_indicator",
  "ab_p_demo_age_z", "ab_g_dyn__cohort_edu__cgs",
  "ab_g_dyn__cohort_income__hhold__3lvl", "ab_g_dyn__design_site",
  "ab_g_dyn__design_mr__serial", "ab_g_dyn__design_mr__manufact",
  "ab_g_dyn__design_mr__software",
  "fc_y_aclt_002_complete", "fc_y_aclt_002__02_complete",
  "fc_y_aclt_002__03_complete", "fc_y_aclt_002__04_complete",
  "ab_p_demo__ntvlang_003_complete", "ab_g_stc__cohort_sex_complete",
  "ab_g_dyn__cohort_edu__cgs_complete", "ab_g_dyn__cohort_income__hhold__3lvl_complete",
  "ab_g_stc__cohort_ethn_complete",
  "nc_y_nihtb__readr__uncor_score_complete_z",
  "nc_y_nihtb__comp__cryst__uncor_score_complete_z",
  "nc_y_nihtb__comp__fluid__uncor_score_complete_z",
  "gn_y_popstruct_pc__01_complete_z", "gn_y_popstruct_pc__02_complete_z",
  "gn_y_popstruct_pc__03_complete_z", "gn_y_popstruct_pc__04_complete_z",
  "gn_y_popstruct_pc__05_complete_z", "gn_y_popstruct_pc__06_complete_z",
  "gn_y_popstruct_pc__07_complete_z", "gn_y_popstruct_pc__08_complete_z",
  "gn_y_popstruct_pc__09_complete_z", "gn_y_popstruct_pc__10_complete_z",
  "ab_g_stc__design_id__fam_complete",
  "mr_y_tfmri__sst__beh__ssrt_intgr_cleaned_z",
  "mr_y_tfmri__sst__beh__ssrt_mean_cleaned_z",
  "nc_y_est_interference_rt_z",
  "nc_y_est_interference_acc_z"
)
ds6_balanced_blimp_246_stroop <- ds6_clean_sst %>%
  filter(session_id %in% c(2, 4, 6)) %>%
  select(all_of(stroop_cols)) %>%
  arrange(participant_id, session_id)
# write.csv(ds6_balanced_blimp_246_stroop, "abcd6_0_bilingual_balanced_impute_246_stroop_CORRECTED.csv", row.names=FALSE, na="NA")

# --- Run 3: Card Sort (1,7) ---
# Derived from stroop column set minus the bilingual acculturation items.
vars_to_drop_cardsort <- c(
  "fc_y_aclt_001", "fc_y_aclt_002", "fc_y_aclt_002__01",
  "fc_y_aclt_002__02", "fc_y_aclt_002__03", "fc_y_aclt_002__04",
  "fc_p_aclt_001", "fc_p_aclt_002", "fc_p_aclt_002__01",
  "fc_p_aclt_002__02", "fc_p_aclt_002__03", "fc_p_aclt_002__04"
)
ds6_balanced_blimp_17_cardsort <- ds6_balanced_blimp_246_stroop %>%
  filter(session_id %in% c(1, 7)) %>%
  select(-any_of(vars_to_drop_cardsort))
# write.csv(ds6_balanced_blimp_17_cardsort, "abcd6_0_bilingual_balanced_impute_17_cardsort.csv", row.names=FALSE, na="NA")

# --- Run 4: List Sort (1,5,7) ---
ds6_balanced_blimp_157_list <- ds6_clean_sst %>%
  filter(session_id %in% c(1, 5, 7)) %>%
  arrange(participant_id, session_id) %>%
  select(
    participant_id, session_id,
    ab_g_stc__cohort_sex_complete, ab_g_dyn__design_site, ab_g_stc__design_id__fam_complete,
    ab_p_demo_age_z, ab_g_dyn__cohort_edu__cgs_complete, ab_g_dyn__cohort_income__hhold__3lvl_complete,
    fc_y_aclt_002__02_complete, fc_y_aclt_002__04_complete, fc_y_aclt_002__03_complete,
    starts_with("gn_y_popstruct_pc"),
    nc_y_nihtb__lswmt__uncor_score_z, nc_y_nihtb__readr__uncor_score_complete_z,
    nc_y_nihtb__picvcb__uncor_score_z, nc_y_nihtb__flnkr__uncor_score_z,
    mr_y_tfmri__sst__beh__ssrt_mean_cleaned_z
  )
# write.csv(ds6_balanced_blimp_157_list, "abcd6_0_bilingual_balanced_impute_157_list.csv", row.names=FALSE, na="NA")

# --- bilingual_use composite (pre-imputation) ---
# Mirrors the post-imputation calc_bilingual_use applied to imp_list_* objects.
# fc_y_aclt_002__03_complete = language use with friends (1-5)
# fc_y_aclt_002__04_complete = language use with family  (1-5)
# Raw sum range 2-10; shift by -2 -> 0 (bilingual) to 8 (monolingual)
calc_bilingual_use <- function(df) {
  df %>%
    mutate(
      bilingual_use = (fc_y_aclt_002__03_complete + fc_y_aclt_002__04_complete) - 2
    )
}

ds6_balanced_blimp_1357        <- calc_bilingual_use(ds6_balanced_blimp_1357)
ds6_balanced_blimp_246_stroop  <- calc_bilingual_use(ds6_balanced_blimp_246_stroop)
ds6_balanced_blimp_17_cardsort <- calc_bilingual_use(ds6_balanced_blimp_17_cardsort)
ds6_balanced_blimp_157_list    <- calc_bilingual_use(ds6_balanced_blimp_157_list)

cat("bilingual_use calculated for all pre-imputation data frames.\n")
cat("------------------------------------------------------------------\n")
cat("EXPORTS COMPLETE.\n")
cat("PLEASE STOP HERE AND RUN BLIMP EXTERNALLY BEFORE PROCEEDING.\n")
cat("------------------------------------------------------------------\n")
# *** SAFETY STOP ***
stop("Script paused: Run external Blimp imputations now.")


# ==============================================================================
# PHASE 5: IMPORT BLIMP OUTPUT
# ==============================================================================
# Reads BLIMP stacked output and assigns full column names positionally from
# the original Phase 4 export dataframe (BLIMP truncates long names).
# The export dataframes in Phase 4 are constructed to exactly match the column
# sets sent to BLIMP, so the count check here should always pass.
process_blimp_output <- function(imputed_file, original_df) {
  imp_data <- read_csv(imputed_file, col_names = TRUE, show_col_types = FALSE)
  orig_names <- names(original_df)
  n_data <- ncol(imp_data) - 1L
  if (n_data != length(orig_names)) {
    stop(paste0(
      "Mismatch in ", imputed_file, "!\n",
      "  Imputed file data columns: ", n_data, "\n",
      "  Original df columns:       ", length(orig_names), "\n",
      "  The Phase 4 export dataframe does not match the BLIMP input file.\n",
      "  Check Phase 4 export construction."
    ))
  }
  colnames(imp_data) <- c("imputation", orig_names)
  imp_data
}

# 1. Load imputed data (filenames must match BLIMP output)
imp_vocab_sst_1357 <- process_blimp_output("imputed_vocab_sstintgr_flanker_nback_1357_stacked.csv", ds6_balanced_blimp_1357)
imp_stroop_246 <- process_blimp_output("imputed_246_stroop_stacked.csv", ds6_balanced_blimp_246_stroop)
imp_cardsort_17 <- process_blimp_output("imputed_17_cardsort_stacked.csv", ds6_balanced_blimp_17_cardsort)
imp_listsort_157 <- process_blimp_output("imputed_157_listsort_stacked.csv", ds6_balanced_blimp_157_list)

# 2. Convert to MITML lists
imp_list_vocab_sst_1357 <- as.mitml.list(split(imp_vocab_sst_1357, imp_vocab_sst_1357$imputation))
imp_list_stroop_246 <- as.mitml.list(split(imp_stroop_246, imp_stroop_246$imputation))
imp_list_cardsort_17 <- as.mitml.list(split(imp_cardsort_17, imp_cardsort_17$imputation))
imp_list_listsort_157 <- as.mitml.list(split(imp_listsort_157, imp_listsort_157$imputation))

# ==============================================================================
# PHASE 5b: POST-IMPUTATION SCORING
# ==============================================================================
# Scoring convention: HIGHER = BETTER PERFORMANCE on all outcome variables.
# This ensures a positive bilingual_use coefficient uniformly means monolinguals
# outperform bilinguals, and a negative coefficient means bilinguals outperform.
#
# Variables already higher-is-better (no transformation):
#   Vocabulary, Flanker, Card Sort, List Sort, N-Back accuracy
#
# Variables requiring transformation (higher raw = worse):
#   SSRT:               higher = slower stopping        -> multiply by -1
#   N-Back RT (0b, 2b): higher = slower                 -> multiply by -1
#   N-Back load RT:     computed as (0b_rt - 2b_rt), z-scored  -> less slowing = higher
#   N-Back load Acc:    computed as (2b_acc - 0b_acc), z-scored -> less drop = higher
#   Stroop int. RT:     higher = more cost               -> multiply by -1
#   Stroop int. Acc:    higher = more cost               -> multiply by -1

# --- Step 1: bilingual_use predictor ---
# Runs 1, 2, 4 contain acculturation items; Run 3 (cardsort) does not.
# bilingual_use for cardsort is joined from the Run 2 stroop list (same
# participants, overlapping sessions 1 & 7).
calc_bilingual_use <- function(df) {
  df %>%
    mutate(
      # fc_y_aclt_002__03_complete = language use with friends (1-5)
      # fc_y_aclt_002__04_complete = language use with family  (1-5)
      # Raw sum range 2-10; shift by -2 -> 0 (bilingual) to 8 (monolingual)
      bilingual_use = (fc_y_aclt_002__03_complete + fc_y_aclt_002__04_complete) - 2
    )
}

imp_list_vocab_sst_1357 <- lapply(imp_list_vocab_sst_1357, calc_bilingual_use)
imp_list_stroop_246 <- lapply(imp_list_stroop_246, calc_bilingual_use)
imp_list_cardsort_17 <- lapply(imp_list_cardsort_17, calc_bilingual_use)
imp_list_listsort_157 <- lapply(imp_list_listsort_157, calc_bilingual_use)

cat("bilingual_use calculated for all lists.\n")

# --- Step 2: score Run 1 outcomes (vocab / SST / flanker / N-back) ---
score_run1 <- function(df) {
  df %>%
    mutate(
      # SST: flip so higher = faster stopping (better)
      sst_ssrt_z = mr_y_tfmri__sst__beh__ssrt_intgr_cleaned_z * -1,

      # N-Back RT: flip so higher = faster (better)
      nback_0b_rt_z = mr_y_tfmri__nback__beh__0b_rt_cleaned_z * -1,
      nback_2b_rt_z = mr_y_tfmri__nback__beh__2b_rt_cleaned_z * -1,

      # N-Back accuracy: higher already = better; rename for consistency
      nback_0b_acc_z = mr_y_tfmri__nback__beh__0b_acc_cleaned_z,
      nback_2b_acc_z = mr_y_tfmri__nback__beh__2b_acc_cleaned_z,

      # N-Back load cost RT: (0b_rt - 2b_rt) -> positive = less slowing (better)
      # z-score within each imputed dataset (correct MI procedure)
      nback_load_rt_z = as.numeric(scale(
        mr_y_tfmri__nback__beh__0b_rt_cleaned_z -
          mr_y_tfmri__nback__beh__2b_rt_cleaned_z
      )),

      # N-Back load cost Acc: (2b_acc - 0b_acc) -> positive = less drop (better)
      # z-score within each imputed dataset
      nback_load_acc_z = as.numeric(scale(
        mr_y_tfmri__nback__beh__2b_acc_cleaned_z -
          mr_y_tfmri__nback__beh__0b_acc_cleaned_z
      ))
    )
}

imp_list_vocab_sst_1357 <- lapply(imp_list_vocab_sst_1357, score_run1)
class(imp_list_vocab_sst_1357) <- c("mitml.list", "list")
cat("Run 1 outcomes scored.\n")

# --- Step 3: score Run 2 stroop interference outcomes ---
# interference_rt_z  = (incongr_rt - congr_rt)  -> higher = more cost -> flip
# interference_acc_z = (congr_acc  - incongr_acc) -> higher = more cost -> flip, then winsorize
score_run2 <- function(df) {
  df %>%
    mutate(
      stroop_int_rt_z  = nc_y_est_interference_rt_z * -1,
      # Flip then winsorize at +-3 SD. Two problems justify winsorizing:
      # (1) Ceiling effects on raw accuracy (median congr=0.98, incongr=0.94)
      #     leave very little variance in the difference score, so modest
      #     deviations produce large z-scores.
      # (2) Participants with fully missing stroop sessions receive BLIMP-imputed
      #     interference scores reaching +-15 SD. Winsorizing retains ordinal
      #     information while preventing these cases from unduly influencing models.
      stroop_int_acc_z = pmin(pmax(nc_y_est_interference_acc_z * -1, -3), 3)
    )
}

imp_list_stroop_246 <- lapply(imp_list_stroop_246, score_run2)
class(imp_list_stroop_246) <- c("mitml.list", "list")
cat("Run 2 (stroop) outcomes scored.\n")

# --- Step 4: factor conversion (all lists) ---
convert_factors <- function(df) {
  df %>%
    mutate(
      ab_g_stc__cohort_sex_complete                 = factor(ab_g_stc__cohort_sex_complete),
      ab_g_dyn__cohort_edu__cgs_complete            = as.ordered(ab_g_dyn__cohort_edu__cgs_complete),
      ab_g_dyn__cohort_income__hhold__3lvl_complete = as.ordered(ab_g_dyn__cohort_income__hhold__3lvl_complete)
    )
}

imp_list_vocab_sst_1357 <- lapply(imp_list_vocab_sst_1357, convert_factors)
class(imp_list_vocab_sst_1357) <- c("mitml.list", "list")
imp_list_stroop_246 <- lapply(imp_list_stroop_246, convert_factors)
class(imp_list_stroop_246) <- c("mitml.list", "list")
imp_list_cardsort_17 <- lapply(imp_list_cardsort_17, convert_factors)
class(imp_list_cardsort_17) <- c("mitml.list", "list")
imp_list_listsort_157 <- lapply(imp_list_listsort_157, convert_factors)
class(imp_list_listsort_157) <- c("mitml.list", "list")
cat("Factor conversion complete.\n")

# ==============================================================================
# PHASE 5c: MODELING
# ==============================================================================
# Helper: build lmer formula string
build_formula <- function(outcome, preds, covars, rand) {
  paste(outcome, "~", preds, "+", paste(covars, collapse = " + "), "+", rand)
}

# Helper: print pooled estimates + 95% CIs
summarize_with_ci <- function(res_obj, model_name) {
  cat("\n=======================================================\n")
  cat("MODEL:", model_name, "\n")
  cat("=======================================================\n")
  print(res_obj)
  cat("\n--- 95% Confidence Intervals ---\n")
  print(confint(res_obj))
  cat("\n")
}

# Plotting function — higher y always = better performance
plot_model_interaction <- function(pooled_result, title_str, y_label) {
  coefs <- pooled_result$estimates
  b0 <- coefs["(Intercept)", "Estimate"]
  b_age <- coefs["ab_p_demo_age_z", "Estimate"]
  b_bil <- coefs["bilingual_use", "Estimate"]
  b_int <- coefs["bilingual_use:ab_p_demo_age_z", "Estimate"]

  pred_data <- expand.grid(
    age_z         = seq(-2.5, 2.5, length.out = 100),
    bilingual_use = 0:8
  )
  pred_data$predicted_y <- b0 +
    b_age * pred_data$age_z +
    b_bil * pred_data$bilingual_use +
    b_int * pred_data$age_z * pred_data$bilingual_use

  ggplot(pred_data, aes(
    x = age_z, y = predicted_y,
    group = factor(bilingual_use),
    color = factor(bilingual_use)
  )) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    scale_color_viridis_d(
      option = "viridis", direction = 1,
      name = "Language Use\n(0=Bilingual\n8=Monolingual)"
    ) +
    theme_minimal() +
    labs(
      title = title_str,
      subtitle = "Higher = Better Performance",
      x = "Age (Z-Score)",
      y = y_label
    ) +
    theme(legend.position = "right")
}

# --- Random effects structures ---
# Random slope: for outcomes where age trajectory varies by participant
random_effects_slope <- paste(
  "(1 | ab_g_dyn__design_site)",
  "(1 | ab_g_stc__design_id__fam_complete)",
  "(1 + ab_p_demo_age_z | participant_id)",
  sep = " + "
)
# Random intercept: for all other outcomes (avoids singular fits)
random_effects_intercept <- paste(
  "(1 | ab_g_dyn__design_site)",
  "(1 | ab_g_stc__design_id__fam_complete)",
  "(1 | participant_id)",
  sep = " + "
)

# --- Universal covariates ---
covariates_universal <- c(
  "ab_g_stc__cohort_sex_complete",
  "ab_g_dyn__cohort_edu__cgs_complete",
  "ab_g_dyn__cohort_income__hhold__3lvl_complete",
  "nc_y_nihtb__readr__uncor_score_complete_z",
  paste0("gn_y_popstruct_pc__", sprintf("%02d", 1:10), "_complete_z")
)

# --- Optimizer ---
strict_control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ==============================================================================
# DIAGNOSTIC: OUTCOME DISTRIBUTIONS (single randomly sampled imputation)
# ==============================================================================
# Samples one imputed dataset at random from each list, extracts all scored
# outcome variables, and plots histograms in a single combined figure.
# Purpose: visual check that distributions look reasonable post-scoring.

samp_run1 <- imp_list_vocab_sst_1357[[sample(length(imp_list_vocab_sst_1357), 1)]]
samp_run2 <- imp_list_stroop_246[[sample(length(imp_list_stroop_246), 1)]]
samp_run3 <- imp_list_cardsort_17[[sample(length(imp_list_cardsort_17), 1)]]
samp_run4 <- imp_list_listsort_157[[sample(length(imp_list_listsort_157), 1)]]

outcome_samples <- bind_rows(
  tibble(outcome = "Vocabulary", value = samp_run1$nc_y_nihtb__picvcb__uncor_score_z),
  tibble(outcome = "List Sort WM", value = samp_run4$nc_y_nihtb__lswmt__uncor_score_z),
  tibble(outcome = "SST (flipped)", value = samp_run1$sst_ssrt_z),
  tibble(outcome = "Flanker", value = samp_run1$nc_y_nihtb__flnkr__uncor_score_z),
  tibble(outcome = "Card Sort", value = samp_run3$nc_y_nihtb__crdst__uncor_score_z),
  tibble(outcome = "Stroop Int. RT\n(flipped)", value = samp_run2$stroop_int_rt_z),
  tibble(outcome = "Stroop Int. Acc\n(flipped)", value = samp_run2$stroop_int_acc_z),
  tibble(outcome = "N-Back 0b Acc", value = samp_run1$nback_0b_acc_z),
  tibble(outcome = "N-Back 0b RT\n(flipped)", value = samp_run1$nback_0b_rt_z),
  tibble(outcome = "N-Back 2b Acc", value = samp_run1$nback_2b_acc_z),
  tibble(outcome = "N-Back 2b RT\n(flipped)", value = samp_run1$nback_2b_rt_z),
  tibble(outcome = "N-Back Load RT\n(flipped, z)", value = samp_run1$nback_load_rt_z),
  tibble(outcome = "N-Back Load Acc\n(z)", value = samp_run1$nback_load_acc_z)
) %>%
  filter(!is.na(value)) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome)))

p_histograms <- ggplot(outcome_samples, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#4682B4", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  facet_wrap(~outcome, scales = "free", ncol = 4) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "Outcome Distributions (single sampled imputation, all outcomes higher = better)",
    subtitle = "Dashed line = 0. Distributions should be roughly unimodal and near-normal.",
    x        = "Score (Z)",
    y        = "Count"
  ) +
  theme(
    strip.text    = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )

print(p_histograms)

# ==============================================================================
# MODELS
# All outcomes: higher = better.
# Positive bilingual_use coefficient -> monolinguals outperform bilinguals.
# Negative bilingual_use coefficient -> bilinguals outperform monolinguals.
# ==============================================================================

# --- 1. Vocabulary ---
# Higher z-score = larger vocabulary. No transformation needed.
message("Running Vocab Model...")
f_vocab <- build_formula(
  "nc_y_nihtb__picvcb__uncor_score_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_slope
)
fit_vocab <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_vocab), REML = TRUE, control = strict_control)
)
res_vocab <- testEstimates(fit_vocab, var.comp = TRUE)

# --- 2. Working Memory (List Sort) ---
# Higher z-score = better working memory. No transformation needed.
message("Running List Sort Model...")
f_listsort <- build_formula(
  "nc_y_nihtb__lswmt__uncor_score_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_listsort <- with(
  imp_list_listsort_157,
  lmer(as.formula(f_listsort), REML = TRUE, control = strict_control)
)
res_listsort <- testEstimates(fit_listsort, var.comp = TRUE)

# --- 3. Inhibitory Control (SST) ---
# sst_ssrt_z = ssrt_intgr_z * -1: higher = faster stopping (better).
message("Running SST Model...")
f_sst <- build_formula(
  "sst_ssrt_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_sst <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_sst), REML = TRUE, control = strict_control)
)
res_sst <- testEstimates(fit_sst, var.comp = TRUE)

# --- 4. Attention (Flanker) ---
# Higher z-score = better attention. No transformation needed.
message("Running Flanker Model...")
f_flanker <- build_formula(
  "nc_y_nihtb__flnkr__uncor_score_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_slope
)
fit_flanker <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_flanker), REML = TRUE, control = strict_control)
)
res_flanker <- testEstimates(fit_flanker, var.comp = TRUE)

# --- 5. Cognitive Flexibility (Card Sort) ---
# Higher z-score = better flexibility. No transformation needed.
message("Running Card Sort Model...")
f_cardsort <- build_formula(
  "nc_y_nihtb__crdst__uncor_score_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_cardsort <- with(
  imp_list_cardsort_17,
  lmer(as.formula(f_cardsort), REML = TRUE, control = strict_control)
)
res_cardsort <- testEstimates(fit_cardsort, var.comp = TRUE)

# --- 6. Stroop Interference RT ---
# stroop_int_rt_z = interference_rt_z * -1: higher = less RT cost (better).
message("Running Stroop Interference RT Model...")
f_stroop_int_rt <- build_formula(
  "stroop_int_rt_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_stroop_int_rt <- with(
  imp_list_stroop_246,
  lmer(as.formula(f_stroop_int_rt), REML = TRUE, control = strict_control)
)
res_stroop_int_rt <- testEstimates(fit_stroop_int_rt, var.comp = TRUE)

# --- 7. Stroop Interference Accuracy ---
# stroop_int_acc_z = interference_acc_z * -1: higher = less accuracy cost (better).
# Simplified RE (participant intercept only) to avoid singular fit.
message("Running Stroop Interference Accuracy Model...")
f_stroop_int_acc <- build_formula(
  "stroop_int_acc_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, "(1 | participant_id)"
)
fit_stroop_int_acc <- with(
  imp_list_stroop_246,
  lmer(as.formula(f_stroop_int_acc), REML = TRUE, control = strict_control)
)
res_stroop_int_acc <- testEstimates(fit_stroop_int_acc, var.comp = TRUE)

# --- 8. N-Back 0-Back Accuracy ---
# Higher = more accurate at baseline attention load. No transformation needed.
message("Running N-Back 0-Back Accuracy Model...")
f_nback_0b_acc <- build_formula(
  "nback_0b_acc_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_0b_acc <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_0b_acc), REML = TRUE, control = strict_control)
)
res_nback_0b_acc <- testEstimates(fit_nback_0b_acc, var.comp = TRUE)

# --- 9. N-Back 0-Back RT ---
# nback_0b_rt_z = 0b_rt_z * -1: higher = faster at baseline (better).
message("Running N-Back 0-Back RT Model...")
f_nback_0b_rt <- build_formula(
  "nback_0b_rt_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_0b_rt <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_0b_rt), REML = TRUE, control = strict_control)
)
res_nback_0b_rt <- testEstimates(fit_nback_0b_rt, var.comp = TRUE)

# --- 10. N-Back 2-Back Accuracy ---
# Higher = more accurate under working memory load. No transformation needed.
message("Running N-Back 2-Back Accuracy Model...")
f_nback_2b_acc <- build_formula(
  "nback_2b_acc_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_2b_acc <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_2b_acc), REML = TRUE, control = strict_control)
)
res_nback_2b_acc <- testEstimates(fit_nback_2b_acc, var.comp = TRUE)

# --- 11. N-Back 2-Back RT ---
# nback_2b_rt_z = 2b_rt_z * -1: higher = faster under working memory load (better).
message("Running N-Back 2-Back RT Model...")
f_nback_2b_rt <- build_formula(
  "nback_2b_rt_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_2b_rt <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_2b_rt), REML = TRUE, control = strict_control)
)
res_nback_2b_rt <- testEstimates(fit_nback_2b_rt, var.comp = TRUE)

# --- 12. N-Back Load Cost RT ---
# nback_load_rt_z = z-score(0b_rt - 2b_rt): higher = less RT slowing under load (better).
message("Running N-Back Load Cost RT Model...")
f_nback_load_rt <- build_formula(
  "nback_load_rt_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_load_rt <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_load_rt), REML = TRUE, control = strict_control)
)
res_nback_load_rt <- testEstimates(fit_nback_load_rt, var.comp = TRUE)

# --- 13. N-Back Load Cost Accuracy ---
# nback_load_acc_z = z-score(2b_acc - 0b_acc): higher = less accuracy drop under load (better).
message("Running N-Back Load Cost Accuracy Model...")
f_nback_load_acc <- build_formula(
  "nback_load_acc_z",
  "bilingual_use * ab_p_demo_age_z",
  covariates_universal, random_effects_intercept
)
fit_nback_load_acc <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_nback_load_acc), REML = TRUE, control = strict_control)
)
res_nback_load_acc <- testEstimates(fit_nback_load_acc, var.comp = TRUE)

# ==============================================================================
# SUMMARY OUTPUTS
# All outcomes: higher = better. Positive bilingual_use = monolingual advantage.
# ==============================================================================
cat("\n--- ESTABLISHED COGNITIVE TASKS ---\n")
summarize_with_ci(res_vocab, "Vocabulary          [higher = larger vocabulary]")
summarize_with_ci(res_listsort, "List Sort WM        [higher = better working memory]")
summarize_with_ci(res_sst, "SST Inhibition      [higher = faster stopping]")
summarize_with_ci(res_flanker, "Flanker Attention   [higher = better attention]")
summarize_with_ci(res_cardsort, "Card Sort           [higher = better flexibility]")

cat("\n--- STROOP: INHIBITORY CONTROL ---\n")
summarize_with_ci(res_stroop_int_rt, "Stroop Int. RT  [higher = less RT cost]")
summarize_with_ci(res_stroop_int_acc, "Stroop Int. Acc [higher = less accuracy cost]")

cat("\n--- N-BACK: BASELINE (0-BACK) ---\n")
summarize_with_ci(res_nback_0b_acc, "N-Back 0b Accuracy [higher = more accurate]")
summarize_with_ci(res_nback_0b_rt, "N-Back 0b RT       [higher = faster]")

cat("\n--- N-BACK: WORKING MEMORY LOAD (2-BACK) ---\n")
summarize_with_ci(res_nback_2b_acc, "N-Back 2b Accuracy [higher = more accurate under load]")
summarize_with_ci(res_nback_2b_rt, "N-Back 2b RT       [higher = faster under load]")

cat("\n--- N-BACK: LOAD COST ---\n")
summarize_with_ci(res_nback_load_rt, "N-Back Load RT  [higher = less slowing under load]")
summarize_with_ci(res_nback_load_acc, "N-Back Load Acc [higher = less accuracy drop under load]")

# ==============================================================================
# PLOTTING
# All plots: higher y = better. Positive monolingual slope = monolingual advantage.
# ==============================================================================
print(plot_model_interaction(res_vocab, "Vocabulary", "Score (Z) [Higher = Better]"))
print(plot_model_interaction(res_listsort, "Working Memory", "Score (Z) [Higher = Better]"))
print(plot_model_interaction(res_sst, "Inhibitory Control (SST)", "Score (Z) [Higher = Faster Stopping]"))
print(plot_model_interaction(res_flanker, "Attention (Flanker)", "Score (Z) [Higher = Better]"))
print(plot_model_interaction(res_cardsort, "Flexibility (Card Sort)", "Score (Z) [Higher = Better]"))
print(plot_model_interaction(res_stroop_int_rt, "Stroop Interference RT", "Score (Z) [Higher = Less RT Cost]"))
print(plot_model_interaction(res_stroop_int_acc, "Stroop Interference Acc", "Score (Z) [Higher = Less Acc Cost]"))
print(plot_model_interaction(res_nback_0b_acc, "N-Back 0-Back Accuracy", "Score (Z) [Higher = More Accurate]"))
print(plot_model_interaction(res_nback_0b_rt, "N-Back 0-Back RT", "Score (Z) [Higher = Faster]"))
print(plot_model_interaction(res_nback_2b_acc, "N-Back 2-Back Accuracy", "Score (Z) [Higher = More Accurate]"))
print(plot_model_interaction(res_nback_2b_rt, "N-Back 2-Back RT", "Score (Z) [Higher = Faster]"))
print(plot_model_interaction(res_nback_load_rt, "N-Back Load Cost RT", "Score (Z) [Higher = Less Slowing]"))
print(plot_model_interaction(res_nback_load_acc, "N-Back Load Cost Acc", "Score (Z) [Higher = Less Acc Drop]"))


# ==============================================================================
# SUPPLEMENTAL ANALYSIS: RT-ACCURACY COUPLING UNDER WORKING MEMORY LOAD
#
# Two-way: Does RT-accuracy coupling weaken with age under load?
# Three-way: Does the decoupling rate differ by bilingual experience?
#
# Models:
#   fit_sat     — two-way: load cost RT x age -> load cost accuracy
#   fit_sat_bil — three-way: load cost RT x age x bilingual_use -> load cost accuracy
#
# Figures:
#   sat_combined_by_age.tiff   — 3-row panel (0-back / 2-back / load cost)
#   sat_bilingual_by_age.tiff  — load cost coupling by age tertile x language group
# ==============================================================================

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------------

# --- Two-way: Load Cost RT x Age -> Load Cost Accuracy ---
# Key term: nback_load_rt_z:ab_p_demo_age_z
# Negative = RT-accuracy coupling weakens with age (decoupling)
f_sat <- build_formula(
  "nback_load_acc_z",
  "nback_load_rt_z * ab_p_demo_age_z",
  covariates_universal,
  random_effects_intercept
)
fit_sat <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_sat), REML = TRUE, control = strict_control)
)
res_sat <- testEstimates(fit_sat, var.comp = TRUE)
summarize_with_ci(res_sat, "Two-Way: Load Cost RT x Age -> Load Cost Acc")

# --- Three-way: Load Cost RT x Age x Bilingual Use -> Load Cost Accuracy ---
# Key term: nback_load_rt_z:ab_p_demo_age_z:bilingual_use
# Positive = decoupling is weaker as bilingual_use increases toward monolingual
# (i.e., bilinguals decouple faster than monolinguals)
f_sat_bil <- build_formula(
  "nback_load_acc_z",
  "nback_load_rt_z * ab_p_demo_age_z * bilingual_use",
  covariates_universal,
  random_effects_intercept
)
fit_sat_bil <- with(
  imp_list_vocab_sst_1357,
  lmer(as.formula(f_sat_bil), REML = TRUE, control = strict_control)
)
res_sat_bil <- testEstimates(fit_sat_bil, var.comp = TRUE)
summarize_with_ci(res_sat_bil, "Three-Way: Load Cost RT x Age x Bilingual Use -> Load Cost Acc")

# ------------------------------------------------------------------------------
# SHARED SETUP: Sample one imputation and compute age tertiles
# Used by both figures below.
# Note: unflipped raw RT z-scores (mr_y_tfmri__nback__beh__*_rt_cleaned_z) are
# used for the scatter plots so that RT-accuracy correlations are negative and
# interpretable as classical speed-accuracy relationships (slower = less accurate).
# The flipped nback_*_rt_z variables are higher=better and produce artifactually
# positive correlations with accuracy that increase with age, obscuring the
# decoupling signal.
# ------------------------------------------------------------------------------
samp <- imp_list_vocab_sst_1357[[sample(length(imp_list_vocab_sst_1357), 1)]]

samp$age_tertile <- cut(
  samp$ab_p_demo_age_z,
  breaks         = quantile(samp$ab_p_demo_age_z, probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE),
  labels         = c("Younger", "Middle", "Older"),
  include.lowest = TRUE
)

base_fs <- 16

# ------------------------------------------------------------------------------
# HELPER: Build one row of three faceted scatter panels with correlation
# annotations. Uses unflipped RT variables so correlations are negative
# (slower RT -> lower accuracy) and decoupling shows as weakening negative r.
# ------------------------------------------------------------------------------
make_sat_row <- function(dat, x_var, y_var, x_label, y_label, row_title) {
  plot_dat <- dat %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]), !is.na(age_tertile)) %>%
    mutate(
      across(all_of(c(x_var, y_var)), ~ pmin(pmax(.x, -3), 3))
    )

  cors <- plot_dat %>%
    group_by(age_tertile) %>%
    summarise(
      r = cor(.data[[x_var]], .data[[y_var]], use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("r = %.2f", r))

  ggplot(plot_dat, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.08, size = 0.8, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 1.0) +
    geom_text(
      data     = cors,
      aes(label = label),
      x        = Inf, y = Inf,
      hjust    = 1.1, vjust = 1.5,
      size     = base_fs / 3.5,
      fontface = "bold"
    ) +
    facet_wrap(~age_tertile, ncol = 3) +
    labs(title = row_title, x = x_label, y = y_label) +
    theme_minimal(base_size = base_fs) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_fs),
      strip.text = element_text(face = "bold", size = base_fs),
      axis.title = element_text(size = base_fs - 1),
      axis.text  = element_text(size = base_fs - 3)
    )
}

# ------------------------------------------------------------------------------
# FIGURE 1: Three-row panel
# Row 1: 0-back  — baseline, no load; expect stable negative r across age
# Row 2: 2-back  — under load; expect weakening negative r with age
# Row 3: load cost — demand-specific difference; clearest decoupling signal
# All rows use unflipped RT (higher = slower) for consistent interpretation
# ------------------------------------------------------------------------------
p_0b <- make_sat_row(
  dat       = samp,
  x_var     = "mr_y_tfmri__nback__beh__0b_rt_cleaned_z", # unflipped: higher = slower
  y_var     = "nback_0b_acc_z", # higher = more accurate
  x_label   = "0-Back RT (Z) [Higher = Slower]",
  y_label   = "0-Back Acc (Z)\n[Higher = More Accurate]",
  row_title = "Baseline (0-Back): No Load"
)

p_2b <- make_sat_row(
  dat       = samp,
  x_var     = "mr_y_tfmri__nback__beh__2b_rt_cleaned_z", # unflipped: higher = slower
  y_var     = "nback_2b_acc_z", # higher = more accurate
  x_label   = "2-Back RT (Z) [Higher = Slower]",
  y_label   = "2-Back Acc (Z)\n[Higher = More Accurate]",
  row_title = "Under Load (2-Back): Working Memory Demand"
)

p_lc <- make_sat_row(
  dat       = samp,
  x_var     = "nback_load_rt_z", # higher = less slowing (already higher=better)
  y_var     = "nback_load_acc_z", # higher = less accuracy drop
  x_label   = "Load Cost RT (Z) [Higher = Less Slowing]",
  y_label   = "Load Cost Acc (Z)\n[Higher = Less Accuracy Drop]",
  row_title = "Load Cost (2-Back minus 0-Back): Demand-Specific Component"
)

# Combine into single figure
p_combined <- p_0b / p_2b / p_lc +
  plot_annotation(
    title = "RT\u2013Accuracy Coupling by Age and Working Memory Load",
    subtitle = paste0(
      "Coupling weakens with age specifically for load-specific costs,\n",
      "consistent with emerging speed-accuracy flexibility"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = base_fs + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_fs - 1)
    )
  )

print(p_combined)
ggsave("sat_combined_by_age.tiff", p_combined,
  width = 12, height = 14, dpi = 600, compression = "lzw"
)
ggsave("sat_combined_by_age.png", p_combined,
  width = 12, height = 14, dpi = 600
)

# ------------------------------------------------------------------------------
# FIGURE 2: Load cost coupling by age tertile x language group
# Tests whether RT-accuracy decoupling rate differs by bilingual experience.
# Bilingual = bilingual_use 0-2; Monolingual = bilingual_use 6-8.
# Middle range excluded to sharpen the group contrast visually.
# Uses load cost variables (higher=better on both axes) for this figure since
# the group comparison is the focus, not the raw RT-accuracy direction.
# ------------------------------------------------------------------------------
samp_bil <- samp %>%
  mutate(
    lang_group = case_when(
      bilingual_use <= 2 ~ "Bilingual (0-2)",
      bilingual_use >= 6 ~ "Monolingual (6-8)",
      TRUE ~ NA_character_
    ),
    nback_load_rt_z = pmin(pmax(nback_load_rt_z, -3), 3),
    nback_load_acc_z = pmin(pmax(nback_load_acc_z, -3), 3)
  ) %>%
  filter(
    !is.na(lang_group),
    !is.na(nback_load_rt_z),
    !is.na(nback_load_acc_z),
    !is.na(age_tertile)
  )

# Per-group per-tertile correlations
cors_bil <- samp_bil %>%
  group_by(lang_group, age_tertile) %>%
  summarise(
    r = cor(nback_load_rt_z, nback_load_acc_z, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label     = sprintf("r = %.2f", r),
    vjust_val = ifelse(lang_group == "Bilingual (0-2)", 1.5, 3.2)
  )

p_sat_bil <- ggplot(
  samp_bil,
  aes(
    x = nback_load_rt_z,
    y = nback_load_acc_z,
    color = lang_group,
    fill = lang_group
  )
) +
  geom_point(alpha = 0.06, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  geom_text(
    data = cors_bil,
    aes(label = label, vjust = vjust_val),
    x = Inf, y = Inf,
    hjust = 1.1,
    size = base_fs / 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~age_tertile, ncol = 3) +
  scale_color_manual(
    values = c("Bilingual (0-2)" = "#0a2fa8", "Monolingual (6-8)" = "#8b0000"),
    name   = "Language Group"
  ) +
  scale_fill_manual(
    values = c("Bilingual (0-2)" = "#0a2fa8", "Monolingual (6-8)" = "#8b0000"),
    name   = "Language Group"
  ) +
  labs(
    title = "RT\u2013Accuracy Coupling Under Load by Age and Language Group",
    subtitle = expression(paste(
      "Bilinguals show stronger decoupling with age; ",
      "three-way interaction: \u03b2 = 0.008, uncorrected ",
      italic(p), " = .028"
    )),
    x = "Load Cost RT (Z) [Higher = Less Slowing]",
    y = "Load Cost Acc (Z)\n[Higher = Less Accuracy Drop]"
  ) +
  theme_minimal(base_size = base_fs) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = base_fs),
    plot.subtitle   = element_text(hjust = 0.5, size = base_fs - 2),
    strip.text      = element_text(face = "bold", size = base_fs),
    axis.title      = element_text(size = base_fs - 1),
    axis.text       = element_text(size = base_fs - 3),
    legend.position = "bottom",
    legend.text     = element_text(size = base_fs - 2),
    legend.title    = element_text(size = base_fs - 1)
  )

print(p_sat_bil)
ggsave("sat_bilingual_by_age.tiff", p_sat_bil,
  width = 12, height = 6, dpi = 600, compression = "lzw"
)
ggsave("sat_bilingual_by_age.png", p_sat_bil,
  width = 12, height = 6, dpi = 600
)


library(MuMIn)
library(purrr)
library(mitml)
library(xtable)
# 1. HELPER: Pool marginal R2 change across imputations for a given predictor.
get_rsp_r2 <- function(fit_obj, predictor) {
  message(paste("Calculating pooled rsp (R2 method) for:", predictor))
  r2_diffs <- map_dbl(fit_obj, function(m) {
    r2_full <- MuMIn::r.squaredGLMM(m)[1, 1]
    m_reduced <- update(m, as.formula(paste(". ~ . -", predictor)))
    r2_reduced <- MuMIn::r.squaredGLMM(m_reduced)[1, 1]
    max(0, r2_full - r2_reduced)
  })
  sqrt(mean(r2_diffs))
}
# 2. EXTRACTION FUNCTION
extract_stats <- function(res_obj, outcome_name, predictor, rsp_unsigned) {
  estimates_table <- as.data.frame(res_obj$estimates)
  all_preds <- rownames(estimates_table)
  actual_pred <- all_preds[grep(predictor, all_preds)]
  if (length(actual_pred) == 0) {
    return(NULL)
  }
  est_row <- estimates_table[actual_pred[1], ]
  ci_row <- as.data.frame(confint(res_obj))[actual_pred[1], ]
  beta <- as.numeric(est_row[1])
  se <- as.numeric(est_row[2])
  t_val <- as.numeric(est_row[3])
  p_raw <- as.numeric(est_row[5])
  rsp_val <- rsp_unsigned * sign(beta)
  p_fmt <- if (!is.na(p_raw) && p_raw < 0.001) "$<$ .001" else sprintf("%.3f", p_raw)
  data.frame(
    Outcome = outcome_name,
    Beta = beta,
    SE = se,
    t = t_val,
    p_fmt = p_fmt,
    CI_L = as.numeric(ci_row[1]),
    CI_U = as.numeric(ci_row[2]),
    rsp = rsp_val,
    p_raw = p_raw,
    stringsAsFactors = FALSE
  )
}
# 3. EXECUTION BLOCK
{
  # 13 models total (7 original + 6 N-back)
  model_names <- c(
    "NIH Vocabulary",
    "NIH List Sort WM",
    "SSRT",
    "NIH Flanker",
    "NIH Card Sort",
    "Stroop Interference RT",
    "Stroop Interference Accuracy",
    "EN-Back 0b Accuracy",
    "EN-Back 0b RT",
    "EN-Back 2b Accuracy",
    "EN-Back 2b RT",
    "EN-Back Load RT",
    "EN-Back Load Accuracy"
  )
  fit_list <- setNames(
    list(
      fit_vocab, fit_listsort, fit_sst, fit_flanker,
      fit_cardsort, fit_stroop_int_rt, fit_stroop_int_acc,
      fit_nback_0b_acc, fit_nback_0b_rt,
      fit_nback_2b_acc, fit_nback_2b_rt,
      fit_nback_load_rt, fit_nback_load_acc
    ),
    model_names
  )
  res_list <- setNames(
    list(
      res_vocab, res_listsort, res_sst, res_flanker,
      res_cardsort, res_stroop_int_rt, res_stroop_int_acc,
      res_nback_0b_acc, res_nback_0b_rt,
      res_nback_2b_acc, res_nback_2b_rt,
      res_nback_load_rt, res_nback_load_acc
    ),
    model_names
  )
  pred_age <- "ab_p_demo_age_z"
  pred_use <- "bilingual_use$"
  pred_int <- "bilingual_use:ab_p_demo_age_z"
  # A. Compute pooled unsigned rsp for all 39 cells (13 models x 3 predictors)
  message("--- Computing R2-based rsp values (this may take a few minutes) ---")
  rsp_age <- map(fit_list, ~ get_rsp_r2(.x, pred_age))
  rsp_use <- map(fit_list, ~ get_rsp_r2(.x, "bilingual_use"))
  rsp_int <- map(fit_list, ~ get_rsp_r2(.x, pred_int))
  # B. Extract pooled estimates and attach rsp values
  age_raw <- map2_df(res_list, model_names, ~ extract_stats(.x, .y, pred_age, rsp_age[[.y]]))
  use_raw <- map2_df(res_list, model_names, ~ extract_stats(.x, .y, pred_use, rsp_use[[.y]]))
  int_raw <- map2_df(res_list, model_names, ~ extract_stats(.x, .y, pred_int, rsp_int[[.y]]))
  # C. FDR correction across full family of 39 tests (Benjamini-Hochberg)
  all_data <- rbind(age_raw, use_raw, int_raw)
  all_data$q_val <- p.adjust(all_data$p_raw, method = "fdr")
  all_data$q_fmt <- ifelse(
    !is.na(all_data$q_val) & all_data$q_val < 0.001,
    "$<$ .001",
    sprintf("%.3f", all_data$q_val)
  )
  final_table <- all_data[, c("Outcome", "Beta", "SE", "t", "p_fmt", "CI_L", "CI_U", "rsp", "q_fmt")]
  n <- length(model_names) # 13
  age_f <- final_table[1:n, ]
  use_f <- final_table[(n + 1):(2 * n), ]
  int_f <- final_table[(2 * n + 1):(3 * n), ]
  # D. LaTeX output
  digs <- c(0, 0, 3, 3, 3, 0, 3, 3, 3, 0)
  print_section <- function(df, header) {
    cat(paste0("\\midrule\n\\multicolumn{9}{l}{\\textbf{", header, "}} \\\\\n\\midrule\n"))
    print(
      xtable(df, digits = digs),
      include.rownames = FALSE,
      include.colnames = FALSE,
      only.contents = TRUE,
      hline.after = NULL,
      comment = FALSE,
      sanitize.text.function = identity
    )
  }
  cat("\\begin{table*}[t]\n\\centering\n\\begin{minipage}{\\textwidth}\n\\centering\n")
  cat("\\caption{Associations of Outcomes with Age, Bilingual Use, and Their Interaction}\n")
  cat("\\small\n\\begin{tabular}{lcccccccc}\n\\toprule\n")
  cat("Outcome & $\\beta$ & $SE$ & \\textit{t} & \\textit{p} & 95\\% CI Lower & 95\\% CI Upper & $r_{sp}$ & FDR \\textit{q} \\\\\n")
  print_section(age_f, "Main Effect of Age")
  print_section(use_f, "Main Effect of Bilingual Use")
  print_section(int_f, "Interaction (Age $\\times$ Bilingual Use)")
  cat("\\bottomrule\n\\end{tabular}\n")
  cat("\\label{tab:outcomes}\n\\smallskip\n\\footnotesize\n")
  cat("\\noindent \\caption*{Note. Total analytic sample $n = 11{,}868$. Estimates are pooled across 50 imputations via \\textit{mitml}. ")
  cat("NIH Vocabulary = NIH Toolbox Picture Vocabulary; NIH List Sort (WM) = NIH Toolbox List Sorting Working Memory; ")
  cat("NIH Flanker = NIH Toolbox Flanker; NIH Card Sort = NIH Toolbox Dimensional Change Card Sort; ")
  cat("SSRT = Stop-Signal Reaction Time. Stroop Interference RT and Acc = Reaction Time and Accuracy cost scores ")
  cat("(Incongruent $-$ Congruent); Stroop Interference Accuracy winsorized at $\\pm$3 SD prior to modeling ")
  cat("due to ceiling effects and implausible imputed values in the tails. ")
  cat("N-Back outcomes: 0b = 0-back (baseline); 2b = 2-back (working memory load); ")
  cat("Load RT and Accuracy = z-scored difference scores (0b $-$ 2b for RT; 2b $-$ 0b for Accuracy), ")
  cat("where higher scores indicate less performance cost under load. ")
  cat("All analyses controlled for participant, family, and site as random effects, and sex, highest family education, family income, ")
  cat("NIH Toolbox Oral Reading Recognition, and the top 10 genetic principal components. ")
  cat("$\\beta$ = Standardized beta. SE = Standard error of $\\beta$. 95\\% CIs = 95\\% Confidence Intervals of the $\\beta$. ")
  cat("The semipartial correlation ($r_{sp}$) reflects unique variance explained (marginal $R^2$ change upon predictor removal), ")
  cat("pooled across imputations and signed by the direction of $\\beta$. ")
  cat("FDR \\textit{q} values represent False Discovery Rate adjusted $p$-values across all 39 tests using the Benjamini-Hochberg procedure.}\n")
  cat("\\end{minipage}\n\\end{table*}\n")
}

# ==============================================================================
# SAVE: Model objects and table data
# Run this once after Phase 5 completes. Load block below replaces re-running.
# ==============================================================================
save(
  # Fit objects (needed for rsp / R2 calculations in Phase 6 table)
  fit_vocab, fit_listsort, fit_sst, fit_flanker, fit_cardsort,
  fit_stroop_int_rt, fit_stroop_int_acc,
  fit_nback_0b_acc, fit_nback_0b_rt,
  fit_nback_2b_acc, fit_nback_2b_rt,
  fit_nback_load_rt, fit_nback_load_acc,
  # Pooled result objects (needed for estimates, CIs, summaries)
  res_vocab, res_listsort, res_sst, res_flanker, res_cardsort,
  res_stroop_int_rt, res_stroop_int_acc,
  res_nback_0b_acc, res_nback_0b_rt,
  res_nback_2b_acc, res_nback_2b_rt,
  res_nback_load_rt, res_nback_load_acc,
  file = "bilingual_model_objects.RData"
)
cat("Model objects saved to bilingual_model_objects.RData\n")

# ==============================================================================
# SAVE: Final table data (age_raw, use_raw, int_raw, final_table)
# Saved separately so the table can be re-exported without re-running R2 calcs.
# These are produced inside the execution block above; save after it completes.
# ==============================================================================
save(age_raw, use_raw, int_raw, final_table,
  file = "bilingual_table_data.RData"
)
cat("Table data saved to bilingual_table_data.RData\n")

# Load the model objects (fits and pooled results)
# load("bilingual_model_objects.RData")

# Load the final table data
# load("bilingual_table_data.RData")


################################
#                              #
# PHASE 6: EQUIVALENCY TESTING #
#                              #
################################
library(ggplot2)

# 1. Define the mapping for clean labels
outcome_mapping <- c(
  "NIH Vocabulary" = "NIH Toolbox Vocabulary",
  "NIH List Sort WM" = "NIH Toolbox List Sorting",
  "SSRT" = "Stop-Signal RT",
  "NIH Flanker" = "NIH Toolbox Flanker",
  "NIH Card Sort" = "NIH Toolbox Card Sort",
  "Stroop Interference Accuracy" = "Stroop Interference Accuracy",
  "Stroop Interference RT" = "Stroop Interference RT",
  "N-Back 0b Accuracy" = "EN-Back 0-Back Accuracy",
  "N-Back 0b RT" = "EN-Back 0-Back RT",
  "N-Back 2b Accuracy" = "EN-Back 2-Back Accuracy",
  "N-Back 2b RT" = "EN-Back 2-Back RT",
  "N-Back Load Accuracy" = "EN-Back Load Cost Accuracy",
  "N-Back Load RT" = "EN-Back Load Cost RT"
)

# 2. Prepare plot data
plot_data <- use_raw

# Apply the new names
plot_data$Outcome_Clean <- outcome_mapping[plot_data$Outcome]

# Set factor levels based on the order in the table (reversed for coord_flip)
plot_data$Outcome_Clean <- factor(plot_data$Outcome_Clean, levels = rev(unname(outcome_mapping)))

alpha <- 0.05
# TOST 90% CIs
plot_data$CI_Lower <- plot_data$Beta - qnorm(1 - alpha) * plot_data$SE
plot_data$CI_Upper <- plot_data$Beta + qnorm(1 - alpha) * plot_data$SE

# Equivalence decision
plot_data$Decision <- ifelse(
  plot_data$CI_Lower > -0.05 & plot_data$CI_Upper < 0.05,
  "Equivalent",
  "Not Equivalent"
)
plot_data$Decision <- factor(plot_data$Decision, levels = c("Not Equivalent", "Equivalent"))

# Colors
palette_silver <- "#C0C0C0"
palette_bronze <- "#CD7F32"
ravenclaw_blue <- "#222F5B"
plot_colors <- c("Equivalent" = ravenclaw_blue, "Not Equivalent" = palette_bronze)

# 3. Create Plot
equiv_plot <- ggplot(plot_data, aes(x = Outcome_Clean, y = Beta, color = Decision)) +
  geom_pointrange(
    aes(ymin = CI_Lower, ymax = CI_Upper),
    size = 0.5,
    linewidth = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = palette_silver, linewidth = 1.2) +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = palette_bronze, linewidth = 1.2) +
  coord_flip() +
  scale_color_manual(values = plot_colors, name = "Equivalence Decision (\u00B10.05)") +
  scale_y_continuous(
    limits = c(-0.15, 0.15),
    breaks = seq(-0.15, 0.15, 0.05),
    labels = function(x) sprintf("%.2f", x)
  ) +
  labs(
    y        = expression(Standardized ~ beta),
    title    = "Equivalence Test: Main Effect of Bilingual Use",
    subtitle = "90% CI vs. Equivalence Bounds (\u00B10.05 Bronze; \u00B10.10 Silver)"
  ) +
  theme_bw() +
  theme(
    text              = element_text(size = 16),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.title.y      = element_blank(),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5),
    legend.position   = "bottom"
  )

# 4. Save
ggsave("bilingual_use_equiv.png", equiv_plot, width = 9, height = 6, dpi = 600)
print(equiv_plot)


################################
#                              #
# PHASE 7: RUN AS GAMMs        #
#                              #
################################
# ==============================================================================
# GAMM: POOLED RESULTS + LATEX TABLE
#
# Strategy:
#   - Smooth age terms:    mean EDF + mean F + mean p across 50 imputations
#   - Parametric terms:    Rubin's rules for Beta, SE, t, p, CI
#   - r_sp for age:        Linear predictor variance method from LMM fits
#                          part-R2 = var(beta_age * age) / Y_total
#                          per Stoffel et al. (2021) / Nakagawa et al. (2017)
#                          Uses LMM because r.squaredGLMM is unreliable for GAMs
#                          with random slopes (random slope absorbs age variance,
#                          making delta-R2m near zero despite large fixed effect)
#   - r_sp for bilingual + interaction: delta-R2m from GAMM fits
#                          sqrt(mean(delta_R2)) * sign(beta), pooled across imputations
# ==============================================================================
library(mgcv)
library(MuMIn)
library(mitml)
library(lme4)
library(purrr)
library(dplyr)
library(xtable)
library(parallel)
# Use half of available cores to reduce peak forked-memory demand.
ncores <- max(1, floor(detectCores() * 0.5))

# ==============================================================================
# LOADER: Skip fitting and extraction if saved objects exist
# Comment out the stop() line on first run.
# On subsequent runs, load objects and jump to wherever you need to work.
# ==============================================================================
gamm_fits_path <- file.path(dir_proj, "gamm_fit_objects.RData")
gamm_results_path <- file.path(dir_proj, "gamm_results_and_table.RData")

if (file.exists(gamm_fits_path) && file.exists(gamm_results_path)) {
  message("Saved GAMM objects found. Loading instead of refitting...")
  load(gamm_fits_path)
  load(gamm_results_path)
  message("Loaded: gamm_fit_objects.RData")
  message("Loaded: gamm_results_and_table.RData")
  message("Skip to whichever section you need — models and results are ready.")
  # stop("Loader stop: comment out this line to continue past the loader.")
} else {
  message("No saved objects found — running full Phase 7 pipeline.")
}

# ==============================================================================
# SECTION 0: GAMM MODEL FITTING
# ==============================================================================

# --- Formula Components ---
gamm_terms_hybrid <- "s(ab_p_demo_age_z, k=5) + bilingual_use + bilingual_use:ab_p_demo_age_z"
gamm_rand_slope <- paste(
  "s(ab_g_dyn__design_site, bs='re')",
  "s(ab_g_stc__design_id__fam_complete, bs='re')",
  "s(participant_id, bs='re')",
  "s(participant_id, ab_p_demo_age_z, bs='re')",
  sep = " + "
)
gamm_rand_intercept <- paste(
  "s(ab_g_dyn__design_site, bs='re')",
  "s(ab_g_stc__design_id__fam_complete, bs='re')",
  "s(participant_id, bs='re')",
  sep = " + "
)
gamm_covars <- paste(c(
  "ab_g_stc__cohort_sex_complete",
  "ab_g_dyn__cohort_edu__cgs_complete",
  "ab_g_dyn__cohort_income__hhold__3lvl_complete",
  "nc_y_nihtb__readr__uncor_score_complete_z",
  paste0("gn_y_popstruct_pc__", sprintf("%02d", 1:10), "_complete_z")
), collapse = " + ")

# --- Session-Level Descriptives (first imputation) ---
{
  print_session_stats <- function(df_list, outcome_col, label) {
    df <- df_list[[1]]
    stats <- df %>%
      group_by(session_id) %>%
      summarise(
        N    = n(),
        Mean = mean(!!sym(outcome_col), na.rm = TRUE),
        SD   = sd(!!sym(outcome_col), na.rm = TRUE)
      ) %>%
      ungroup()
    cat("\n---", label, "---\n")
    print(as.data.frame(stats))
  }
  cat("AVERAGE SCORES BY SESSION (Z-SCORES)\n")
  cat("====================================\n")
  print_session_stats(imp_list_vocab_sst_1357, "nc_y_nihtb__picvcb__uncor_score_z", "NIH Vocabulary")
  print_session_stats(imp_list_vocab_sst_1357, "nc_y_nihtb__flnkr__uncor_score_z", "NIH Flanker")
  print_session_stats(imp_list_vocab_sst_1357, "mr_y_tfmri__sst__beh__ssrt_intgr_cleaned_z", "SST SSRT")
  print_session_stats(imp_list_listsort_157, "nc_y_nihtb__lswmt__uncor_score_z", "List Sort (WM)")
  print_session_stats(imp_list_cardsort_17, "nc_y_nihtb__crdst__uncor_score_z", "Card Sort")
  print_session_stats(imp_list_stroop_246, "nc_y_est_interference_rt_z", "Stroop RT")
  print_session_stats(imp_list_stroop_246, "nc_y_est_interference_acc_z", "Stroop Accuracy")
  print_session_stats(imp_list_vocab_sst_1357, "nback_0b_acc_z", "N-Back 0-Back Acc")
  print_session_stats(imp_list_vocab_sst_1357, "nback_0b_rt_z", "N-Back 0-Back RT")
  print_session_stats(imp_list_vocab_sst_1357, "nback_2b_acc_z", "N-Back 2-Back Acc")
  print_session_stats(imp_list_vocab_sst_1357, "nback_2b_rt_z", "N-Back 2-Back RT")
  print_session_stats(imp_list_vocab_sst_1357, "nback_load_rt_z", "N-Back Load Cost RT")
  print_session_stats(imp_list_vocab_sst_1357, "nback_load_acc_z", "N-Back Load Cost Acc")
}

# --- Fit Hybrid GAMMs ---

# 1. Vocabulary (random slope)
f_vocab_hyb <- paste(
  "nc_y_nihtb__picvcb__uncor_score_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_slope
)
message("Fitting Hybrid GAMM: Vocabulary...")
fit_hyb_vocab <- with(imp_list_vocab_sst_1357, gam(as.formula(f_vocab_hyb), method = "REML"))

# 2. List Sort (intercept only)
f_listsort_hyb <- paste(
  "nc_y_nihtb__lswmt__uncor_score_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: List Sort...")
fit_hyb_listsort <- with(imp_list_listsort_157, gam(as.formula(f_listsort_hyb), method = "REML"))

# 3. Flanker (random slope)
f_flanker_hyb <- paste(
  "nc_y_nihtb__flnkr__uncor_score_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_slope
)
message("Fitting Hybrid GAMM: Flanker...")
fit_hyb_flanker <- with(imp_list_vocab_sst_1357, gam(as.formula(f_flanker_hyb), method = "REML"))

# 4. Card Sort (intercept only)
f_cardsort_hyb <- paste(
  "nc_y_nihtb__crdst__uncor_score_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: Card Sort...")
fit_hyb_cardsort <- with(imp_list_cardsort_17, gam(as.formula(f_cardsort_hyb), method = "REML"))

# 5. SST — refit on flipped scored variable (matching LMM)
f_sst_hyb <- paste(
  "sst_ssrt_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: SST...")
fit_hyb_sst <- with(imp_list_vocab_sst_1357, gam(as.formula(f_sst_hyb), method = "REML"))

# 6. Stroop RT — refit on flipped scored variable (matching LMM)
f_stroop_rt_hyb <- paste(
  "stroop_int_rt_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: Stroop RT...")
fit_hyb_stroop_rt <- with(imp_list_stroop_246, gam(as.formula(f_stroop_rt_hyb), method = "REML"))

# 7. Stroop Interference Accuracy (simple RE)
# NOTE: Uses the unflipped, unwinsorized interference accuracy z-score here.
# AIC comparison (Section below) shows linear age is preferred for this outcome,
# so this GAMM is used only for AIC comparison and is NOT used in final tables or
# plots — those use the LMM with the flipped, winsorized stroop_int_acc_z.
f_stroop_acc_hyb <- paste(
  "nc_y_est_interference_acc_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+ s(participant_id, bs='re')"
)
message("Fitting Hybrid GAMM: Stroop Acc...")
fit_hyb_stroop_acc <- with(imp_list_stroop_246, gam(as.formula(f_stroop_acc_hyb), method = "REML"))

# 8. N-Back 0-Back Accuracy (intercept only)
f_nback_0b_acc_hyb <- paste(
  "nback_0b_acc_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back 0-Back Acc...")
fit_hyb_nback_0b_acc <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_0b_acc_hyb), method = "REML"))

# 9. N-Back 0-Back RT (intercept only)
f_nback_0b_rt_hyb <- paste(
  "nback_0b_rt_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back 0-Back RT...")
fit_hyb_nback_0b_rt <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_0b_rt_hyb), method = "REML"))

# 10. N-Back 2-Back Accuracy (intercept only)
f_nback_2b_acc_hyb <- paste(
  "nback_2b_acc_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back 2-Back Acc...")
fit_hyb_nback_2b_acc <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_2b_acc_hyb), method = "REML"))

# 11. N-Back 2-Back RT (intercept only)
f_nback_2b_rt_hyb <- paste(
  "nback_2b_rt_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back 2-Back RT...")
fit_hyb_nback_2b_rt <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_2b_rt_hyb), method = "REML"))

# 12. N-Back Load Cost RT (intercept only)
f_nback_load_rt_hyb <- paste(
  "nback_load_rt_z ~",
  gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back Load Cost RT...")
fit_hyb_nback_load_rt <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_load_rt_hyb), method = "REML"))

# 13. N-Back Load Cost Accuracy (intercept only)
# f_nback_load_acc_hyb <- paste("nback_load_acc_z ~",
#                               gamm_terms_hybrid, "+", gamm_covars, "+", gamm_rand_intercept)
# message("Fitting Hybrid GAMM: N-Back Load Cost Acc...")
# fit_hyb_nback_load_acc <- with(imp_list_vocab_sst_1357, gam(as.formula(f_nback_load_acc_hyb), method = "REML"))

# 13. N-Back Load Cost Accuracy — LINEAR age (smooth not justified by AIC)
f_nback_load_acc_hyb <- paste(
  "nback_load_acc_z ~",
  "ab_p_demo_age_z + bilingual_use + bilingual_use:ab_p_demo_age_z +",
  gamm_covars, "+", gamm_rand_intercept
)
message("Fitting Hybrid GAMM: N-Back Load Cost Acc (linear age)...")
fit_hyb_nback_load_acc <- with(
  imp_list_vocab_sst_1357,
  gam(as.formula(f_nback_load_acc_hyb), method = "REML")
)

message("All GAMMs fitted.")

# --- AIC Comparison: Linear vs Smooth Age (first imputation only) ---
{
  all_hyb_fits <- list(
    "Vocab"           = fit_hyb_vocab,
    "WM"              = fit_hyb_listsort,
    "Flanker"         = fit_hyb_flanker,
    "Card Sort"       = fit_hyb_cardsort,
    "SST"             = fit_hyb_sst,
    "Stroop RT"       = fit_hyb_stroop_rt,
    "Stroop Acc"      = fit_hyb_stroop_acc,
    "N-Back 0b Acc"   = fit_hyb_nback_0b_acc,
    "N-Back 0b RT"    = fit_hyb_nback_0b_rt,
    "N-Back 2b Acc"   = fit_hyb_nback_2b_acc,
    "N-Back 2b RT"    = fit_hyb_nback_2b_rt,
    "N-Back Load RT"  = fit_hyb_nback_load_rt,
    "N-Back Load Acc" = fit_hyb_nback_load_acc
  )
  cat("\n--- Model Comparison: Linear Age vs. Smooth Age (AIC, first imputation) ---\n")
  aic_results <- purrr::map2_df(all_hyb_fits, names(all_hyb_fits), function(m_list, name) {
    m_smooth <- m_list[[1]]
    f_linear <- update(
      formula(m_smooth),
      . ~ . - s(ab_p_demo_age_z, k = 5) + ab_p_demo_age_z
    )
    m_linear <- gam(f_linear, data = m_smooth$model, method = "REML")
    delta_aic <- AIC(m_linear) - AIC(m_smooth)
    edf_age <- summary(m_smooth)$s.table["s(ab_p_demo_age_z)", "edf"]
    data.frame(
      Outcome   = name,
      EDF_Age   = round(edf_age, 2),
      Delta_AIC = round(delta_aic, 1),
      Decision  = ifelse(delta_aic > 10, "Use Smooth Age", "Use Linear Age")
    )
  })
  print(as.data.frame(aic_results))
}



# Example console output from the AIC comparison above (results from the
# authors' data; your output will differ):
#
#--- Model Comparison: Linear Age vs. Smooth Age (AIC, first imputation) ---
#           Outcome EDF_Age Delta_AIC       Decision
# 1            Vocab    3.83      81.7 Use Smooth Age
# 2               WM    3.87      95.1 Use Smooth Age
# 3          Flanker    3.81     225.1 Use Smooth Age
# 4        Card Sort    3.55      17.4 Use Smooth Age
# 5              SST    3.58      60.4 Use Smooth Age
# 6        Stroop RT    1.63      -0.3 Use Linear Age
# 7       Stroop Acc    1.04      -0.2 Use Linear Age
# 8    N-Back 0b Acc    3.80     211.0 Use Smooth Age
# 9     N-Back 0b RT    3.67      91.2 Use Smooth Age
# 10   N-Back 2b Acc    3.91     323.5 Use Smooth Age
# 11    N-Back 2b RT    2.99      10.3 Use Smooth Age
# 12  N-Back Load RT    3.28      40.2 Use Smooth Age
# 13 N-Back Load Acc    1.31      -0.4 Use Linear Age

# ==============================================================================
# SAVE GAMM FIT OBJECTS
# Add immediately after the "All GAMMs fitted." message in Section 0
# Saves all 13 fit lists so models never need to be refit
# ==============================================================================
save(
  fit_hyb_vocab,
  fit_hyb_listsort,
  fit_hyb_flanker,
  fit_hyb_cardsort,
  fit_hyb_sst,
  fit_hyb_stroop_rt,
  fit_hyb_stroop_acc,
  fit_hyb_nback_0b_acc,
  fit_hyb_nback_0b_rt,
  fit_hyb_nback_2b_acc,
  fit_hyb_nback_2b_rt,
  fit_hyb_nback_load_rt,
  fit_hyb_nback_load_acc,
  file = file.path(dir_proj, "gamm_fit_objects.RData")
)
message("GAMM fit objects saved to gamm_fit_objects.RData")

# Load the model objects (gamm fits)
# load("gamm_fit_objects.RData")

# ==============================================================================
# GAMM TABLE EXTRACTION — COMPLETE SELF-CONTAINED SCRIPT
# Run in a fresh R session. Do not load any RData files beforehand.
#
# Requires:
#   bilingual_model_objects.RData  — LMM fits + res objects
#   gamm_fit_objects.RData         — GAMM fits
#
# Produces:
#   gamm_results_and_table.RData   — results list + latex_table string
#   gamm_results_checkpoint.rds    — checkpoint after each GAMM outcome
# ==============================================================================
# ==============================================================================
# GAMM TABLE EXTRACTION — COMPLETE SELF-CONTAINED SCRIPT
# ==============================================================================
library(mgcv)
library(MuMIn)
library(mitml)
library(lme4)
library(purrr)
library(dplyr)
library(parallel)
setwd(dir_proj)
ncores <- 2

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

get_smooth_stats <- function(fit) {
  s_tab <- summary(fit)$s.table
  row <- s_tab[grepl("ab_p_demo_age_z", rownames(s_tab)), , drop = FALSE]
  if (nrow(row) == 0) {
    return(list(edf = NA, F = NA, p = NA))
  }
  list(
    edf = as.numeric(row[1, "edf"]),
    F   = as.numeric(row[1, "F"]),
    p   = as.numeric(row[1, "p-value"])
  )
}

get_param_stats <- function(fit, term_pattern) {
  if (inherits(fit, "lmerMod")) {
    p_tab <- summary(fit)$coefficients
    matched <- rownames(p_tab)[grep(term_pattern, rownames(p_tab))]
    if (length(matched) == 0) {
      return(NULL)
    }
    row <- p_tab[matched[1], , drop = FALSE]
    beta <- as.numeric(row[1, "Estimate"])
    se <- as.numeric(row[1, "Std. Error"])
    tval <- as.numeric(row[1, "t value"])
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
  } else {
    p_tab <- summary(fit)$p.table
    matched <- rownames(p_tab)[grep(term_pattern, rownames(p_tab))]
    if (length(matched) == 0) {
      return(NULL)
    }
    row <- p_tab[matched[1], , drop = FALSE]
    beta <- as.numeric(row[1, "Estimate"])
    se <- as.numeric(row[1, "Std. Error"])
    tval <- as.numeric(row[1, "t value"])
    pval <- as.numeric(row[1, "Pr(>|t|)"])
  }
  list(
    beta = beta, se = se, t = tval, p = pval,
    ci_l = beta - 1.96 * se,
    ci_u = beta + 1.96 * se,
    term = matched[1]
  )
}

drop_terms_from_formula <- function(f, terms_to_remove) {
  strip_term <- function(node, to_remove) {
    if (!is.call(node)) {
      dep <- paste(deparse(node, width.cutoff = 500L), collapse = " ")
      if (any(sapply(to_remove, function(tr) grepl(tr, dep, fixed = TRUE)))) {
        return(NULL)
      }
      return(node)
    }
    op <- paste(deparse(node[[1]], width.cutoff = 500L), collapse = " ")
    dep <- paste(deparse(node, width.cutoff = 500L), collapse = " ")
    if (op == "+") {
      left <- strip_term(node[[2]], to_remove)
      right <- strip_term(node[[3]], to_remove)
      if (is.null(left) && is.null(right)) {
        return(NULL)
      }
      if (is.null(left)) {
        return(right)
      }
      if (is.null(right)) {
        return(left)
      }
      call("+", left, right)
    } else {
      if (any(sapply(to_remove, function(tr) grepl(tr, dep, fixed = TRUE)))) {
        return(NULL)
      }
      node
    }
  }
  rhs_new <- strip_term(f[[3]], terms_to_remove)
  if (is.null(rhs_new)) {
    stop("All terms removed from formula.")
  }
  f[[3]] <- rhs_new
  f
}

get_delta_r2_gam <- function(fit, term_to_drop) {
  r2_full <- tryCatch(MuMIn::r.squaredGLMM(fit)[1, 1],
    error = function(e) {
      message("r.squaredGLMM failed (full)")
      NA_real_
    }
  )
  if (is.na(r2_full)) {
    return(NA_real_)
  }
  f_red <- tryCatch(drop_terms_from_formula(formula(fit), term_to_drop),
    error = function(e) {
      message("Formula drop failed: ", e$message)
      NULL
    }
  )
  if (is.null(f_red)) {
    return(NA_real_)
  }
  fit_red <- tryCatch(mgcv::gam(f_red, data = fit$model, method = "REML"),
    error = function(e) {
      message("Refit failed: ", e$message)
      NULL
    }
  )
  if (is.null(fit_red)) {
    return(NA_real_)
  }
  r2_red <- tryCatch(MuMIn::r.squaredGLMM(fit_red)[1, 1], error = function(e) NA_real_)
  max(0, r2_full - r2_red)
}

pool_parametric <- function(fit_list, term_pattern) {
  stats_list <- Filter(
    Negate(is.null),
    lapply(fit_list, get_param_stats, term_pattern = term_pattern)
  )
  m <- length(stats_list)
  betas <- sapply(stats_list, `[[`, "beta")
  ses <- sapply(stats_list, `[[`, "se")
  beta_pool <- mean(betas)
  w_var <- mean(ses^2)
  b_var <- var(betas)
  t_var <- w_var + (1 + 1 / m) * b_var
  se_pool <- sqrt(t_var)
  r_m <- (1 + 1 / m) * b_var / w_var
  df <- max((m - 1) * (1 + 1 / r_m)^2, 1)
  t_pool <- beta_pool / se_pool
  p_pool <- 2 * pt(abs(t_pool), df = df, lower.tail = FALSE)
  list(
    beta = beta_pool, se = se_pool, t = t_pool, p = p_pool,
    ci_l = beta_pool - qt(0.975, df) * se_pool,
    ci_u = beta_pool + qt(0.975, df) * se_pool,
    df = df
  )
}

pool_smooth <- function(fit_list) {
  stats_list <- lapply(fit_list, get_smooth_stats)
  list(
    edf = mean(sapply(stats_list, `[[`, "edf"), na.rm = TRUE),
    F   = mean(sapply(stats_list, `[[`, "F"), na.rm = TRUE),
    p   = mean(sapply(stats_list, `[[`, "p"), na.rm = TRUE)
  )
}

# --- Pool r_sp for LMM parametric terms via t-statistic method ---
# Used for bilingual_use and interaction in lmm_only outcomes (indices 6, 7, 13)
# to avoid lme4 marginality constraint inflating delta-R2 when bilingual_use
# main effect is dropped alongside its interaction term.
# r_sp = t / sqrt(t^2 + df); large-sample df = 10000 (conservative at n~11k)
pool_rsp_lmm_t <- function(fit_list, term_pattern) {
  stats_list <- Filter(
    Negate(is.null),
    lapply(fit_list, get_param_stats, term_pattern = term_pattern)
  )
  if (length(stats_list) == 0) {
    return(NA_real_)
  }
  t_vals <- sapply(stats_list, `[[`, "t")
  betas <- sapply(stats_list, `[[`, "beta")
  rsp_vals <- t_vals / sqrt(t_vals^2 + 10000)
  mean(rsp_vals, na.rm = TRUE)
}

pool_rsp_param <- function(fit_list, term_to_drop, term_pattern) {
  is_lmm <- inherits(fit_list[[1]], "lmerMod")
  delta_r2_vals <- as.numeric(mcmapply(function(f) {
    on.exit(gc())
    if (is_lmm) {
      r2_full <- tryCatch(MuMIn::r.squaredGLMM(f)[1, 1], error = function(e) NA_real_)
      if (is.na(r2_full)) {
        return(NA_real_)
      }
      f_red <- tryCatch(drop_terms_from_formula(formula(f), term_to_drop),
        error = function(e) NULL
      )
      if (is.null(f_red)) {
        return(NA_real_)
      }
      fit_red <- tryCatch(
        lme4::lmer(f_red,
          data = f@frame, REML = TRUE,
          control = lme4::lmerControl(optimizer = "bobyqa")
        ),
        error = function(e) NULL
      )
      if (is.null(fit_red)) {
        return(NA_real_)
      }
      r2_red <- tryCatch(MuMIn::r.squaredGLMM(fit_red)[1, 1], error = function(e) NA_real_)
      max(0, r2_full - r2_red)
    } else {
      get_delta_r2_gam(f, term_to_drop)
    }
  }, fit_list, mc.cores = ncores, SIMPLIFY = TRUE))
  betas <- as.numeric(mcmapply(function(f) {
    on.exit(gc())
    s <- get_param_stats(f, term_pattern)
    if (is.null(s)) NA_real_ else s$beta
  }, fit_list, mc.cores = ncores, SIMPLIFY = TRUE))
  sqrt(mean(delta_r2_vals, na.rm = TRUE)) * sign(mean(betas, na.rm = TRUE))
}

get_rsp_age_lp <- function(fit) {
  beta_age <- fixef(fit)["ab_p_demo_age_z"]
  if (is.na(beta_age)) {
    return(NA_real_)
  }
  age_vals <- fit@frame[["ab_p_demo_age_z"]]
  Y_age <- var(beta_age * age_vals)
  Y_fixed <- var(as.vector(model.matrix(fit) %*% fixef(fit)))
  vc <- as.data.frame(VarCorr(fit))
  Y_RE <- sum(vc$vcov[vc$grp != "Residual"])
  Y_R <- sigma(fit)^2
  Y_total <- Y_fixed + Y_RE + Y_R
  sqrt(max(0, Y_age / Y_total)) * sign(beta_age)
}

pool_rsp_age_lmm <- function(lmm_fit_list) {
  rsp_vals <- as.numeric(mcmapply(get_rsp_age_lp, lmm_fit_list,
    mc.cores = ncores, SIMPLIFY = TRUE
  ))
  mean(rsp_vals, na.rm = TRUE)
}

# ==============================================================================
# SECTION 2 + 3: SEQUENTIAL EXTRACTION
# ==============================================================================

names_registry <- c(
  "NIH Vocabulary", # 1  — GAMM smooth age
  "NIH List Sort (WM)", # 2  — GAMM smooth age
  "SST Inhibition", # 3  — GAMM smooth age
  "NIH Flanker", # 4  — GAMM smooth age
  "NIH Card Sort", # 5  — GAMM smooth age
  "Stroop Interference RT", # 6  — LMM only (linear age)
  "Stroop Interference Acc", # 7  — LMM only (linear age)
  "N-Back 0-Back Acc", # 8  — GAMM smooth age
  "N-Back 0-Back RT", # 9  — GAMM smooth age
  "N-Back 2-Back Acc", # 10 — GAMM smooth age
  "N-Back 2-Back RT", # 11 — GAMM smooth age
  "N-Back Load Cost RT", # 12 — GAMM smooth age
  "N-Back Load Cost Acc" # 13 — LMM only (linear age, AIC favored)
)

lmm_only_indices <- c(6, 7, 13)
gamm_indices <- c(1, 2, 3, 4, 5, 8, 9, 10, 11, 12)

gamm_name_map <- c(
  "1"  = "fit_hyb_vocab",
  "2"  = "fit_hyb_listsort",
  "3"  = "fit_hyb_sst",
  "4"  = "fit_hyb_flanker",
  "5"  = "fit_hyb_cardsort",
  "8"  = "fit_hyb_nback_0b_acc",
  "9"  = "fit_hyb_nback_0b_rt",
  "10" = "fit_hyb_nback_2b_acc",
  "11" = "fit_hyb_nback_2b_rt",
  "12" = "fit_hyb_nback_load_rt"
)

results <- vector("list", 13)

# ------------------------------------------------------------------------------
# PASS 1: LMM file only
# ------------------------------------------------------------------------------
message("=== PASS 1: LMM fits ===")
local({
  e <- new.env()
  load("bilingual_model_objects.RData", envir = e)
  lmm_fit_lists <- list(
    e$fit_vocab, # 1
    e$fit_listsort, # 2
    e$fit_sst, # 3
    e$fit_flanker, # 4
    e$fit_cardsort, # 5
    e$fit_stroop_int_rt, # 6
    e$fit_stroop_int_acc, # 7
    e$fit_nback_0b_acc, # 8
    e$fit_nback_0b_rt, # 9
    e$fit_nback_2b_acc, # 10
    e$fit_nback_2b_rt, # 11
    e$fit_nback_load_rt, # 12
    e$fit_nback_load_acc # 13
  )
  rm(e)
  gc()

  for (i in seq_along(lmm_fit_lists)) {
    message("  [", i, "/13] ", names_registry[i])
    lfl <- lmm_fit_lists[[i]]
    age_rsp <- tryCatch(pool_rsp_age_lmm(lfl),
      error = function(e) {
        message("  age_rsp failed: ", e$message)
        NA_real_
      }
    )

    if (i %in% lmm_only_indices) {
      age_p <- tryCatch(pool_parametric(lfl, "^ab_p_demo_age_z$"), error = function(e) NULL)
      bil_p <- tryCatch(pool_parametric(lfl, "^bilingual_use$"), error = function(e) NULL)
      int_p <- tryCatch(pool_parametric(lfl, "bilingual_use:ab_p_demo_age_z"),
        error = function(e) NULL
      )

      # Use t-based r_sp for bilingual_use and interaction to avoid lme4
      # marginality constraint inflating delta-R2 when main effect is dropped
      bil_r <- tryCatch(pool_rsp_lmm_t(lfl, "^bilingual_use$"),
        error = function(e) NA_real_
      )
      int_r <- tryCatch(pool_rsp_lmm_t(lfl, "bilingual_use:ab_p_demo_age_z"),
        error = function(e) NA_real_
      )

      results[[i]] <<- list(
        label = names_registry[i],
        age = list(
          type = "linear",
          edf = if (!is.null(age_p)) age_p$beta else NA,
          se = if (!is.null(age_p)) age_p$se else NA,
          stat = if (!is.null(age_p)) age_p$t else NA,
          p = if (!is.null(age_p)) age_p$p else NA,
          ci_l = if (!is.null(age_p)) age_p$ci_l else NA,
          ci_u = if (!is.null(age_p)) age_p$ci_u else NA,
          rsp = age_rsp
        ),
        bil = list(
          beta = if (!is.null(bil_p)) bil_p$beta else NA,
          se = if (!is.null(bil_p)) bil_p$se else NA,
          t = if (!is.null(bil_p)) bil_p$t else NA,
          p = if (!is.null(bil_p)) bil_p$p else NA,
          ci_l = if (!is.null(bil_p)) bil_p$ci_l else NA,
          ci_u = if (!is.null(bil_p)) bil_p$ci_u else NA,
          rsp = bil_r
        ),
        int = list(
          beta = if (!is.null(int_p)) int_p$beta else NA,
          se = if (!is.null(int_p)) int_p$se else NA,
          t = if (!is.null(int_p)) int_p$t else NA,
          p = if (!is.null(int_p)) int_p$p else NA,
          ci_l = if (!is.null(int_p)) int_p$ci_l else NA,
          ci_u = if (!is.null(int_p)) int_p$ci_u else NA,
          rsp = int_r
        )
      )
      message("    Full extraction complete (lmm_only)")
    } else {
      results[[i]] <<- list(label = names_registry[i], age_rsp = age_rsp)
      message("    age_rsp stored, GAMM pass will complete")
    }
    rm(lfl)
    gc()
  }
  message("=== PASS 1 complete ===")
})
gc()
saveRDS(results, "gamm_results_checkpoint_pass1.rds")
message("Pass 1 checkpoint saved.")

# ------------------------------------------------------------------------------
# PASS 2: GAMM file — one fit list at a time
# ------------------------------------------------------------------------------
message("=== PASS 2: GAMM fits (one at a time) ===")
for (i in gamm_indices) {
  obj_name <- gamm_name_map[as.character(i)]
  message("  [", i, "/13] ", names_registry[i], " (", obj_name, ")")
  local({
    e <- new.env()
    load("gamm_fit_objects.RData", envir = e)
    fl <- e[[obj_name]]
    rm(e)
    gc()

    age_s <- tryCatch(pool_smooth(fl), error = function(e) NULL)
    bil_p <- tryCatch(pool_parametric(fl, "^bilingual_use$"), error = function(e) NULL)
    bil_r <- tryCatch(pool_rsp_param(fl, "bilingual_use", "^bilingual_use$"),
      error = function(e) NA_real_
    )
    int_p <- tryCatch(pool_parametric(fl, "bilingual_use:ab_p_demo_age_z"),
      error = function(e) NULL
    )
    int_r <- tryCatch(
      pool_rsp_param(
        fl, "bilingual_use:ab_p_demo_age_z",
        "bilingual_use:ab_p_demo_age_z"
      ),
      error = function(e) NA_real_
    )
    age_rsp_stored <- results[[i]]$age_rsp

    results[[i]] <<- list(
      label = names_registry[i],
      age = list(
        type = "smooth",
        edf = if (!is.null(age_s)) age_s$edf else NA,
        se = NA,
        stat = if (!is.null(age_s)) age_s$F else NA,
        p = if (!is.null(age_s)) age_s$p else NA,
        ci_l = NA,
        ci_u = NA,
        rsp = age_rsp_stored
      ),
      bil = list(
        beta = if (!is.null(bil_p)) bil_p$beta else NA,
        se = if (!is.null(bil_p)) bil_p$se else NA,
        t = if (!is.null(bil_p)) bil_p$t else NA,
        p = if (!is.null(bil_p)) bil_p$p else NA,
        ci_l = if (!is.null(bil_p)) bil_p$ci_l else NA,
        ci_u = if (!is.null(bil_p)) bil_p$ci_u else NA,
        rsp = bil_r
      ),
      int = list(
        beta = if (!is.null(int_p)) int_p$beta else NA,
        se = if (!is.null(int_p)) int_p$se else NA,
        t = if (!is.null(int_p)) int_p$t else NA,
        p = if (!is.null(int_p)) int_p$p else NA,
        ci_l = if (!is.null(int_p)) int_p$ci_l else NA,
        ci_u = if (!is.null(int_p)) int_p$ci_u else NA,
        rsp = int_r
      )
    )
    message("    Complete")
  })
  gc()
  saveRDS(results, "gamm_results_checkpoint.rds")
  message("    Checkpoint saved [", i, "/13]")
}
message("=== PASS 2 complete. All extractions done. ===")

missing_idx <- which(sapply(results, function(r) is.null(r$age$p)))
if (length(missing_idx) > 0) {
  message("WARNING: incomplete entries at indices: ", paste(missing_idx, collapse = ", "))
} else {
  message("All 13 entries complete.")
}

# ==============================================================================
# SECTION 4: FDR CORRECTION
# ==============================================================================
all_p <- c(
  sapply(results, function(r) r$age$p),
  sapply(results, function(r) r$bil$p),
  sapply(results, function(r) r$int$p)
)
all_q <- p.adjust(all_p, method = "fdr")
n_outcomes <- length(results)
for (i in seq_along(results)) {
  results[[i]]$age$q <- all_q[i]
  results[[i]]$bil$q <- all_q[n_outcomes + i]
  results[[i]]$int$q <- all_q[2 * n_outcomes + i]
}

# ==============================================================================
# SECTION 5: LATEX TABLE
# ==============================================================================
fmt_p <- function(p) {
  if (is.na(p) || is.null(p)) {
    return("--")
  }
  if (p < 0.001) "$<$ .001" else sprintf("%.3f", p)
}
fmt_ci <- function(x) {
  if (is.na(x) || is.null(x)) {
    return("--")
  }
  s <- sprintf("%.3f", x)
  if (s == "-0.000") "0.000" else s
}
fmt_rsp <- function(x) {
  if (is.na(x) || is.null(x)) {
    return("--")
  }
  if (abs(x) < 0.0005) "0.000" else sprintf("%.3f", x)
}

age_rows <- paste(sapply(results, function(r) {
  ag <- r$age
  if (ag$type == "smooth") {
    sprintf(
      "  %s & %.2f & -- & %.3f & %s & -- & -- & %s & %s \\\\",
      r$label, ag$edf, ag$stat, fmt_p(ag$p), fmt_rsp(ag$rsp), fmt_p(ag$q)
    )
  } else {
    sprintf(
      "  %s & %.3f & %.3f & %.3f & %s & %s & %s & %s & %s \\\\",
      r$label, ag$edf, ag$se, ag$stat, fmt_p(ag$p),
      fmt_ci(ag$ci_l), fmt_ci(ag$ci_u), fmt_rsp(ag$rsp), fmt_p(ag$q)
    )
  }
}), collapse = "\n")

bil_rows <- paste(sapply(results, function(r) {
  bl <- r$bil
  sprintf(
    "  %s & %.3f & %.3f & %.3f & %s & %s & %s & %s & %s \\\\",
    r$label, bl$beta, bl$se, bl$t, fmt_p(bl$p),
    fmt_ci(bl$ci_l), fmt_ci(bl$ci_u), fmt_rsp(bl$rsp), fmt_p(bl$q)
  )
}), collapse = "\n")

int_rows <- paste(sapply(results, function(r) {
  it <- r$int
  sprintf(
    "  %s & %.3f & %.3f & %.3f & %s & %s & %s & %s & %s \\\\",
    r$label, it$beta, it$se, it$t, fmt_p(it$p),
    fmt_ci(it$ci_l), fmt_ci(it$ci_u), fmt_rsp(it$rsp), fmt_p(it$q)
  )
}), collapse = "\n")

latex_table <- paste0(
  "% ============================================================
% GAMM HYBRID TABLE: 13 outcomes x 3 predictors = 39 FDR tests
% ============================================================
\\begin{table*}[t]
\\centering
\\begin{minipage}{\\textwidth}
\\centering
\\caption{Associations of Cognitive Outcomes with Age and Bilingual Use}
\\small
\\begin{tabular}{lcccccccc}
\\toprule
Outcome & Est/EDF & $SE$ & \\textit{t/F} & \\textit{p} & 95\\% CI Lower & 95\\% CI Upper & $r_{sp}$ & FDR \\textit{q} \\\\
\\midrule
\\multicolumn{9}{l}{\\textbf{Main Effect of Age}} \\\\
\\midrule
", age_rows, "
\\midrule
\\multicolumn{9}{l}{\\textbf{Main Effect of Bilingual Use}} \\\\
\\midrule
", bil_rows, "
\\midrule
\\multicolumn{9}{l}{\\textbf{Interaction: Age $\\times$ Bilingual Use}} \\\\
\\midrule
", int_rows, "
\\bottomrule
\\end{tabular}
\\label{tab:hybrid_gamm}
\\smallskip
\\footnotesize
\\noindent\\caption*{\\textit{Note.} Total analytic sample $n = 11{,}868$. Estimates are pooled across 50 imputations via Rubin's rules. Age effects for NIH Toolbox and SST outcomes use a penalized regression spline ($k = 5$; reported as mean estimated degrees of freedom [EDF; numbers closer to 1.0 indicate a linear pattern] and mean $F$ across 50 imputations). $p$-values are the mean of the approximate Wald $p$-values across imputations. Stroop Interference and N-Back Load Cost Accuracy outcomes were modeled as linear mixed-effects models (LMMs) and report $\\beta$. \\textit{t} values are reported for $\\beta$. NIH Vocabulary = NIH Toolbox Picture Vocabulary; NIH List Sort WM = NIH Toolbox List Sorting Working Memory Task; NIH Flanker = NIH Toolbox Flanker; NIH Card Sort = NIH Toolbox Dimensional Change Card Sort; SSRT = Stop-Signal Reaction Time. Stroop Interference RT and Acc = Reaction Time and Accuracy from the Emotional Stroop Interference -- Congruent calculation. N-Back 0- and 2-Back RT and Acc = Reaction Time and Accuracy from the N-Back 0-Back and 2-Back Tasks. N-Back Load Cost RT = 0-Back minus 2-Back Reaction Time (z-standardized); N-Back Load Cost Acc = 2-Back minus 0-Back Accuracy (z-standardized); both scored so that higher = better. All analyses controlled for participant, family, and site as random effects, and sex, highest family education, family income, NIH Toolbox Oral Reading Recognition, and the top 10 genetic principal components as fixed effects. $\\beta$ = Standardized beta. SE = Standard error of $\\beta$. 95\\% CIs = 95\\% Confidence Intervals of the $\\beta$. $r_{sp}$ = Semipartial correlation. For smooth age terms, $r_{sp}$ uses the fixed-effect linear predictor variance method (Stoffel et al., 2021; Nakagawa et al., 2017). For bilingual use and interaction terms in GAMM outcomes, $r_{sp}$ uses $\\sqrt{\\overline{\\Delta R^2_m}}$ pooled across imputations, signed by $\\beta$. For bilingual use and interaction terms in LMM outcomes (Stroop Interference and N-Back Load Cost Accuracy), $r_{sp}$ is computed from pooled Wald $t$-statistics to avoid marginality constraint inflation in delta-$R^2$ estimation. FDR \\textit{q}-values are Benjamini-Hochberg adjusted across all 39 tests (13 outcomes $\\times$ 3 predictors).}
\\end{minipage}
\\end{table*}"
)

cat(latex_table, "\n")
message("Done.")

# ==============================================================================
# SAVE
# ==============================================================================
save(results, latex_table,
  file = "gamm_results_and_table.RData"
)
message("Saved: gamm_results_and_table.RData")


# ==============================================================================
# FINAL MASTER PLOTTING SCRIPT: ORIGINAL UNITS (AGE 9-18)
#
# RT outcomes (SST, Stroop RT, N-Back RT) use scale_y_reverse() so the plot
# reads as faster = higher visually while preserving true ms values on the axis.
# Accuracy and NIH Toolbox outcomes plot normally (higher = better).
# Load cost variables stay in z-score units (difference scores, no raw unit).
# Loads GAMM/LMM fits one at a time to avoid memory crash.
#
# Back-transformation for RT outcomes:
#   Model was fit on flipped z-scores (higher = faster = better).
#   Back-transform as: raw_ms = orig_mean - (pred_z * orig_sd)
#   This recovers true positive ms values. scale_y_reverse() then displays
#   faster (smaller ms) at top of axis.
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(lme4)

# --- 1. SETUP GLOBALS ---
age_mean_val <- mean(ds6_plotting_ref$ab_p_demo_age, na.rm = TRUE)
age_sd_val <- sd(ds6_plotting_ref$ab_p_demo_age, na.rm = TRUE)
z_at_9 <- (9 - age_mean_val) / age_sd_val
z_at_18 <- (18 - age_mean_val) / age_sd_val
age_seq <- seq(z_at_9, z_at_18, length.out = 50)

target_levels <- as.character(0:8)
legend_labels <- c(
  "0" = "0 (Strongly Bilingual)", "1" = "1", "2" = "2", "3" = "3",
  "4" = "4", "5" = "5", "6" = "6", "7" = "7", "8" = "8 (Monolingual)"
)
heat_palette <- c(
  "0" = "#0a2fa8", "1" = "#1a6ef5", "2" = "#56aaff", "3" = "#a8d4ff",
  "4" = "#f7f7f7", "5" = "#ffb380", "6" = "#ff6633", "7" = "#d42020", "8" = "#8b0000"
)

# raw:     column name in ds6_plotting_ref for back-transformation.
#          raw = NULL: plot in z-score units (load cost only).
# flip_y:  TRUE applies scale_y_reverse() so lower RT appears higher on axis.
#          Also triggers negated back-transformation: orig_mean - (pred_z * orig_sd)
#          so that model predictions (in flipped z-score space) map correctly
#          back to positive ms values where smaller = faster = better.
# cap_one: TRUE clamps y_plot to [0, 1] for proportion/accuracy outcomes.
stats_lookup <- list(
  "vocab" = list(
    raw     = "nc_y_nihtb__picvcb__uncor_score",
    ylab    = "NIH Toolbox Vocabulary",
    flip_y  = FALSE,
    cap_one = FALSE
  ),
  "listsort" = list(
    raw     = "nc_y_nihtb__lswmt__uncor_score",
    ylab    = "NIH Toolbox List Sorting",
    flip_y  = FALSE,
    cap_one = FALSE
  ),
  "flanker" = list(
    raw     = "nc_y_nihtb__flnkr__uncor_score",
    ylab    = "NIH Toolbox Flanker",
    flip_y  = FALSE,
    cap_one = FALSE
  ),
  "cardsort" = list(
    raw     = "nc_y_nihtb__crdst__uncor_score",
    ylab    = "NIH Toolbox Card Sort",
    flip_y  = FALSE,
    cap_one = FALSE
  ),
  "sst" = list(
    raw       = "mr_y_tfmri__sst__beh__ssrt_intgr_cleaned",
    ylab      = "Stop-Signal RT (ms)",
    flip_y    = TRUE,
    cap_one   = FALSE
  ),
  "stroop_rt" = list(
    raw       = "nc_y_est_interference_rt",
    ylab      = "Stroop Interference RT (ms)",
    flip_y    = TRUE,
    cap_one   = FALSE
  ),
  "stroop_acc" = list(
    raw     = "nc_y_est_interference_acc",
    ylab    = "Stroop Interference Accuracy",
    flip_y  = FALSE,
    cap_one = TRUE
  ),
  "nback_0b_acc" = list(
    raw     = "mr_y_tfmri__nback__beh__0b_acc_cleaned",
    ylab    = "EN-Back 0-Back Accuracy",
    flip_y  = FALSE,
    cap_one = TRUE
  ),
  "nback_0b_rt" = list(
    raw     = "mr_y_tfmri__nback__beh__0b_rt_cleaned",
    ylab    = "EN-Back 0-Back RT (ms)",
    flip_y  = TRUE,
    cap_one = FALSE
  ),
  "nback_2b_acc" = list(
    raw     = "mr_y_tfmri__nback__beh__2b_acc_cleaned",
    ylab    = "EN-Back 2-Back Accuracy",
    flip_y  = FALSE,
    cap_one = TRUE
  ),
  "nback_2b_rt" = list(
    raw     = "mr_y_tfmri__nback__beh__2b_rt_cleaned",
    ylab    = "EN-Back 2-Back RT (ms)",
    flip_y  = TRUE,
    cap_one = FALSE
  ),
  "nback_load_rt" = list(
    raw       = NULL,
    ylab      = "EN-Back Cost RT (Z) [Higher = Less Slowing]",
    flip_y    = FALSE,
    cap_one   = FALSE,
    ylab_size = 17
  ),
  "nback_load_acc" = list(
    raw       = NULL,
    ylab      = "EN-Back Cost Accuracy (Z) [Higher = Less Drop]",
    flip_y    = FALSE,
    cap_one   = FALSE,
    ylab_size = 17
  )
)

# obj_name: object name inside the RData file
# rdata:    "gamm" = gamm_fit_objects.RData
#           "lmm"  = bilingual_model_objects.RData
model_configs <- list(
  list(obj_name = "fit_hyb_vocab", name = "vocab", rdata = "gamm"),
  list(obj_name = "fit_hyb_listsort", name = "listsort", rdata = "gamm"),
  list(obj_name = "fit_hyb_flanker", name = "flanker", rdata = "gamm"),
  list(obj_name = "fit_hyb_cardsort", name = "cardsort", rdata = "gamm"),
  list(obj_name = "fit_hyb_sst", name = "sst", rdata = "gamm"),
  list(obj_name = "fit_stroop_int_rt", name = "stroop_rt", rdata = "lmm"),
  list(obj_name = "fit_stroop_int_acc", name = "stroop_acc", rdata = "lmm"),
  list(obj_name = "fit_hyb_nback_0b_acc", name = "nback_0b_acc", rdata = "gamm"),
  list(obj_name = "fit_hyb_nback_0b_rt", name = "nback_0b_rt", rdata = "gamm"),
  list(obj_name = "fit_hyb_nback_2b_acc", name = "nback_2b_acc", rdata = "gamm"),
  list(obj_name = "fit_hyb_nback_2b_rt", name = "nback_2b_rt", rdata = "gamm"),
  list(obj_name = "fit_hyb_nback_load_rt", name = "nback_load_rt", rdata = "gamm"),
  list(obj_name = "fit_nback_load_acc", name = "nback_load_acc", rdata = "lmm")
)

gamm_path <- file.path(dir_proj, "gamm_fit_objects.RData")
lmm_path <- file.path(dir_proj, "bilingual_model_objects.RData")

# --- 2. EXECUTION LOOP ---
for (cfg in model_configs) {
  message("\n-----------------------------------------------")
  message("PROCESSING: ", cfg$name)

  # Load only the one fit list needed for this outcome
  fit_list <- local({
    e <- new.env()
    path <- if (cfg$rdata == "gamm") gamm_path else lmm_path
    load(path, envir = e)
    fl <- e[[cfg$obj_name]]
    rm(e)
    gc()
    fl
  })

  if (is.null(fit_list)) {
    message("Object not found. Skipping.")
    next
  }

  meta <- stats_lookup[[cfg$name]]
  model_obj <- fit_list[[1]]
  full_data <- model.frame(model_obj)
  use_raw <- !is.null(meta$raw)

  # A. Scaling stats for back-transformation
  # Raw variables are in original units — no pre-flipping.
  # orig_mean and orig_sd are always computed from the unflipped raw variable.
  # For RT outcomes the negated back-transform in step E handles the direction.
  if (use_raw) {
    raw_vals <- ds6_plotting_ref[[meta$raw]]
    orig_mean <- mean(raw_vals, na.rm = TRUE)
    orig_sd <- sd(raw_vals, na.rm = TRUE)
  }

  # B. POPULATION TREND GRID
  ref_row <- full_data[1, , drop = FALSE]
  grp_means <- expand_grid(ab_p_demo_age_z = age_seq, bilingual_use = 0:8)
  other_vars <- setdiff(
    names(full_data),
    c("ab_p_demo_age_z", "bilingual_use", names(full_data)[1])
  )
  for (v in other_vars) grp_means[[v]] <- ref_row[[v]]

  # C. SHIMMER DATA
  set.seed(42)
  df_shimmer <- full_data %>%
    mutate(bilingual_round = round(pmax(0, pmin(8, bilingual_use)))) %>%
    group_by(bilingual_round) %>%
    filter(participant_id %in% sample(
      unique(participant_id),
      min(200, n_distinct(participant_id))
    )) %>%
    ungroup()

  # D. PREDICTIONS
  if (inherits(model_obj, "gam")) {
    re_labels <- sapply(model_obj$smooth, function(x) {
      if (inherits(x, "random.effect") || inherits(x, "re.smooth")) x$label else NULL
    }) %>% unlist()
    grp_means$pred_z <- predict(model_obj,
      newdata = grp_means,
      exclude = re_labels, newdata.guaranteed = TRUE
    )
    df_shimmer$pred_z <- predict(model_obj,
      newdata = df_shimmer,
      newdata.guaranteed = TRUE
    )
  } else {
    grp_means$pred_z <- predict(model_obj,
      newdata = grp_means,
      re.form = NA, allow.new.levels = TRUE
    )
    df_shimmer$pred_z <- predict(model_obj,
      newdata = df_shimmer,
      re.form = NULL, allow.new.levels = TRUE
    )
  }

  # E. BACK-TRANSFORMATION & AGE FILTERING
  # RT outcomes (flip_y = TRUE): model fit on flipped z-scores (higher = faster).
  #   back-transform: orig_mean - (pred_z * orig_sd)
  #   -> recovers positive ms; larger pred_z (faster) gives smaller ms value.
  #   scale_y_reverse() then displays smaller ms (faster) at top of axis.
  # All other outcomes: standard (pred_z * orig_sd) + orig_mean
  back_transform <- function(pz) {
    if (isTRUE(meta$flip_y)) {
      orig_mean - (pz * orig_sd)
    } else {
      (pz * orig_sd) + orig_mean
    }
  }

  if (use_raw) {
    df_shimmer <- df_shimmer %>%
      mutate(
        age_years       = (ab_p_demo_age_z * age_sd_val) + age_mean_val,
        y_plot          = back_transform(pred_z),
        bilingual_use_f = factor(as.character(bilingual_round), levels = target_levels)
      ) %>%
      mutate(
        y_plot = if (isTRUE(meta$cap_one)) pmin(pmax(y_plot, 0), 1) else y_plot
      ) %>%
      filter(age_years >= 9)

    grp_means_smooth <- grp_means %>%
      mutate(
        age_years       = (ab_p_demo_age_z * age_sd_val) + age_mean_val,
        y_plot          = back_transform(pred_z),
        bilingual_use_f = factor(as.character(bilingual_use), levels = target_levels)
      ) %>%
      mutate(
        y_plot = if (isTRUE(meta$cap_one)) pmin(pmax(y_plot, 0), 1) else y_plot
      )
  } else {
    # Load cost outcomes: no back-transformation, stay in z-score units
    df_shimmer <- df_shimmer %>%
      mutate(
        age_years       = (ab_p_demo_age_z * age_sd_val) + age_mean_val,
        y_plot          = pred_z,
        bilingual_use_f = factor(as.character(bilingual_round), levels = target_levels)
      ) %>%
      filter(age_years >= 9)

    grp_means_smooth <- grp_means %>%
      mutate(
        age_years       = (ab_p_demo_age_z * age_sd_val) + age_mean_val,
        y_plot          = pred_z,
        bilingual_use_f = factor(as.character(bilingual_use), levels = target_levels)
      )
  }

  # F. Y-AXIS LIMITS
  trend_range <- diff(range(grp_means_smooth$y_plot, na.rm = TRUE))
  y_lower <- min(grp_means_smooth$y_plot) - (1.0 * trend_range)
  y_upper <- max(grp_means_smooth$y_plot) + (1.0 * trend_range)

  # G. PLOT CONSTRUCTION
  p <- ggplot() +
    geom_line(
      data = df_shimmer,
      aes(x = age_years, y = y_plot, group = participant_id, color = bilingual_use_f),
      alpha = 0.20, linewidth = 0.25, show.legend = FALSE
    ) +
    geom_line(
      data = grp_means_smooth,
      aes(x = age_years, y = y_plot, color = bilingual_use_f, group = bilingual_use_f),
      linewidth = 1.4
    ) +
    scale_color_manual(
      values = heat_palette, labels = legend_labels,
      name = "Bilingual Use", breaks = target_levels,
      guide = guide_legend(reverse = TRUE)
    ) +
    scale_x_continuous(limits = c(9, 18), breaks = seq(9, 17, by = 2)) +
    theme_minimal(base_size = 20) +
    labs(x = "Age (Years)", y = meta$ylab) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      axis.line          = element_line(color = "black"),
      legend.position    = "right"
    )

  # Override y-axis label font size for outcomes with long labels
  if (!is.null(meta$ylab_size)) {
    p <- p + theme(axis.title.y = element_text(size = meta$ylab_size))
  }

  # Apply reversed y-axis for RT outcomes and set limits
  if (isTRUE(meta$flip_y)) {
    p <- p + scale_y_reverse(limits = c(y_upper, y_lower))
  } else {
    p <- p + scale_y_continuous(limits = c(y_lower, y_upper))
  }

  print(p)
  ggsave(paste0("PLOT_RAW_UNITS_", cfg$name, ".tiff"),
    plot = p, width = 10, height = 6, dpi = 300, compression = "lzw"
  )
  message("Saved: PLOT_RAW_UNITS_", cfg$name, ".tiff")

  # Release fit list before next iteration
  rm(fit_list, model_obj, full_data, df_shimmer, grp_means, grp_means_smooth, p)
  gc()
}

message("\nAll plots complete.")




# ==============================================================================
# SIMPLIFIED BILINGUAL EFFECT TABLE: Per Unit and Total Range (Score 0 to Score 8)
# Beta (z-score units) * orig_sd = effect in original units per 1-unit change
# Multiply by 8 for total range (bilingual=0 to monolingual=8)
# Load cost outcomes excluded — z-score difference scores with no raw unit
# ==============================================================================

sd_lookup <- list(
  "NIH Vocabulary"               = sd(ds6_plotting_ref$nc_y_nihtb__picvcb__uncor_score, na.rm = TRUE),
  "NIH List Sorting"             = sd(ds6_plotting_ref$nc_y_nihtb__lswmt__uncor_score, na.rm = TRUE),
  "SST Inhibition"               = sd(ds6_plotting_ref$mr_y_tfmri__sst__beh__ssrt_intgr_cleaned, na.rm = TRUE),
  "NIH Flanker"                  = sd(ds6_plotting_ref$nc_y_nihtb__flnkr__uncor_score, na.rm = TRUE),
  "NIH Card Sort"                = sd(ds6_plotting_ref$nc_y_nihtb__crdst__uncor_score, na.rm = TRUE),
  "Stroop Interference RT"       = sd(ds6_plotting_ref$nc_y_est_interference_rt, na.rm = TRUE),
  "Stroop Interference Accuracy" = sd(ds6_plotting_ref$nc_y_est_interference_acc, na.rm = TRUE),
  "EN-Back 0-Back Accuracy"      = sd(ds6_plotting_ref$mr_y_tfmri__nback__beh__0b_acc_cleaned, na.rm = TRUE),
  "EN-Back 0-Back RT"            = sd(ds6_plotting_ref$mr_y_tfmri__nback__beh__0b_rt_cleaned, na.rm = TRUE),
  "EN-Back 2-Back Accuracy"      = sd(ds6_plotting_ref$mr_y_tfmri__nback__beh__2b_acc_cleaned, na.rm = TRUE),
  "EN-Back 2-Back RT"            = sd(ds6_plotting_ref$mr_y_tfmri__nback__beh__2b_rt_cleaned, na.rm = TRUE)
)

units_lookup <- list(
  "NIH Vocabulary"               = "",
  "NIH List Sorting"             = "",
  "SST Inhibition"               = " ms",
  "NIH Flanker"                  = "",
  "NIH Card Sort"                = "",
  "Stroop Interference RT"       = " ms",
  "Stroop Interference Accuracy" = "",
  "EN-Back 0-Back Accuracy"      = "",
  "EN-Back 0-Back RT"            = " ms",
  "EN-Back 2-Back Accuracy"      = "",
  "EN-Back 2-Back RT"            = " ms"
)

# Display labels for the table (separate from internal r$label values)
display_label_lookup <- list(
  "NIH Vocabulary"               = "NIH Toolbox Vocabulary",
  "NIH List Sorting"             = "NIH Toolbox List Sorting",
  "SST Inhibition"               = "Stop-Signal RT",
  "NIH Flanker"                  = "NIH Toolbox Flanker",
  "NIH Card Sort"                = "NIH Toolbox Card Sort",
  "Stroop Interference RT"       = "Stroop Interference RT",
  "Stroop Interference Accuracy" = "Stroop Interference Accuracy",
  "EN-Back 0-Back Accuracy"      = "EN-Back 0-Back Accuracy",
  "EN-Back 0-Back RT"            = "EN-Back 0-Back RT",
  "EN-Back 2-Back Accuracy"      = "EN-Back 2-Back Accuracy",
  "EN-Back 2-Back RT"            = "EN-Back 2-Back RT"
)

# Display order for rows in the output table
display_order <- c(
  "NIH Vocabulary",
  "NIH List Sorting",
  "SST Inhibition",
  "NIH Flanker",
  "NIH Card Sort",
  "Stroop Interference Accuracy",
  "Stroop Interference RT",
  "EN-Back 0-Back Accuracy",
  "EN-Back 0-Back RT",
  "EN-Back 2-Back Accuracy",
  "EN-Back 2-Back RT"
)

exclude_labels <- c("EN-Back Load Cost RT", "EN-Back Load Cost Accuracy")

# Build a name-keyed list of results for the included outcomes, then loop in display_order
results_by_label <- setNames(
  Filter(function(r) !r$label %in% exclude_labels, results),
  sapply(Filter(function(r) !r$label %in% exclude_labels, results), function(r) r$label)
)

table_rows <- sapply(display_order, function(internal_label) {
  r <- results_by_label[[internal_label]]
  if (is.null(r)) stop(sprintf("Missing results entry for: %s", internal_label))
  beta <- r$bil$beta
  orig_sd <- sd_lookup[[internal_label]]
  units <- units_lookup[[internal_label]]
  display_label <- display_label_lookup[[internal_label]]
  if (is.null(orig_sd)) stop(sprintf("Missing sd_lookup entry for: %s", internal_label))
  if (is.null(units)) stop(sprintf("Missing units_lookup entry for: %s", internal_label))
  if (is.null(display_label)) stop(sprintf("Missing display_label_lookup entry for: %s", internal_label))
  per_unit <- beta * orig_sd
  total <- per_unit * 8
  sprintf(
    "  %s & %.3f%s & %.3f%s \\\\",
    display_label, per_unit, units, total, units
  )
})

latex_simple_table <- paste0(
  "\\begin{table*}[h!]
\\centering
\\caption{Predicted Differences on Cognitive Outcomes Across the Bilingual Use Range (Score 0 to Score 8)}
\\label{tab:simplified-bilingual}
\\begin{tabular}{lcc}
\\hline
\\textbf{Measure} & \\textbf{Effect per Unit Change} & \\textbf{Total Effect (0 compared to 8)} \\\\ \\hline
", paste(table_rows, collapse = "\n"), "
\\hline
\\end{tabular}
\\caption*{\\textit{Note.} Bilingual Use was coded on a scale of 0 (exclusively non-English use across friends and family) to 8 (exclusively English use across both contexts). Effects are presented in the original units of each outcome and follow the higher $=$ better sign convention used throughout: positive values indicate higher predicted scores among more monolingually-experienced children, and negative values indicate higher predicted scores among more bilingually-experienced children. Units are NIH Toolbox standardized scores for NIH Toolbox outcomes (mean 100, $SD = 15$), milliseconds for RT outcomes, and proportion correct for accuracy outcomes. NIH Toolbox Vocabulary = NIH Toolbox Picture Vocabulary Test; NIH Toolbox List Sorting = NIH Toolbox List Sorting Working Memory Test; NIH Toolbox Flanker = NIH Toolbox Flanker Inhibitory Control and Attention Test; NIH Toolbox Card Sort = NIH Toolbox Dimensional Change Card Sort; Stop-Signal RT = Stop-Signal Reaction Time (SSRT) from the Stop-Signal Task. Stroop Interference RT and Stroop Interference Accuracy = the incongruent-minus-congruent reaction time and accuracy difference scores from the Emotional Stroop task. EN-Back 0-Back and 2-Back RT and Accuracy = reaction time and accuracy on the 0-back and 2-back conditions of the Emotional $n$-back task. EN-Back Load Cost outcomes are excluded because they are computed from $z$-standardized components and have no meaningful original-unit scale.}
\\end{table*}"
)

cat(latex_simple_table)

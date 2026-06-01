# ===========================================================================
# Pipeline configuration
# Single source of truth for paths, sessions, variable-name mappings, and
# model settings used across the targets pipeline. Sourced via tar_source().
# ===========================================================================

pipeline_config <- function() {
  list(
    # --- Sessions / waves (mirror python/config.py SESSIONS) ----------------
    sessions = c("ses-02A", "ses-06A"),
    session_labels = c(
      "ses-02A" = "yr2", 
      "ses-04A" = "yr4", 
      "ses-06A" = "yr6"
    ),

    # --- Directories --------------------------------------------------------
    handoff_dir = "data/handoff",
    # Pre-assembled NBDC/ABCD dataset (read by NBDCtools; not read directly here)
    dir_abcd    = "dev/data/fitbit-summaries/7_0/0_1/data",
    out_fig_dir = "outputs/figures",
    out_tab_dir = "outputs/tables",

    # --- Variables of interest (passed straight to NBDCtools `vars`) --------
    # List the variable names you want. To instead pull a whole curated list
    # from a CSV, inline-read its `name` column here, e.g.:
    vars_mh = readr::read_csv(
      "data/vars-of-interest/psychiatric_disorders-diagnosis-sim04.csv"
    )$name,
    # vars_mh = c(
    #   "mh_y_ksads__dep__mdd__pres_dx",       # major depressive disorder
    #   "mh_y_ksads__dep__unspec__pres_dx",    # unspecified depressive disorder
    #   "mh_y_ksads__dmdd__dmdd__pres_dx",     # disruptive mood dysregulation
    #   "mh_y_ksads__panic__oth__pres_dx",     # panic / anxiety
    #   "mh_y_ksads__bpd__bpd1__curdep__pres_dx",  # bipolar I, current depressed
    #   "mh_y_ksads__bpd__bpd2__rcnt__hypomix__pres_dx"  # bipolar II
    # ),
    # Metabolic labs — list the variables once finalized (empty = skip).
    vars_metabolic = c(
      "ph_y_bld__rslt__hdl_qnt",  # HDL cholesterol
      "ph_y_bld__rslt__chol_qnt", # Total cholesterol
      "ph_y_bld__rslt__a1c_qnt"   # Hemoglobin a1c
    ),

    # --- Demographic / covariate NBDCtools variable names -------------------
    # These mirror scripts/01_gather_metadata.R.
    var_sex     = "ab_g_stc__cohort_sex",
    var_site    = "ab_g_dyn__design_site",
    var_age     = "ab_g_dyn__visit_age",
    # Pubertal stage: sex-specific approximate Tanner categories (PDS-derived).
    # Each participant has a value in exactly one (the other is NA); the loader
    # coalesces them into a single `pubertal_stage` (males -> m, females -> f).
    var_puberty_m = "ph_y_pds__m_categ",
    var_puberty_f = "ph_y_pds__f_categ",

    # --- Outcome model settings --------------------------------------------
    # Mental health: KSADS variables are categorical present/past diagnosis
    # flags. Primary outcome = total current-diagnosis burden (count of
    # "present" diagnoses); secondary = binary "any depression diagnosis".
    mh_present_suffix = "pres_dx",       # which KSADS columns count toward burden
    mh_depression_grep = "ksads__dep",   # pattern for the depression-specific binary

    # Metabolic outcome columns as they appear in master (i.e. the NBDCtools
    # column names from vars_metabolic above). Split into:
    #   metabolic_outcome_vars : all metabolic outcomes used in models / CCA
    #   metabolic_log_vars     : subset to log-transform before modelling
    metabolic_outcome_vars = c(
      "ph_y_bld__rslt__hdl_qnt",
      "ph_y_bld__rslt__chol_qnt",
      "ph_y_bld__rslt__a1c_qnt"
    ),
    metabolic_log_vars = character(0),  # none of these require log-transform

    # Primary behavioral predictors for the feature-based outcome models.
    # Within-wave z-scored features (consistent meaning across waves), so
    # `feature * wave_age` interactions are interpretable. Pre-specify a small
    # set to limit multiple comparisons.
    primary_features = c(
      "z_min_slp_onset_sd",          # sleep-onset variability (irregularity)
      "z_steps_total_mean_daily",    # activity volume
      "z_steps_total_mvpa_fraction", # moderate-to-vigorous activity
      "z_min_slp_duration_mean"      # sleep duration
    ),

    # --- QC -----------------------------------------------------------------
    fig_dpi = 300
  )
}

# Map a session_id vector to wave labels (yr2 / yr6 / ...)
session_to_wave <- function(session_id, cfg = pipeline_config()) {
  unname(cfg$session_labels[as.character(session_id)])
}

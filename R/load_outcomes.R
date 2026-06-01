# ===========================================================================
# Load outcome & covariate data via NBDCtools
# Mirrors the access pattern in scripts/01_gather_metadata.R. These functions
# pull from the pre-assembled ABCD dataset (cfg$dir_abcd) and restrict to the
# active sessions, mapping session_id -> wave_label to match the handoff files.
# ===========================================================================

.filter_sessions <- function(df, cfg) {
  df |>
    dplyr::filter(.data$session_id %in% cfg$sessions) |>
    dplyr::mutate(wave_label = session_to_wave(.data$session_id, cfg),
                  .after = "session_id")
}

# --- Mental health (KSADS diagnoses) ---------------------------------------
# KSADS variables are categorical present/past diagnosis flags. We derive:
#   mh_burden_present : count of "present" current diagnoses (primary outcome)
#   mh_any_depression : binary, any depression diagnosis present (secondary)
load_mh_outcomes <- function(cfg = pipeline_config()) {
  if (length(cfg$vars_mh) == 0) {
    stop("cfg$vars_mh is empty — list the mental-health variables of interest ",
         "in R/config_pipeline.R.", call. = FALSE)
  }

  raw <- NBDCtools::create_dataset_abcd(dir = cfg$dir_abcd, vars = cfg$vars_mh) |>
    .filter_sessions(cfg)

  present_cols <- grep(paste0(cfg$mh_present_suffix, "$"), names(raw), value = TRUE)
  dep_cols     <- grep(cfg$mh_depression_grep, present_cols, value = TRUE)

  # KSADS present/absent encoded as character; treat "1"/"present"/"yes" as present.
  is_present <- function(x) as.integer(as.character(x) %in% c("1", "present", "yes", "TRUE"))

  raw |>
    dplyr::mutate(
      mh_burden_present = rowSums(
        dplyr::across(dplyr::all_of(present_cols), is_present), na.rm = TRUE
      ),
      mh_any_depression = as.integer(
        rowSums(dplyr::across(dplyr::all_of(dep_cols), is_present), na.rm = TRUE) > 0
      )
    ) |>
    dplyr::select(participant_id, session_id, wave_label,
                  mh_burden_present, mh_any_depression,
                  dplyr::all_of(present_cols))
}

# --- Metabolic labs --------------------------------------------------------
# Uses cfg$vars_metabolic (a character vector). Returns the raw lab columns;
# log-transforms and BMI-SDS handling happen in the model fn. If no metabolic
# variables are listed yet, returns just the keys so the join still works and
# the downstream metabolic models skip gracefully.
load_metabolic_outcomes <- function(cfg = pipeline_config()) {
  if (length(cfg$vars_metabolic) == 0) {
    message("load_metabolic_outcomes: cfg$vars_metabolic is empty — ",
            "returning keys only (metabolic models will be skipped).")
    return(tibble::tibble(participant_id = character(),
                          session_id = character(),
                          wave_label = character()))
  }

  NBDCtools::create_dataset_abcd(dir = cfg$dir_abcd, vars = cfg$vars_metabolic) |>
    .filter_sessions(cfg)
}

# --- Demographics + pubertal stage -----------------------------------------
load_demographics <- function(cfg = pipeline_config()) {
  raw <- NBDCtools::create_dataset_abcd(
    dir = cfg$dir_abcd,
    vars = c(cfg$var_sex, cfg$var_site, cfg$var_age,
             cfg$var_puberty_m, cfg$var_puberty_f,
             "ab_g_stc__cohort_ethn", "ab_g_stc__cohort_race__nih",
             "ab_g_dyn__cohort_income__hhold__3lvl", "ab_g_dyn__cohort_edu__cgs")
  ) |>
    format_sociodemo_vars() |>     # from R/create_alluvial_plot.R
    .filter_sessions(cfg)

  # The two sex-specific Tanner-category columns must be present. If absent,
  # the derived score has not been added — run ABCDscores::add_scores() first.
  for (v in c(cfg$var_puberty_m, cfg$var_puberty_f)) {
    if (!v %in% names(raw)) {
      stop(sprintf(paste0("Puberty score column '%s' not found. These are ",
                          "derived PDS scores; ensure they are computed (e.g. ",
                          "via ABCDscores::add_scores()) before loading."), v),
           call. = FALSE)
    }
  }

  # Coalesce the sex-specific approximate Tanner stages into one ordinal scale.
  # Each participant has a value in exactly one column; 777/999 (declined/
  # don't-know) are treated as missing.
  out <- raw |>
    dplyr::mutate(
      wave_age = .data[[cfg$var_age]],
      pubertal_stage = dplyr::coalesce(
        dplyr::na_if(dplyr::na_if(as.numeric(.data[[cfg$var_puberty_m]]), 777), 999),
        dplyr::na_if(dplyr::na_if(as.numeric(.data[[cfg$var_puberty_f]]), 777), 999)
      )
    )

  # Pubertal stage is a required confounder — halt if entirely missing.
  if (all(is.na(out$pubertal_stage))) {
    stop(paste0("Pubertal stage is entirely missing after coalescing ",
                cfg$var_puberty_m, " and ", cfg$var_puberty_f,
                ". Confirm the scores are present before running outcome models."),
         call. = FALSE)
  }

  out |>
    dplyr::select(participant_id, session_id, wave_label,
                  sex, wave_age, pubertal_stage,
                  ethn, race_nih, hhold_3lvl, edu)
}

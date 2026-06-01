# ===========================================================================
# Phase 6.2 — Metabolic outcome models (plan §6.2)
# Same wave-nested cluster specification as the MH models (cluster effects are
# estimated within wave; no cross-wave label equivalence assumed). Skewed labs
# (config metabolic_log_vars) are log-transformed. One tidy fixed-effects table
# per available outcome.
# ===========================================================================

.fixef_table_metab <- function(model) {
  cf <- summary(model)$coefficients
  tibble::as_tibble(cf, rownames = "term")
}

fit_metabolic_models <- function(master, cfg = pipeline_config()) {
  outcomes <- intersect(cfg$metabolic_outcome_vars, names(master))

  if (length(outcomes) == 0) {
    message("fit_metabolic_models: no configured metabolic outcomes present in ",
            "master — skipping. Define cfg$vars_metabolic_file and column names.")
    return(list())
  }

  fit_one <- function(var) {
    df <- master
    if (length(cfg$metabolic_log_vars) > 0 && var %in% cfg$metabolic_log_vars) {
      df[[var]] <- log(df[[var]])
    }
    m <- lmerTest::lmer(
      stats::reformulate(
        c("wave_label", "wave_label:cluster_aligned", "sex", "pubertal_stage",
          "wear_time_fraction", "(1 | participant_id)"),
        response = var
      ),
      data = df, REML = TRUE
    )
    .fixef_table_metab(m)
  }

  stats::setNames(lapply(outcomes, fit_one), outcomes)
}

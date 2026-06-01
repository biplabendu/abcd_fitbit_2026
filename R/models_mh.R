# ===========================================================================
# Phase 6.1 — Mental health outcome models (plan §6.1)
# KSADS outcomes are categorical, so we model:
#   primary   : diagnosis burden (count of present current diagnoses) via LMM
#   secondary : any depression diagnosis (binary) via logistic GLMM
#
# WAVE-NESTED cluster specification.  Because clusters are derived independently
# per wave and the cross-wave alignment is only a best-match (see Supp. cluster
# alignment), a pooled `cluster * wave_age` interaction would assume a cluster
# label means the same phenotype at both ages — which it does not. We therefore
# estimate cluster effects *within* each wave by nesting cluster in wave
# (`wave_label + wave_label:cluster_aligned`). No cross-wave equivalence is
# assumed; the participant random intercept still absorbs repeated measures.
#
# The developmental (× age) question is handled by the feature-based models
# (models_features.R) and the trajectory-group model, whose predictors keep a
# consistent meaning across waves.
# ===========================================================================

.fixef_table <- function(model) {
  cf <- summary(model)$coefficients
  tibble::as_tibble(cf, rownames = "term")
}

# Tidy fixed-effects tables for both MH models.
fit_mh_models <- function(master) {
  out <- list()

  if ("mh_burden_present" %in% names(master)) {
    m_burden <- lmerTest::lmer(
      mh_burden_present ~ wave_label + wave_label:cluster_aligned +
        sex + pubertal_stage + wear_time_fraction + (1 | participant_id),
      data = master, REML = TRUE
    )
    out$burden <- .fixef_table(m_burden)
  }

  if ("mh_any_depression" %in% names(master)) {
    m_dep <- lme4::glmer(
      mh_any_depression ~ wave_label + wave_label:cluster_aligned +
        sex + pubertal_stage + wear_time_fraction + (1 | participant_id),
      data = master, family = stats::binomial
    )
    out$depression <- .fixef_table(m_dep)
  }

  out
}

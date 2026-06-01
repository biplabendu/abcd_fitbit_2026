# ===========================================================================
# Phase 6 (primary) — Feature-based outcome models
#
# These are the PRIMARY inferential models. Unlike the per-wave cluster labels,
# the behavioral features have a consistent meaning at every wave, so a
# `feature * wave_age` interaction is well-defined: it tests whether the
# association of a behavior with the outcome changes across adolescence
# (developmental moderation). Within-wave z-scored features are used so the
# predictor is "relative position among same-age peers" at both waves.
#
# Fixed effects: feature * wave_age + sex + pubertal_stage + wear_time_fraction
# Random effect: participant intercept (repeated measures).
# ===========================================================================

.fixef_tbl <- function(model) {
  tibble::as_tibble(summary(model)$coefficients, rownames = "term")
}

# Fit one feature × age model per feature for a single (continuous) outcome.
# Returns a named list of tidy fixed-effects tables.
fit_feature_outcome_models <- function(master, outcome, cfg = pipeline_config()) {
  if (!outcome %in% names(master)) {
    message("fit_feature_outcome_models: outcome '", outcome,
            "' not in master — skipping.")
    return(list())
  }
  feats <- intersect(cfg$primary_features, names(master))
  if (length(feats) == 0) {
    message("fit_feature_outcome_models: no primary features present — skipping.")
    return(list())
  }

  fit_one <- function(f) {
    m <- lmerTest::lmer(
      stats::reformulate(
        c(sprintf("%s * wave_age", f), "sex", "pubertal_stage",
          "wear_time_fraction", "(1 | participant_id)"),
        response = outcome
      ),
      data = master, REML = TRUE
    )
    .fixef_tbl(m)
  }

  stats::setNames(lapply(feats, fit_one), feats)
}

# Compact one-row-per-feature summary: the feature main effect and the
# feature × age interaction (estimate + p), with Benjamini–Hochberg FDR across
# features applied to the interaction p-values.
summarise_feature_models <- function(models) {
  if (length(models) == 0) return(tibble::tibble())

  pick <- function(tbl, pattern, col) {
    row <- grep(pattern, tbl$term)
    if (length(row) == 0) return(NA_real_)
    tbl[[col]][row[1]]
  }
  pcol <- function(tbl) if ("Pr(>|t|)" %in% names(tbl)) "Pr(>|t|)" else "Pr(>|z|)"

  out <- purrr::imap_dfr(models, function(tbl, feat) {
    pc <- pcol(tbl)
    tibble::tibble(
      feature     = feat,
      beta_main   = pick(tbl, paste0("^", feat, "$"), "Estimate"),
      se_main     = pick(tbl, paste0("^", feat, "$"), "Std. Error"),
      p_main      = pick(tbl, paste0("^", feat, "$"), pc),
      beta_x_age  = pick(tbl, ":wave_age$", "Estimate"),
      se_x_age    = pick(tbl, ":wave_age$", "Std. Error"),
      p_x_age     = pick(tbl, ":wave_age$", pc)
    )
  })
  out$p_x_age_fdr <- stats::p.adjust(out$p_x_age, method = "BH")
  out
}

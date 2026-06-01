# ===========================================================================
# Phase 6.3 — Multivariate outcome analysis (plan §6.3)
# Tests whether behavioral features jointly predict the outcome set:
#   - canonical correlation between z-scored features and outcomes
#   - MANOVA with cluster as grouping variable across all outcomes
# Operates on complete cases over the chosen feature/outcome columns.
#
# Uses base stats::cancor() (no extra dependencies) rather than the CCA
# package, which pulls in fields/maps and is awkward to build from source.
# ===========================================================================

.outcome_columns <- function(master, cfg) {
  c("mh_burden_present",
    intersect(cfg$metabolic_outcome_vars, names(master))) |>
    intersect(names(master))
}

multivariate_analysis <- function(master, cfg = pipeline_config()) {
  z_cols   <- grep("^z_", names(master), value = TRUE)
  out_cols <- .outcome_columns(master, cfg)

  if (length(z_cols) < 2 || length(out_cols) < 2) {
    message("multivariate_analysis: need >=2 feature and >=2 outcome columns ",
            "(have ", length(z_cols), " / ", length(out_cols), ") — skipping.")
    return(NULL)
  }

  cc_dat <- master |>
    dplyr::select(dplyr::all_of(c(z_cols, out_cols, "cluster_aligned"))) |>
    tidyr::drop_na()

  X <- as.matrix(cc_dat[z_cols])
  Y <- as.matrix(cc_dat[out_cols])

  # Canonical correlation (base R). Canonical variates are the centred data
  # projected onto the canonical coefficients.
  cc <- stats::cancor(X, Y)
  xscores <- sweep(X, 2, cc$xcenter) %*% cc$xcoef
  yscores <- sweep(Y, 2, cc$ycenter) %*% cc$ycoef

  manova_fit <- stats::manova(Y ~ cc_dat$cluster_aligned)

  list(
    canonical_cor = cc$cor,
    xcoef         = cc$xcoef,
    ycoef         = cc$ycoef,
    xscores       = xscores,
    yscores       = yscores,
    manova        = summary(manova_fit),
    n_complete    = nrow(cc_dat)
  )
}

# Fig 6 — canonical variates scatter (first pair).
plot_canonical_variates <- function(mv_result) {
  if (is.null(mv_result)) {
    message("plot_canonical_variates: no multivariate result — skipping.")
    return(NULL)
  }
  scores <- tibble::tibble(
    u1 = mv_result$xscores[, 1],
    v1 = mv_result$yscores[, 1]
  )
  ggplot2::ggplot(scores, ggplot2::aes(u1, v1)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    ggplot2::labs(x = "Behavioral canonical variate 1",
                  y = "Outcome canonical variate 1",
                  title = sprintf("Canonical correlation = %.3f",
                                  mv_result$canonical_cor[1])) +
    ggplot2::theme_linedraw(14)
}

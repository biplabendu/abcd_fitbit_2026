# ===========================================================================
# Phase 5 — Longitudinal trajectory analysis (plan §5)
# Works on participants with a valid cluster assignment at every active wave.
# With 2 waves the trajectory is yr2 -> yr6; extends automatically to 3 waves.
# ===========================================================================

# Participants present (clustered) at all active waves, wide by wave.
.complete_trajectories <- function(clusters, cfg = pipeline_config()) {
  waves <- unname(cfg$session_labels[cfg$sessions])
  wide <- clusters |>
    dplyr::select(participant_id, wave_label, cluster_aligned) |>
    tidyr::pivot_wider(names_from = wave_label, values_from = cluster_aligned)
  tidyr::drop_na(wide, dplyr::all_of(waves))
}

# Table 3 — cluster transition matrix + test vs. independence (first->last wave).
transition_matrix <- function(clusters, cfg = pipeline_config()) {
  waves <- unname(cfg$session_labels[cfg$sessions])
  wide  <- .complete_trajectories(clusters, cfg)
  from  <- wide[[waves[1]]]
  to    <- wide[[waves[length(waves)]]]

  counts <- table(from = from, to = to)
  list(
    counts       = as.data.frame.matrix(counts),
    probs        = as.data.frame.matrix(prop.table(counts, margin = 1)),
    chisq        = suppressWarnings(stats::chisq.test(counts)),
    n_complete   = nrow(wide)
  )
}

# Fig 3 — alluvial plot of cluster flows across waves.
plot_cluster_alluvial <- function(clusters, cfg = pipeline_config()) {
  waves <- unname(cfg$session_labels[cfg$sessions])
  long <- .complete_trajectories(clusters, cfg) |>
    tibble::rowid_to_column("traj_id") |>
    tidyr::pivot_longer(dplyr::all_of(waves),
                        names_to = "wave", values_to = "cluster") |>
    dplyr::mutate(wave = factor(wave, levels = waves))

  ggplot2::ggplot(long,
    ggplot2::aes(x = wave, stratum = cluster, alluvium = traj_id,
                 fill = cluster, label = cluster)) +
    ggalluvial::geom_flow(alpha = 0.6) +
    ggalluvial::geom_stratum() +
    ggplot2::labs(x = NULL, y = "Participants", fill = "Cluster") +
    ggplot2::theme_linedraw(14)
}

# §5.3 — feature-level change scores between consecutive waves.
delta_features <- function(features, cfg = pipeline_config()) {
  waves <- unname(cfg$session_labels[cfg$sessions])
  key   <- c("steps_total_mean_daily", "min_slp_duration_mean",
             "min_slp_onset_sd", "mets_mean_daily")
  key   <- intersect(key, names(features))

  wide <- features |>
    dplyr::select(participant_id, wave_label, dplyr::all_of(key)) |>
    tidyr::pivot_wider(names_from = wave_label, values_from = dplyr::all_of(key),
                       names_glue = "{.value}__{wave_label}")

  first <- waves[1]; last <- waves[length(waves)]
  for (k in key) {
    a <- paste0(k, "__", last); b <- paste0(k, "__", first)
    if (all(c(a, b) %in% names(wide))) {
      wide[[paste0("delta_", k, "_", first, "_", last)]] <- wide[[a]] - wide[[b]]
    }
  }
  dplyr::select(wide, participant_id, dplyr::starts_with("delta_"))
}

# §5.2 — group-based trajectory modeling on a key continuous feature.
# Fits lcmm::hlme for ng = 2..max_groups; returns models + BIC comparison.
fit_gbtm <- function(master, feature = "min_slp_onset_sd", max_groups = 4) {
  if (!feature %in% names(master)) {
    message("fit_gbtm: feature '", feature, "' not in master — skipping.")
    return(NULL)
  }
  dat <- master |>
    dplyr::transmute(
      participant_id,
      participant_id_num = as.integer(factor(participant_id)),
      wave_label,
      wave_age,
      y = .data[[feature]]
    ) |>
    tidyr::drop_na()

  # ng = 1 is fitted first; its parameter estimates seed the multi-group models
  # via the B argument (required by lcmm when ng > 1).
  m1 <- lcmm::hlme(y ~ wave_age, subject = "participant_id_num", ng = 1, data = dat)

  models <- c(
    list(m1),
    lapply(seq(2, max_groups), function(ng) {
      lcmm::hlme(y ~ wave_age, subject = "participant_id_num", ng = ng,
                 mixture = ~wave_age, B = m1, data = dat)
    })
  )

  bic      <- vapply(models, function(m) m$BIC, numeric(1))
  best_pos <- which.min(bic)            # 1-based position in models list
  best_ng  <- seq_len(max_groups)[best_pos]  # ng value (equals position by construction, but explicit)
  best     <- models[[best_pos]]

  # Characterize each latent class empirically: join the posterior class
  # assignment back to the data, then summarise group size, level, the
  # per-wave means, and the slope of the feature on age.
  dat_cls <- dplyr::left_join(
    dat, best$pprob[, c("participant_id_num", "class")],
    by = "participant_id_num"
  )

  class_summary <- dat_cls |>
    dplyr::group_by(class) |>
    dplyr::summarise(
      n_subj         = dplyr::n_distinct(participant_id_num),
      mean_level     = mean(y, na.rm = TRUE),
      slope_per_year = tryCatch(
        stats::coef(stats::lm(y ~ wave_age))[["wave_age"]],
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct = round(100 * n_subj / sum(n_subj), 1))

  wave_means <- dat_cls |>
    dplyr::group_by(class, wave_label) |>
    dplyr::summarise(mean_y = mean(y, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = wave_label, values_from = mean_y)

  class_summary <- dplyr::left_join(class_summary, wave_means, by = "class")

  # --- Assign interpretable roles to the latent classes (k = 3 case) --------
  # regular  = lowest mean level; worsening = steepest positive slope of the
  # remaining classes; irregular = the other. Falls back to "class{n}" labels
  # for other k so downstream code degrades gracefully.
  roles <- setNames(paste0("class", class_summary$class), class_summary$class)
  if (nrow(class_summary) == 3) {
    reg_cl  <- class_summary$class[which.min(class_summary$mean_level)]
    rest    <- class_summary[class_summary$class != reg_cl, ]
    wors_cl <- rest$class[which.max(rest$slope_per_year)]
    irr_cl  <- rest$class[rest$class != wors_cl]
    roles[as.character(reg_cl)]  <- "regular"
    roles[as.character(wors_cl)] <- "worsening"
    roles[as.character(irr_cl)]  <- "irregular"
  }
  class_summary$role <- roles[as.character(class_summary$class)]

  # --- Participant-keyed membership with role-named posterior probabilities --
  id_map <- dplyr::distinct(dat, participant_id, participant_id_num)
  pp     <- dplyr::left_join(best$pprob, id_map, by = "participant_id_num")
  prob_cols <- grep("^prob", names(pp), value = TRUE)  # prob1, prob2, ...

  membership <- tibble::tibble(
    participant_id = pp$participant_id,
    traj_group     = roles[as.character(pp$class)],
    max_posterior  = apply(pp[prob_cols], 1, max)
  )
  # rename probK -> prob_<role>
  for (k in class_summary$class) {
    membership[[paste0("prob_", roles[as.character(k)])]] <- pp[[paste0("prob", k)]]
  }

  quality <- tibble::tibble(
    mean_max_posterior = mean(membership$max_posterior, na.rm = TRUE),
    pct_confident_0.8  = round(100 * mean(membership$max_posterior > 0.8, na.rm = TRUE), 1)
  )

  list(models     = models,
       bic        = tibble::tibble(ng = seq_len(max_groups), BIC = bic),
       best_ng    = best_ng,
       feature    = feature,
       class_summary = class_summary,
       membership = membership,
       quality    = quality)
}

# --- Trajectory-group outcome model ---------------------------------------
# Uses the GBTM soft membership probabilities (reference = "regular") as
# person-level predictors of an outcome, interacted with age. Probabilities
# (not hard labels) are used to respect classification uncertainty.
fit_trajectory_outcome <- function(master, gbtm, outcome = "mh_burden_present") {
  if (is.null(gbtm) || is.null(gbtm$membership)) {
    message("fit_trajectory_outcome: no GBTM membership — skipping.")
    return(NULL)
  }
  if (!outcome %in% names(master)) {
    message("fit_trajectory_outcome: outcome '", outcome, "' absent — skipping.")
    return(NULL)
  }

  dat <- dplyr::left_join(master, gbtm$membership, by = "participant_id")

  # Non-reference soft-probability predictors. Take ONLY the trajectory
  # membership columns (not e.g. the GMM `prob_cluster_*` columns in master),
  # and drop the reference group (prob_regular).
  traj_prob_cols <- setdiff(names(gbtm$membership),
                            c("participant_id", "traj_group", "max_posterior"))
  prob_preds <- setdiff(traj_prob_cols, "prob_regular")
  if (length(prob_preds) == 0) {
    message("fit_trajectory_outcome: no non-reference probability columns — skipping.")
    return(NULL)
  }

  rhs <- c(paste0("(", paste(prob_preds, collapse = " + "), ") * wave_age"),
           "sex", "pubertal_stage", "wear_time_fraction", "(1 | participant_id)")
  m <- lmerTest::lmer(stats::reformulate(rhs, response = outcome),
                      data = dat, REML = TRUE)
  list(table = tibble::as_tibble(summary(m)$coefficients, rownames = "term"),
       predictors = prob_preds,
       n_obs = nrow(stats::model.frame(m)))
}

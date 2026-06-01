# ===========================================================================
# Phase 4 — Cluster characterization (plan §4.3)
# Tables and figures describing each cluster's behavioral profile, by wave.
# Cluster names are assigned descriptively from profiles BEFORE outcome models.
# ===========================================================================

# Table 2 — mean ± SD of each (raw) feature by cluster × wave.
cluster_profiles_table <- function(master) {
  feat_cols <- master |>
    dplyr::select(dplyr::matches("_(mean_daily|sd_daily|cv|amplitude|acrophase|IS|IV|L5|M10|relative_amplitude|sedentary_fraction|mvpa_fraction|duration_mean|onset_mean|onset_sd|efficiency_mean|social_jet_lag)$")) |>
    dplyr::select(-dplyr::starts_with("z_")) |>
    names()

  master |>
    dplyr::group_by(wave_label, cluster_aligned) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(feat_cols),
                    list(mean = ~mean(.x, na.rm = TRUE),
                         sd   = ~stats::sd(.x, na.rm = TRUE)),
                    .names = "{.col}__{.fn}"),
      .groups = "drop"
    )
}

# Table 1 — participant characteristics by cluster × wave.
participant_characteristics_table <- function(master) {
  master |>
    dplyr::group_by(wave_label, cluster_aligned) |>
    dplyr::summarise(
      n              = dplyr::n(),
      pct_female     = mean(sex == "Female", na.rm = TRUE) * 100,
      mean_age       = mean(wave_age, na.rm = TRUE),
      mean_puberty   = mean(pubertal_stage, na.rm = TRUE),
      mean_wear_frac = mean(wear_time_fraction, na.rm = TRUE),
      .groups = "drop"
    )
}

# Fig 1 — UMAP scatter coloured by cluster, faceted by wave.
plot_umap_clusters <- function(umap, clusters) {
  df <- umap |>
    dplyr::left_join(
      dplyr::select(clusters, participant_id, session_id, cluster_aligned),
      by = c("participant_id", "session_id")
    )

  ggplot2::ggplot(df, ggplot2::aes(umap_1, umap_2, colour = cluster_aligned)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggplot2::facet_wrap(~wave_label) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") +
    ggplot2::theme_linedraw(14) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3)))
}

# Supp Table 1 — wear quality by cluster (companion to the confound check).
wear_by_cluster_table <- function(master) {
  master |>
    dplyr::group_by(wave_label, cluster_aligned) |>
    dplyr::summarise(
      n                   = dplyr::n(),
      mean_n_valid_days   = mean(n_valid_days, na.rm = TRUE),
      mean_wear_fraction  = mean(wear_time_fraction, na.rm = TRUE),
      mean_dispersion     = mean(valid_day_dispersion, na.rm = TRUE),
      .groups = "drop"
    )
}

# Fig 2 — mean diurnal curves per cluster.
# Reads the slot-level handoff data/handoff/diurnal_profiles.csv exported by
# python/04_feature_extraction.py (participant_id, session_id, wave_label,
# signal, start_hour, value). Returns NULL with a message if absent.
plot_diurnal_curves <- function(clusters, cfg = pipeline_config()) {
  path <- file.path(cfg$handoff_dir, "diurnal_profiles.csv")
  if (!file.exists(path)) {
    message("plot_diurnal_curves: ", path, " not found — skipping (see TODO).")
    return(NULL)
  }
  prof <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::left_join(
      dplyr::select(clusters, participant_id, session_id, cluster_aligned),
      by = c("participant_id", "session_id")
    ) |>
    dplyr::group_by(wave_label = session_to_wave(session_id, cfg),
                    cluster_aligned, signal, start_hour) |>
    dplyr::summarise(mean = mean(value, na.rm = TRUE),
                     se   = stats::sd(value, na.rm = TRUE) / sqrt(dplyr::n_distinct(participant_id)),
                     .groups = "drop")

  ggplot2::ggplot(prof, ggplot2::aes(start_hour, mean, colour = cluster_aligned,
                                     fill = cluster_aligned)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                         alpha = 0.2, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(signal ~ wave_label, scales = "free_y") +
    ggplot2::labs(x = "Hour of day", y = "Mean signal", colour = "Cluster", fill = "Cluster") +
    ggplot2::theme_linedraw(14)
}

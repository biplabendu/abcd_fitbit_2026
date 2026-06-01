# ===========================================================================
# Compliance confound check (plan §4.2)
# Before any scientific analysis: verify clusters are not driven by wear
# compliance. Kruskal-Wallis test of wear_time_fraction (and n_valid_days)
# across cluster_aligned, run per wave. A significant result means wear
# compliance is a confound and wear_time_fraction must enter all models.
# ===========================================================================

compliance_confound_check <- function(master) {
  master |>
    dplyr::group_by(wave_label) |>
    dplyr::group_modify(~ {
      kw_wear <- stats::kruskal.test(wear_time_fraction ~ cluster_aligned, data = .x)
      kw_days <- stats::kruskal.test(n_valid_days       ~ cluster_aligned, data = .x)
      tibble::tibble(
        metric    = c("wear_time_fraction", "n_valid_days"),
        statistic = c(kw_wear$statistic, kw_days$statistic),
        df        = c(kw_wear$parameter, kw_days$parameter),
        p_value   = c(kw_wear$p.value, kw_days$p.value)
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(confound_flag = p_value < 0.05)
}

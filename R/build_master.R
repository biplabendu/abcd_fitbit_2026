# ===========================================================================
# Build the master analysis frame
# Joins all handoff + outcome tables on participant_id + session_id.
# cluster_assignments is the spine (one row per participant × wave with a
# cluster label); everything else is left-joined onto it.
# ===========================================================================

build_master <- function(clusters, gmm, features, wear,
                         mh, metabolic, demographics) {
  join_keys <- c("participant_id", "session_id", "wave_label")

  # features and wear both contain n_valid_days (same values, different source).
  # Drop it from wear before joining to avoid .x / .y suffixes in master.
  wear_clean <- dplyr::select(wear, -dplyr::any_of("n_valid_days"))

  master <- clusters |>
    dplyr::left_join(gmm,          by = join_keys) |>
    dplyr::left_join(features,     by = join_keys) |>
    dplyr::left_join(wear_clean,   by = join_keys) |>
    dplyr::left_join(mh,           by = join_keys) |>
    dplyr::left_join(metabolic,    by = join_keys) |>
    dplyr::left_join(demographics, by = join_keys)

  # Basic integrity checks
  stopifnot(!any(duplicated(master[join_keys])))
  master
}

# ===========================================================================
# Load Python handoff CSVs  (data/handoff/)
# These functions are the R side of the Python→R contract. Each asserts the
# columns it depends on before returning, per the plan's guidance.
# ===========================================================================

.assert_cols <- function(df, cols, what) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("[%s] missing expected columns: %s",
                 what, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(df)
}

# Features: one CSV per wave (features_<label>.csv); bind into one frame.
load_features <- function(feature_files) {
  df <- purrr::map_dfr(feature_files, readr::read_csv, show_col_types = FALSE)
  .assert_cols(df, c("participant_id", "session_id", "wave_label", "n_valid_days"),
               "features")
  df
}

load_clusters <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("participant_id", "session_id", "wave_label",
                     "cluster_hard", "cluster_aligned"), "cluster_assignments")
  df |>
    dplyr::mutate(cluster_aligned = factor(cluster_aligned))
}

load_gmm <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("participant_id", "session_id", "wave_label"),
               "gmm_probabilities")
  df
}

load_wear <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("participant_id", "session_id", "wave_label",
                     "n_valid_days", "wear_time_fraction",
                     "valid_day_dispersion", "exclude_clustering"),
               "wear_quality")
  df
}

load_umap <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("participant_id", "session_id", "wave_label",
                     "umap_1", "umap_2"), "umap_coords")
  df
}

load_centroids <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("cluster", "session_id", "wave_label"), "cluster_centroids")
  df
}

load_alignment <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  .assert_cols(df, c("ref_wave", "target_wave", "aligned_cluster",
                     "cosine_similarity"), "cluster_alignment")
  df
}

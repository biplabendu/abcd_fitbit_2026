# ===========================================================================
# _targets.R — R analysis pipeline (plan Phases 4–7)
#
# Consumes the Python handoff CSVs (data/handoff/) plus outcome/covariate data
# (via NBDCtools) and runs: load → master → QC → characterization →
# longitudinal → outcome models → figures/tables.
#
# Usage:
#   targets::tar_make()              # run the pipeline
#   targets::tar_visnetwork()        # view the DAG
#   targets::tar_manifest()          # list targets without running
#
# All configuration lives in R/config_pipeline.R (pipeline_config()).
# ===========================================================================

library(targets)
library(tarchetypes)

# Source all helper/function files in R/ (config, loaders, models, plots).
tar_source("R")

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2",
    "ggalluvial", "lme4", "lmerTest", "lcmm", "NBDCtools"
  ),
  format = "rds"
)

list(

  # ---- Configuration -------------------------------------------------------
  tar_target(
    name = cfg, 
    command = pipeline_config()
  ),
  # ---- File inputs (handoff CSVs) -----------------------------------------
  # `features_files` auto-extends when sessions are added in config.
  tar_target(
    name = features_files,
    command = file.path(
      cfg$handoff_dir,
      paste0("features_", unname(cfg$session_labels[cfg$sessions]), ".csv")
    ),
    format = "file"
  ),
  tar_target(
    name = clusters_file,  
    command = file.path(cfg$handoff_dir, "cluster_assignments.csv"), 
    format = "file"
  ),
  tar_target(
    name = gmm_file,
    command = file.path(cfg$handoff_dir, "gmm_probabilities.csv"),
    format = "file"
  ),
  tar_target(
    name = wear_file,
    command = file.path(cfg$handoff_dir, "wear_quality.csv"),
    format = "file"
  ),
  tar_target(
    name = umap_file,
    command = file.path(cfg$handoff_dir, "umap_coords.csv"),
    format = "file"
  ),
  tar_target(
    name = centroids_file,
    command = file.path(cfg$handoff_dir, "cluster_centroids.csv"),
    format = "file"
  ),
  tar_target(
    name = alignment_file,
    command = file.path(cfg$handoff_dir, "cluster_alignment.csv"),
    format = "file"
  ),

  # ---- Phase 4.1 Load ------------------------------------------------------
  tar_target(
    name = features,
    command = load_features(features_files)
  ),
  tar_target(
    name = clusters,
    command = load_clusters(clusters_file)
  ),
  tar_target(
    name = gmm,
    command = load_gmm(gmm_file)
  ),
  tar_target(
    name = wear,
    command = load_wear(wear_file)
  ),
  tar_target(
    name = umap,
    command = load_umap(umap_file)
  ),
  tar_target(
    name = centroids,
    command = load_centroids(centroids_file)
  ),
  tar_target(
    name = alignment,
    command = load_alignment(alignment_file)
  ),

  # Outcome / covariate data (NBDCtools reads a directory, not a tracked file)
  tar_target(
    name = mh,
    command = load_mh_outcomes(cfg)
  ),
  tar_target(
    name = metabolic,
    command = load_metabolic_outcomes(cfg)
  ),
  tar_target(
    name = demographics,
    command = load_demographics(cfg)
  ),

  # ---- Phase 4.1 Master frame ---------------------------------------------
  tar_target(
    name = master,
    command = build_master(
      clusters, gmm, features, wear, mh, metabolic, demographics
    )
  ),

  # ---- Phase 4.2 Compliance confound check --------------------------------
  tar_target(
    name = qc_compliance,
    command = compliance_confound_check(master)
  ),
  tar_target(
    name = qc_compliance_csv,
    command = save_table(qc_compliance, "supp1_compliance_confound", cfg),
    format = "file"
  ),

  # ---- Phase 4.3 Cluster characterization ---------------------------------
  tar_target(
    name = tbl_profiles,
    command = cluster_profiles_table(master)
  ),
  tar_target(
    name = tbl_profiles_csv,
    command = save_table(tbl_profiles, "table2_cluster_profiles", cfg),
    format = "file"
  ),
  tar_target(
    name = tbl_participants,
    command = participant_characteristics_table(master)
  ),
  tar_target(
    name = tbl_participants_csv,
    command = save_table(tbl_participants, "table1_participant_characteristics", cfg),
    format = "file"
  ),
  tar_target(
    name = tbl_wear_cluster,
    command = wear_by_cluster_table(master)
  ),
  tar_target(
    name = tbl_wear_cluster_csv,
    command = save_table(tbl_wear_cluster, "supp1b_wear_by_cluster", cfg),
    format = "file"
  ),
  tar_target(
    name = fig_umap,
    command = plot_umap_clusters(umap, clusters)
  ),
  tar_target(
    name = fig_umap_files,
    command = save_figure(fig_umap, "fig1_umap_clusters", cfg),
    format = "file"
  ),
  # Fig 2 — diurnal curves (NULL until diurnal_profiles.csv is exported)
  tar_target(
    name = fig_diurnal,
    command = plot_diurnal_curves(clusters, cfg)
  ),

  # ---- Phase 5 Longitudinal trajectories ----------------------------------
  tar_target(
    name = transitions,
    command = transition_matrix(clusters, cfg)
  ),
  tar_target(
    name = transitions_csv,
    command = save_table(
      tibble::rownames_to_column(transitions$counts, "from"),
      "table3_transition_counts", cfg
    ),
    format = "file"
  ),
  tar_target(
    name = fig_alluvial,
    command = plot_cluster_alluvial(clusters, cfg)
  ),
  tar_target(
    name = fig_alluvial_files,
    command = save_figure(fig_alluvial, "fig3_cluster_alluvial", cfg),
    format = "file"
  ),
  tar_target(
    name = deltas,
    command = delta_features(features, cfg)
  ),
  # §5.2 group-based trajectory modeling (NULL if feature absent)
  tar_target(
    name = gbtm,
    command = fit_gbtm(master)
  ),

  # ---- Phase 6 (primary) Feature-based outcome models ---------------------
  # PRIMARY models: behavioral feature × age (consistent meaning across waves).
  tar_target(
    name = mh_feature_models,
    command = fit_feature_outcome_models(master, "mh_burden_present", cfg)
  ),
  tar_target(
    name = mh_feature_summary,
    command = summarise_feature_models(mh_feature_models)
  ),

  # ---- Phase 6 (secondary) Trajectory-group outcome model -----------------
  # Longitudinally-coherent grouping (GBTM soft probabilities).
  tar_target(
    name = mh_traj_model,
    command = fit_trajectory_outcome(master, gbtm, "mh_burden_present")
  ),

  # ---- Phase 6 Cluster outcome models (within-wave, descriptive) ----------
  # Wave-nested cluster effects (no cross-wave pooling). See models_mh.R.
  tar_target(
    name = mh_models,
    command = fit_mh_models(master)
  ),
  tar_target(
    name = metabolic_models,
    command = fit_metabolic_models(master, cfg)
  ),

  # ---- Phase 6.3 Multivariate behavioral–health axis ----------------------
  tar_target(
    name = multivariate,
    command = multivariate_analysis(master, cfg)
  ),
  tar_target(
    name = fig_canonical,
    command = plot_canonical_variates(multivariate)
  )
)

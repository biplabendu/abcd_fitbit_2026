library(dplyr)
library(timecourseRnaseq)
library(pheatmap)
for (f in list.files("R", full.names = TRUE)) source(f)

dat_full <- readRDS(
  "dev/out/fitbit_median_for_clustering.RDS"
)


# Filter to participants with weekday and weekend -------------------------

ids_to_remove <- dat_full |> 
  group_by(id, sess) |> 
  reframe(
    wkday = if_else(sum(!is_wknd) > 0, "exist", "x"),
    wknd = if_else(sum(is_wknd) > 0, "exist", "x")
  ) |> 
  filter(
    wkday == "x" | wknd == "x"
  ) |> 
  pull(1)

dat_sub1 <- dat_full |> 
  glimpse() |> 
  filter_out(
    id %in% ids_to_remove
  )

# Handle duplicates -------------------------------------------------------

ids_to_remove2 <- dat_sub1 |> 
  count(id, sess, is_wknd, dt_hr) |> 
  filter(n > 1) |> 
  pull(1) |> 
  unique()

dat_sub2 <- dat_sub1 |> 
  filter_out(
    id %in% ids_to_remove2
  )


# Z-score -----------------------------------------------------------------

dat_list <- dat_sub2 |> 
  group_by(
    sess, is_wknd
  ) |> 
  mutate(
    across(
      c(hrate, min_slp, steps_total),
      ~ ((.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)) |> 
        round(2),
      .names = "z_{.col}"
    )
  ) |> 
  ungroup() |> 
  group_split(sess, is_wknd) |> 
  setNames(
    c("yr2-wkday", "yr2-wknd", "yr6-wkday", "yr6-wknd")
  )

# steps -------------------------------------------------------------------

## Year 2 -------------------------------
steps_yr2_wkday <- dat_list[["yr2-wkday"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wkday_T{dt_hr}"
    ),
    values = z_steps_total
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

steps_yr2_wknd <- dat_list[["yr2-wknd"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wknd_T{dt_hr}"
    ),
    values = z_steps_total
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

## Year 6 -------------------------------
steps_yr6_wkday <- dat_list[["yr6-wkday"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wkday_T{dt_hr}"
    ),
    values = z_steps_total
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

steps_yr6_wknd <- dat_list[["yr6-wknd"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wknd_T{dt_hr}"
    ),
    values = z_steps_total
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

# STEPS (Visualize) -------------------------------------------------------------

ids_steps <- intersect(
  intersect(steps_yr2_wkday$id, steps_yr2_wknd$id),
  intersect(steps_yr6_wkday$id, steps_yr6_wknd$id)
)

dat_steps <- steps_yr2_wkday |>
# steps_yr2_wkday |> 
  filter(
    id %in% ids_steps
  ) |> 
  left_join(steps_yr2_wknd, join_by(id)) |> 
  left_join(steps_yr6_wkday, join_by(id)) |> 
  left_join(steps_yr6_wknd, join_by(id)) |> 
  glimpse()

## ALL TOGETHER ------
dat_steps |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "STEPS: Yr 2 (Wkday + Wknds) + Yr 6"
  )

## BY EVENT ---------
dat_steps |> 
  select(id, matches("yr2")) |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "STEPS: Yr 2 (Wkday + Wknds)"
  )

dat_steps |> 
  select(id, matches("yr6")) |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "STEPS: Yr 6 (Wkday + Wknds)"
  )

## Save Data: STEPS ------
dat_steps |>
  saveRDS(
    "dev/out/zscore_steps.RDS"
  )


# sleep -------------------------------------------------------------------

## Year 2 -------------------------------
sleep_yr2_wkday <- dat_list[["yr2-wkday"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wkday_T{dt_hr}"
    ),
    values = z_min_slp
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

sleep_yr2_wknd <- dat_list[["yr2-wknd"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wknd_T{dt_hr}"
    ),
    values = z_min_slp
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

## Year 6 -------------------------------
sleep_yr6_wkday <- dat_list[["yr6-wkday"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wkday_T{dt_hr}"
    ),
    values = z_min_slp
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

sleep_yr6_wknd <- dat_list[["yr6-wknd"]] |>
  transmute(
    id,
    zt = glue::glue(
      "{sess}_wknd_T{dt_hr}"
    ),
    values = z_min_slp
  ) |> 
  tidyr::pivot_wider(
    names_from = zt,
    values_from = values
  ) |> 
  filter_out(
    if_any(
      everything(),
      ~ is.na(.x)
    )
  )

# SLEEP (Visualize) -------------------------------------------------------------

ids_sleep <- intersect(
  intersect(sleep_yr2_wkday$id, sleep_yr2_wknd$id),
  intersect(sleep_yr6_wkday$id, sleep_yr6_wknd$id)
)

dat_sleep <- sleep_yr2_wkday |>
  filter(
    id %in% ids_sleep
  ) |> 
  left_join(sleep_yr2_wknd, join_by(id)) |> 
  left_join(sleep_yr6_wkday, join_by(id)) |> 
  left_join(sleep_yr6_wknd, join_by(id)) |> 
  glimpse()

## ALL TOGETHER ------
dat_sleep |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "STEPS: Yr 2 (Wkday + Wknds) + Yr 6"
  )

## BY EVENT ---------
dat_sleep |> 
  select(id, matches("yr2")) |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "SLEEP: Yr 2 (Wkday + Wknds)"
  )

dat_sleep |> 
  select(id, matches("yr6")) |> 
  tc_plot_heatmap(
    show_rownames = FALSE,
    n_clusters = 4,
    title = "SLEEP: Yr 6 (Wkday + Wknds)"
  )

## Save Data: SLEEP ------
dat_sleep |>
  saveRDS(
    "dev/out/zscore_sleep.RDS"
  )

# heart rate -------------------------------------------------------------------

# NOTE: The data isn't very useful for within-day variability
#         Check plot for an overview.



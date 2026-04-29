library(dplyr)
library(ggplot2)
library(lubridate)

# replace with your directory path
# folder containing the pre-assembled dataset
dir_abcd <- "dev/data/7_0/0_1/data"


# Vars of interest --------------------------------------------------------
vars <- read.csv(
  "data/vars-of-interest/psychiatric_disorders-diagnosis-sim04.csv"
) |> 
  as_tibble()

data <- NBDCtools::create_dataset_abcd(
  dir = dir_abcd,
  vars = vars$name,
  vars_add = c(
    "ab_g_dyn__design_site",
    "ab_g_dyn__visit_age"
  )
)

# Fitbit data -------------------------------------------------------------
dat <- arrow::read_parquet(
  "dev/data/fitbit-summaries/activity_120m.parquet"
) |> 
  as_tibble() |> 
  arrange(
    participant_id,
    session_id,
    day,
    dt,
    start
  ) |> 
  mutate(
    dt_hr = hour(start) + 1
  ) 


# Summarize data ----------------------------------------------------------
summ <- dat |> 
  group_by(
    participant_id,
    session_id
  ) |> 
  reframe(
    n_rows = n(),
    n_days = length(unique(day)),
    has_wknds = if_else(
      sum(dt_wknd, na.rm = TRUE) > 0,
      TRUE,
      FALSE
    )
  ) |> 
  group_by(participant_id) |> 
  mutate(
    n_sess = length(session_id)
  ) |> 
  ungroup()

## Visualize summary -----
summ |> 
  select(
    participant_id,
    n_sess
  ) |> 
  distinct() |> 
  count(n_sess) |> 
  ggplot(
    aes(
      x = n_sess,
      y = n
    )
  ) +
  geom_bar(
    stat = "identity",
    fill = "maroon"
  ) +
  labs(
    x = "N (Sessions)",
    y = "N (Participants)"
  ) +
  theme_linedraw(24)

# Identify sessions of interest -----------------------------------------

summ_daily <- dat |> 
  select(
    participant_id,
    session_id,
    day,
    dt,
    dt_day,
    dt_wknd
  ) |> 
  distinct() |> 
  group_by(
    participant_id,
    session_id
  ) |> 
  reframe(
    # n_rows = n(),
    # n_days = length(unique(day)),
    # n_wkdays = sum(!dt_wknd),
    n_wknds = sum(dt_wknd)
  ) |> 
  mutate(
    session_id = factor(
      session_id,
      levels = c(
        "ses-02A",
        "ses-04A",
        "ses-06A",
        "ses-08A"
      )
    )
  )

summ_daily |> 
  tidyr::pivot_wider(
    names_from = session_id,
    values_from = n_wknds
  ) |> 
  mutate(
    across(
      !participant_id,
      ~ case_when(
        .x >= 2 ~ "2+wknds",
        .default = "-"
      )
    )
  ) |> 
  count(
    `ses-02A`,
    `ses-06A`
  ) |> 
  arrange(
    desc(`ses-02A`),
    desc(`ses-06A`)
  )


# Identify IDs ------------------------------------------------------------
ids_v02_v06 <- summ_daily |> 
  tidyr::pivot_wider(
    names_from = session_id,
    values_from = n_wknds
  ) |> 
  mutate(
    across(
      !participant_id,
      ~ case_when(
        .x >= 2 ~ "keep",
        .default = "-"
      )
    )
  ) |> 
  filter(
    `ses-02A` == "keep" &
      `ses-06A` == "keep"
  )

ids_v02_v06 |> 
  write.csv(
    "data/ids_fitbit_v02_v06.csv",
    row.names = FALSE
  )

# Subset / Format data ----------------------------------------------------

sdat <- dat |> 
  select(
    id = participant_id,
    sess = session_id,
    day,
    dt,
    dt_hr,
    is_wknd = dt_wknd,
    hrate = hrate_rest_fitb,
    min_slp,
    steps_total
  ) |> 
  filter(
    id %in% ids_v02_v06$participant_id,
    sess %in% c("ses-02A", "ses-06A")
  ) |> 
  mutate(
    sess = if_else(sess == "ses-02A", "yr2", "yr6"),
    yr = year(dt),
    mnth = month(dt)
  )

out <- sdat |> 
  group_by(id, sess, yr, is_wknd, dt_hr) |> 
  reframe(
    mnths = paste(unique(mnth), collapse = ", "),
    hrate = median(hrate, na.rm = TRUE),
    min_slp = median(min_slp, na.rm = TRUE),
    steps_total = median(steps_total, na.rm = TRUE)
  )

# out |> 
#   saveRDS(
#     "dev/out/fitbit_median_for_clustering.RDS"
#   )

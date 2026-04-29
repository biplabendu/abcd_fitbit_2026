# Load functions / pkgs ---------------------------------------------------
library(dplyr)
library(ggplot2)
library(lubridate)
for (f in list.files("R", full.names = TRUE)) source(f)


# Which IDs? --------------------------------------------------------------

# replace with your directory path
# folder containing the pre-assembled dataset
dir_abcd <- "dev/data/fitbit-summaries/7_0/0_1/data"

# load IDs for which we have fitbit data
ids_steps <- readRDS(
  "dev/out/zscore_steps.RDS"
) |> 
  pull(id)
ids_sleep <- readRDS(
  "dev/out/zscore_sleep.RDS"
) |> 
  pull(id)

# symdiff(ids_sleep, ids_steps)
# Same participants have both sleep and step data!!
# n = 2,919

ids <- intersect(ids_sleep, ids_steps)
# write.csv(
#   data.frame(
#     participant_id = ids
#   ),
#   "data/ids_fitbit_v02_v06-sleep_steps_data.csv",
#   row.names = FALSE
# )

## Create dataset --------------------------------------------------------
vars <- read.csv(
  "data/vars-of-interest/psychiatric_disorders-diagnosis-sim04.csv"
) |> 
  as_tibble()

data_meta <- NBDCtools::create_dataset_abcd(
  dir = dir_abcd,
  vars = vars$name,
  vars_add = c(
    "ab_g_dyn__design_site",
    "ab_g_dyn__visit_age",
    "ab_g_stc__cohort_sex",
    "ab_g_stc__cohort_ethn",
    "ab_g_stc__cohort_race__nih",
    "ab_g_dyn__cohort_income__hhold__3lvl",
    "ab_g_dyn__cohort_edu__cgs"
  )
) |> 
  format_sociodemo_vars()

# Cohort composition ------------------------------------------

## Prep data ---------
visit_yrs <- c("2-yr", "4-yr", "6-yr")

out_demo <- data_meta |>
  mutate(
    visit = case_when(
      session_id == "ses-02A" ~ "2-yr",
      session_id == "ses-04A" ~ "4-yr",
      session_id == "ses-06A" ~ "6-yr",
      .default = NA_character_
    ) |> 
      factor(
        levels = visit_yrs
      ),
    .after = session_id
  ) |> 
  filter(
    participant_id %in% ids,
    visit %in% visit_yrs
  ) |> 
  select(
    participant_id,
    visit,
    sex,
    ethn,
    race_nih,
    hhold_3lvl,
    edu
  ) |> 
  glimpse()

## Visualize data ---------
create_alluvial_plot(
  data = out_demo,
  visit_yrs = visit_yrs,
  status_col = "sex",
  other_label = "no-fitbit-data",
  legend_title = "Sex"
) +
  theme(
    legend.position = "right"
  )
create_alluvial_plot(
  data = out_demo,
  visit_yrs = visit_yrs,
  status_col = "ethn",
  other_label = "no-fitbit-data",
  legend_title = "Ethnicity"
) +
  theme(
    legend.position = "right"
  )
create_alluvial_plot(
  data = out_demo,
  visit_yrs = visit_yrs,
  status_col = "race_nih",
  other_label = "no-fitbit-data",
  legend_title = "Race (NIH)",
  color_pallete = "Paired"
) +
  theme(
    legend.position = "right"
  )
create_alluvial_plot(
  data = out_demo,
  visit_yrs = visit_yrs,
  other_label = "no-fitbit-data",
  status_col = "hhold_3lvl",
  legend_title = "Household income",
  color_pallete = "Paired"
) +
  theme(
    legend.position = "right"
  )
create_alluvial_plot(
  data = out_demo,
  visit_yrs = visit_yrs,
  other_label = "no-fitbit-data",
  status_col = "edu",
  legend_title = "Caregiver education",
  color_pallete = "Paired"
) +
  theme(
    legend.position = "right"
  )

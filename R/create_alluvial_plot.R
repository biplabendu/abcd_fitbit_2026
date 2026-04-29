create_alluvial_plot <- function(data, 
                                 type = NULL, 
                                 visit_yrs = NULL, 
                                 color_pallete = "Blues",
                                 status_col = "status_compl",
                                 other_label = "Other",
                                 legend_title = "Status") {
  # ###-###-###-###-###-###-###-###-###-###-
  # # DEV ----
  # data = plot_data_compl
  # type = "t1"
  # visit_yrs = visit_yrs
  # status_col = "ethn"
  # ###-###-###-###-###-###-###-###-###-###-
  
  # status_levels_plot <- c(
  #   "Complete", 
  #   "Incomplete", 
  #   "Not administered",
  #   "Other"
  # )
  status_levels_plot <- c(
    data |> pull(!!status_col) |> levels(),
    other_label
  )
  
  if (length(status_levels_plot) == 4) {
    status_colors <- c(
      "#0072B2",
      "#F19221",
      "#F01717",
      "grey70"
    )
  } else {
    n_lvls = length(status_levels_plot)
    status_colors <- c(
      "grey70",
      RColorBrewer::brewer.pal(n_lvls - 1, color_pallete)
    ) |> 
      rev()
  }
  
  # Filter Data 
  if (!is.null(type)) {
    data <- data |> 
      filter(
        what == type
      )
  }
  if (!is.null(visit_yrs)) {
    data <- data |> 
      filter(
        visit %in% visit_yrs
      )
  }
  
  df_long <- data |>
    tidyr::pivot_wider(
      id_cols = participant_id,
      names_from = "visit",
      values_from = status_col
    ) |>
    ggsankey::make_long(
      visit_yrs
    ) |> 
    # create other category
    mutate(
      across(
        c(node, next_node),
        ~ if_else(
          is.na(.x),
          other_label,
          .x
        )
      ),
      next_node = if_else(
        is.na(next_x),
        NA,
        next_node
      )
    )
  
  # --- Plotting with geom_alluvial ---
  ggplot(
    df_long, 
    aes(
      x = x, 
      next_x = next_x, 
      node = node, 
      next_node = next_node,
      fill = node, 
      label = node
    )
  ) +
    # geom_hline(
    #   yintercept = c(0, 5000, 10000),
    #   col = "grey30",
    #   lty = 1,
    #   lwd = 1,
    #   alpha = 0.5
    # ) +
    # Create the flows
    ggsankey::geom_alluvial(
      flow.alpha = 0.5, 
      width = 0.5
    ) +
    # ggsankey::geom_alluvial_text(
    #   size = 4, 
    #   color = "white", 
    #   fontface = "bold"
    # ) +
    scale_fill_manual(
      name = legend_title,
      values = status_colors,
      na.value = NA,
      breaks = status_levels_plot
    ) +
    # Use the built-in alluvial theme
    ggsankey::theme_alluvial(base_size = 20) +
    # Set labels and title
    labs(
      y = "Participants", 
      x = NULL, 
      title = glue::glue("{toupper(type)}")
    ) +
    # Additional theme tweaks
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

format_sociodemo_vars <- function(data) {
  data |> 
    mutate(
      ethn = ab_g_stc__cohort_ethn,
      race_nih = ab_g_stc__cohort_race__nih,
      hhold_3lvl = ab_g_dyn__cohort_income__hhold__3lvl,
      edu = ab_g_dyn__cohort_edu__cgs,
      sex = ab_g_stc__cohort_sex
    ) |> 
    mutate(
      sex = case_match(
        sex,
        "1" ~ "Male",
        "2" ~ "Female",
        .default = "Unknown"
      ) |> 
        factor(
          levels = c(
            "Male",
            "Female",
            "Unknown"
          )
        ),
      ethn = case_match(
        ethn,
        "1" ~ "Hispanic",
        "2" ~ "Non-hispanic",
        .default = "Unknown"
      ) |> 
        factor(
          levels = c(
            "Hispanic",
            "Non-hispanic",
            "Unknown"
          )
        ),
      race_nih = case_match(
        race_nih,
        "2" ~ "White",
        "3" ~ "Black",
        "4" ~ "Asian",
        "5" ~ "American Indian/Alaska Native",
        "6" ~ "Native Hawaiian/Pacific Islander",
        "8" ~ "More than One Race",
        .default = "Unknown"
      ) |> 
        factor(
          levels = c(
            "White",
            "Black",
            "Asian",
            "American Indian/Alaska Native",
            "Native Hawaiian/Pacific Islander",
            "More than One Race",
            "Unknown"
          )
        ),
      hhold_3lvl = case_match(
        hhold_3lvl,
        "1" ~ "< 50k",
        "2" ~ "50k to 100k",
        "3" ~ "> 100k",
        "777" ~ "Decline to answer",
        "999" ~ "Don't know",
        .default = "Unknown"
      ) |> 
        factor(
          levels = c(
            "< 50k",
            "50k to 100k",
            "> 100k",
            "Decline to answer",
            "Don't know",
            "Unknown"
          )
        ),
      edu = case_match(
        edu,
        "1" ~ "Up to high school",
        "2" ~ "High school diploma/GED",
        "3" ~ "Some college",
        "4" ~ "Bachelor’s degree",
        "5" ~ "Grad school/Professional degree",
        .default = "Unknown"
      ) |> 
        factor(
          levels = c(
            "Up to high school",
            "High school diploma/GED",
            "Some college",
            "Bachelor’s degree",
            "Grad school/Professional degree",
            "Unknown"
          )
        )
    )
}

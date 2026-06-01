# ===========================================================================
# I/O helpers — figure and table saving
# Figures are written at 300 DPI in BOTH png and svg to outputs/figures/.
# Tables are written as CSV to outputs/tables/. Each helper returns the file
# paths so targets can track them as file outputs.
# ===========================================================================

save_figure <- function(plot, name, cfg = pipeline_config(),
                        width = 8, height = 6) {
  dir.create(cfg$out_fig_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(cfg$out_fig_dir, paste0(name, ".png"))
  svg_path <- file.path(cfg$out_fig_dir, paste0(name, ".svg"))

  ggplot2::ggsave(png_path, plot, width = width, height = height,
                  dpi = cfg$fig_dpi, bg = "white")
  ggplot2::ggsave(svg_path, plot, width = width, height = height, bg = "white")

  c(png_path, svg_path)
}

save_table <- function(df, name, cfg = pipeline_config()) {
  dir.create(cfg$out_tab_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cfg$out_tab_dir, paste0(name, ".csv"))
  readr::write_csv(df, path)
  path
}

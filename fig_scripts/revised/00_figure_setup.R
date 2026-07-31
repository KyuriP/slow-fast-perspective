# ==============================================================================
# Shared setup for revised manuscript figures
#
# Run all scripts from the repository root:
#   source("fig_scripts/revised/make_revised_figures.R")
# ==============================================================================

if (!file.exists("R/00_setup.R")) {
  stop(
    "Run this script from the slow-fast-perspective repository root. ",
    "R/00_setup.R was not found."
  )
}

required_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "broom",
  "ggplot2",
  "patchwork",
  "scales",
  "qgraph"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following package(s) before running the figure scripts: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Use the same basic palette and symptom order as the current repository.
LOW_COLOR  <- "#1B9E77"
HIGH_COLOR <- "#D95F02"
BLUE_COLOR <- "#2C7FB8"
GREY_COLOR <- "grey62"
DARK_GREY  <- "grey25"
LIGHT_GREY <- "grey90"

GROUP_COLORS <- c(
  "Low Slow Risk"  = LOW_COLOR,
  "High Slow Risk" = HIGH_COLOR,
  "Dutch"          = "#7570B3",
  "Non-Dutch"      = "#E7298A",
  "Younger"        = "#66A61E",
  "Older"          = "#E6AB02"
)

SYMPTOMS <- c("anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui")

SYMPTOM_LABELS <- c(
  anh = "Anhedonia",
  dep = "Depressed mood",
  slp = "Sleep problems",
  ene = "Low energy",
  app = "Appetite change",
  glt = "Guilt",
  con = "Concentration problems",
  mot = "Psychomotor change",
  sui = "Suicidality"
)

FONT_FAMILY <- "Palatino"

theme_revised <- theme_minimal(
  base_size = 12.5,
  base_family = FONT_FAMILY
) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 11.5, color = DARK_GREY),
    axis.title = element_text(size = 11.5),
    axis.text = element_text(size = 10.5, color = "black"),
    strip.text = element_text(face = "bold", size = 11.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(8, 10, 8, 10)
  )

dir.create("figs/revised", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/revised", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figure_data", recursive = TRUE, showWarnings = FALSE)

ensure_analysis_data <- function() {
  if (!exists("analysis_df", envir = .GlobalEnv, inherits = FALSE)) {
    source("R/00_setup.R", local = .GlobalEnv)
    source("R/01_prepare_analysis_data.R", local = .GlobalEnv)
    source("R/02_model_helpers.R", local = .GlobalEnv)
    source("R/03_prepare_binary_data.R", local = .GlobalEnv)
  }

  required_objects <- c("analysis_df", "analysis_df_bin", "symptoms")
  missing_objects <- required_objects[
    !vapply(
      required_objects,
      exists,
      logical(1),
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ]

  if (length(missing_objects) > 0L) {
    stop(
      "The data-preparation scripts did not create: ",
      paste(missing_objects, collapse = ", ")
    )
  }

  invisible(TRUE)
}

save_revised_figure <- function(
    plot,
    filename,
    width,
    height,
    dpi = 320) {

  pdf_path <- file.path("figs", "revised", paste0(filename, ".pdf"))
  png_path <- file.path("figs", "revised", paste0(filename, ".png"))

  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  message("Saved: ", pdf_path)
  message("Saved: ", png_path)

  invisible(c(pdf = pdf_path, png = png_path))
}

require_result_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required result file was not found:\n  ",
      path,
      "\nRun the corresponding analysis script first."
    )
  }

  invisible(path)
}

binomial_ci <- function(x, n, level = 0.95) {
  if (n <= 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  out <- suppressWarnings(
    stats::prop.test(
      x = x,
      n = n,
      conf.level = level,
      correct = FALSE
    )
  )

  c(
    lower = unname(out$conf.int[1]),
    upper = unname(out$conf.int[2])
  )
}

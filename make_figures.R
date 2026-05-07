# ───────────────────────────────────────────────────────────────────────────────
# make_figures.R
# Master script for generating manuscript figures
# ───────────────────────────────────────────────────────────────────────────────

rm(list = ls())

source("R/00_setup.R")
source("R/01_prepare_analysis_data.R")
source("R/02_model_helpers.R")
source("R/03_fit_shared_model.R")

source("fig_scripts/fig_descriptives.R")
ggsave(
  filename = "figs/figure_01_descriptives.pdf",
  plot = fig_descriptives,
  width = 7,
  height = 8,
  units = "in"
)

source("fig_scripts/fig_threshold_validation.R")
ggsave(
  filename = "figs/figure_02_threshold_validation.pdf",
  plot = fig_threshold_validation,
  width = 8.2,
  height = 5.4,
  units = "in"
)

source("fig_scripts/fig_slow_risk_symptom_map.R")
ggsave(
  filename = "figs/figure_03_slow_risk_symptom_map.pdf",
  plot = fig_slow_risk_symptom_map,
  width = 12,
  height = 12,
  units = "in"
)


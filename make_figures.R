# ───────────────────────────────────────────────────────────────────────────────
# make_figures.R
# Master script for generating manuscript figures
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────
# 1) Setup and data preparation
# ───────────────────────────────────────────────────────────────────────────────

source("R/00_setup.R")
source("R/01_prepare_analysis_data.R")
source("R/02_model_helpers.R")
source("R/03_prepare_binary_data.R")


# ───────────────────────────────────────────────────────────────────────────────
# 2) Figure 1: descriptives
# ───────────────────────────────────────────────────────────────────────────────

source("fig_scripts/fig_descriptives.R")

ggsave(
  filename = "figs/figure_01_descriptives.pdf",
  plot = fig_descriptives,
  width = 7,
  height = 8,
  units = "in"
)


# ───────────────────────────────────────────────────────────────────────────────
# 3) Figure 2: threshold validation / baseline activation
#    This script fits the joint Ising model and creates joint_global.
# ───────────────────────────────────────────────────────────────────────────────

source("fig_scripts/fig_threshold_validation.R")

ggsave(
  filename = "figs/figure_02_threshold_validation.pdf",
  plot = fig_threshold_validation,
  width = 8.2,
  height = 5.4,
  units = "in"
)


# ───────────────────────────────────────────────────────────────────────────────
# 4) Figure 3: Slow Risk Load symptom map
#    This script uses joint_global$J from the joint Ising model above.
# ───────────────────────────────────────────────────────────────────────────────

source("fig_scripts/fig_slow_risk_symptom_map.R")

ggsave(
  filename = "figs/figure_03_slow_risk_symptom_map.pdf",
  plot = fig_slow_risk_symptom_map,
  width = 12,
  height = 12,
  units = "in"
)


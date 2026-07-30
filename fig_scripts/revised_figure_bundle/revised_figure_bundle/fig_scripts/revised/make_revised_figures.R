# ==============================================================================
# Master script for all revised manuscript figures
#
# Run from the repository root:
#   source("fig_scripts/revised/make_revised_figures.R")
# ==============================================================================

source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/00_revised_figure_setup.R")

message("\nCreating revised main figures...\n")

source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/figure_01_slow_fast_concept.R")
source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/figure_02_slowrisk_descriptives.R")
source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/figure_03_activation_differences.R")
source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/figure_04_slowrisk_coefficients.R")

message("\nCreating model-comparison table...\n")

source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/table_02_model_comparison.R")

message("\nCreating supplementary figures...\n")

source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/supplement_S1_secondary_descriptives.R")
source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/supplement_S2_activation_sensitivity.R")
source("fig_scripts/revised_figure_bundle/revised_figure_bundle/fig_scripts/revised/supplement_S3_network_comparison.R")

message(
  "\nAll revised figures are complete.\n",
  "Main and supplementary figures: figs/revised/\n",
  "Model-comparison table: tables/revised/\n",
  "Figure source data: results/figure_data/\n"
)

# ==============================================================================
# Final NCT: Low versus High Slow Risk, 10,000 permutations
# ==============================================================================

SEED <- 20260729L
NCT_ITERATIONS <- 10000L
NCT_GAMMA <- 0.25
NCT_EDGE_ADJUSTMENT <- "holm"

source("R/00_setup.R")
source("R/01_prepare_analysis_data.R")
source("R/02_model_helpers.R")
source("R/03_prepare_binary_data.R")

if (!requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  stop("Install NetworkComparisonTest before running this script.")
}

spec <- list(
  label = "Low versus High Slow Risk",
  group_var = "slow_group",
  group_levels = c("Low Slow Risk", "High Slow Risk")
)

required_columns <- c(symptoms, spec$group_var)

dat <- analysis_df_bin[, required_columns, drop = FALSE]
dat <- dat[complete.cases(dat), , drop = FALSE]

dat <- dat[
  dat[[spec$group_var]] %in% spec$group_levels,
  ,
  drop = FALSE
]

dat[[spec$group_var]] <- factor(
  as.character(dat[[spec$group_var]]),
  levels = spec$group_levels
)

for (v in symptoms) {
  dat[[v]] <- as.integer(dat[[v]])
}

group_1 <- dat[
  dat[[spec$group_var]] == spec$group_levels[1],
  symptoms,
  drop = FALSE
]

group_2 <- dat[
  dat[[spec$group_var]] == spec$group_levels[2],
  symptoms,
  drop = FALSE
]

output_dir <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_10000"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(SEED)

nct_result <- NetworkComparisonTest::NCT(
  data1 = group_1,
  data2 = group_2,
  gamma = NCT_GAMMA,
  it = NCT_ITERATIONS,
  binary.data = TRUE,
  paired = FALSE,
  weighted = TRUE,
  AND = TRUE,
  abs = TRUE,
  test.edges = TRUE,
  edges = "all",
  progressbar = TRUE,
  make.positive.definite = TRUE,
  p.adjust.methods = NCT_EDGE_ADJUSTMENT,
  test.centrality = FALSE,
  verbose = TRUE
)

saveRDS(
  nct_result,
  file.path(output_dir, "nct_result_10000.rds")
)

strengths <- as.numeric(nct_result$glstrinv.sep)

nct_summary <- data.frame(
  comparison = spec$label,
  group_1 = spec$group_levels[1],
  group_2 = spec$group_levels[2],
  n_group_1 = nrow(group_1),
  n_group_2 = nrow(group_2),
  network_structure_statistic_M =
    as.numeric(nct_result$nwinv.real),
  network_structure_p =
    as.numeric(nct_result$nwinv.pval),
  global_strength_group_1 = strengths[1],
  global_strength_group_2 = strengths[2],
  global_strength_difference_S =
    as.numeric(nct_result$glstrinv.real),
  global_strength_p =
    as.numeric(nct_result$glstrinv.pval),
  permutations = NCT_ITERATIONS
)

write.csv(
  nct_summary,
  file.path(output_dir, "nct_summary_10000.csv"),
  row.names = FALSE
)

if (!is.null(nct_result$einv.pvals)) {
  write.csv(
    as.data.frame(nct_result$einv.pvals),
    file.path(output_dir, "nct_edge_tests_10000.csv"),
    row.names = FALSE
  )
}

write.csv(
  nct_result$nw1,
  file.path(output_dir, "nct_network_low_risk_10000.csv")
)

write.csv(
  nct_result$nw2,
  file.path(output_dir, "nct_network_high_risk_10000.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message(
  "\nFinal NCT complete. Results saved to:\n",
  normalizePath(output_dir)
)
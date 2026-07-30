# ==============================================================================
# R/05_nct_final_10000.R
#
# Final Network Comparison Test for Low versus High Slow Risk:
#   - binary Ising networks
#   - 10,000 permutations
#   - gamma = 0.25
#   - AND rule
#   - Holm adjustment across all 36 edges
#
# Run from the repository root:
#   source("R/05_nct_final_10000.R")
# ==============================================================================

SEED <- 20260729L
NCT_ITERATIONS <- 10000L
NCT_GAMMA <- 0.25
NCT_EDGE_ADJUSTMENT <- "holm"

source_files <- c(
  "R/00_setup.R",
  "R/01_prepare_analysis_data.R",
  "R/02_model_helpers.R",
  "R/03_prepare_binary_data.R"
)

missing_source_files <- source_files[!file.exists(source_files)]

if (length(missing_source_files) > 0L) {
  stop(
    "Run this script from the repository root. Missing file(s): ",
    paste(missing_source_files, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(source_files, source))

if (!requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  stop(
    "Package 'NetworkComparisonTest' is required. Install it before running.",
    call. = FALSE
  )
}

spec <- list(
  label = "Low versus High Slow Risk",
  group_var = "slow_group",
  group_levels = c("Low Slow Risk", "High Slow Risk")
)

required_columns <- c(symptoms, spec$group_var)
missing_columns <- setdiff(required_columns, names(analysis_df_bin))

if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

dat <- analysis_df_bin[, required_columns, drop = FALSE]
dat <- dat[stats::complete.cases(dat), , drop = FALSE]
dat <- dat[
  dat[[spec$group_var]] %in% spec$group_levels,
  ,
  drop = FALSE
]

dat[[spec$group_var]] <- factor(
  as.character(dat[[spec$group_var]]),
  levels = spec$group_levels
)

for (symptom in symptoms) {
  dat[[symptom]] <- as.integer(dat[[symptom]])

  invalid_values <- setdiff(
    unique(dat[[symptom]]),
    c(0L, 1L)
  )

  if (length(invalid_values) > 0L) {
    stop(
      "Symptom '", symptom,
      "' contains values other than 0 and 1.",
      call. = FALSE
    )
  }
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

if (nrow(group_1) == 0L || nrow(group_2) == 0L) {
  stop("Both Slow Risk groups must contain observations.", call. = FALSE)
}

constant_group_1 <- names(
  which(vapply(group_1, function(x) length(unique(x)) < 2L, logical(1)))
)
constant_group_2 <- names(
  which(vapply(group_2, function(x) length(unique(x)) < 2L, logical(1)))
)

if (length(constant_group_1) > 0L || length(constant_group_2) > 0L) {
  stop(
    "Constant symptoms detected. Low group: ",
    paste(constant_group_1, collapse = ", "),
    "; High group: ",
    paste(constant_group_2, collapse = ", "),
    call. = FALSE
  )
}

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

if (length(strengths) < 2L) {
  stop(
    "NCT output did not contain two group-specific global strengths.",
    call. = FALSE
  )
}

nct_summary <- data.frame(
  comparison = spec$label,
  group_1 = spec$group_levels[1],
  group_2 = spec$group_levels[2],
  n_group_1 = nrow(group_1),
  n_group_2 = nrow(group_2),
  network_structure_statistic_M =
    as.numeric(nct_result$nwinv.real)[1],
  network_structure_p =
    as.numeric(nct_result$nwinv.pval)[1],
  global_strength_group_1 = strengths[1],
  global_strength_group_2 = strengths[2],
  global_strength_difference_S =
    as.numeric(nct_result$glstrinv.real)[1],
  global_strength_p =
    as.numeric(nct_result$glstrinv.pval)[1],
  permutations = NCT_ITERATIONS,
  gamma = NCT_GAMMA,
  edge_adjustment = NCT_EDGE_ADJUSTMENT,
  stringsAsFactors = FALSE
)

utils::write.csv(
  nct_summary,
  file.path(output_dir, "nct_summary_10000.csv"),
  row.names = FALSE
)

if (!is.null(nct_result$einv.pvals)) {
  utils::write.csv(
    as.data.frame(nct_result$einv.pvals),
    file.path(output_dir, "nct_edge_tests_package_output_10000.csv"),
    row.names = FALSE
  )
}

utils::write.csv(
  nct_result$nw1,
  file.path(output_dir, "nct_network_low_risk_10000.csv"),
  row.names = TRUE
)

utils::write.csv(
  nct_result$nw2,
  file.path(output_dir, "nct_network_high_risk_10000.csv"),
  row.names = TRUE
)

utils::write.csv(
  data.frame(
    group = spec$group_levels,
    n = c(nrow(group_1), nrow(group_2)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "group_counts_10000.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message(
  "\nFinal NCT complete.\nResults saved to:\n",
  normalizePath(output_dir, mustWork = FALSE)
)

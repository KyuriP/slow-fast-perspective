# ==============================================================================
# R/07_build_appendix_outputs.R
#
# Builds all appendix outputs used in the manuscript:
#
#   Table A1  Optimization and Hessian diagnostics
#   Table A2  Activation estimates, fully group-specific model
#   Table A3  Activation estimates, common-interaction model
#   Table A4  Edge-specific Network Comparison Test results
#   Figure A1 Sensitivity of activation differences to interaction specification
#
# Run from the repository root:
#
#   source("R/07_build_appendix_outputs.R")
#
# Outputs are written to:
#
#   results/appendix/
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

EXACT_RESULTS_PATH <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "exact_multigroup_ising",
  "exact_multigroup_models.rds"
)

# Preferred final NCT object produced by the 10,000-permutation analysis.
NCT_RESULTS_PATH <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_10000",
  "nct_result_10000.rds"
)

# Backward-compatible fallback for repositories that still store the same
# final NCT object at the older location.
NCT_RESULTS_FALLBACK <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_result.rds"
)

OUTPUT_DIR <- file.path("results", "appendix")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 1. Package and input checks
# ------------------------------------------------------------------------------

required_packages <- c("ggplot2", "knitr")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (!file.exists(EXACT_RESULTS_PATH)) {
  stop(
    "Could not find exact-model results at:\n",
    EXACT_RESULTS_PATH,
    call. = FALSE
  )
}

if (!file.exists(NCT_RESULTS_PATH)) {
  if (file.exists(NCT_RESULTS_FALLBACK)) {
    warning(
      "Preferred final NCT object was not found at:\n",
      NCT_RESULTS_PATH,
      "\nUsing fallback object:\n",
      NCT_RESULTS_FALLBACK,
      "\nConfirm that this is the final 10,000-permutation result."
    )
    NCT_RESULTS_PATH <- NCT_RESULTS_FALLBACK
  } else {
    stop(
      "Could not find the final NCT results at either:\n",
      NCT_RESULTS_PATH,
      "\nor:\n",
      NCT_RESULTS_FALLBACK,
      call. = FALSE
    )
  }
}


# ------------------------------------------------------------------------------
# 2. Load saved model objects
# ------------------------------------------------------------------------------

exact_results <- readRDS(EXACT_RESULTS_PATH)
nct_result <- readRDS(NCT_RESULTS_PATH)

required_exact_elements <- c("fits", "symptoms", "pair_index")

missing_exact_elements <- setdiff(
  required_exact_elements,
  names(exact_results)
)

if (length(missing_exact_elements) > 0L) {
  stop(
    "The exact-results object is missing: ",
    paste(missing_exact_elements, collapse = ", "),
    call. = FALSE
  )
}

fits <- exact_results$fits
symptoms <- exact_results$symptoms
pair_index <- exact_results$pair_index

p <- length(symptoms)
q <- nrow(pair_index)

expected_fit_names <- c(
  "free_J_free_h",
  "equal_J_free_h",
  "free_J_equal_h",
  "equal_J_equal_h"
)

missing_fits <- setdiff(expected_fit_names, names(fits))

if (length(missing_fits) > 0L) {
  stop(
    "The exact-results object is missing these fitted models: ",
    paste(missing_fits, collapse = ", "),
    call. = FALSE
  )
}

symptom_labels <- c(
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

if (!all(symptoms %in% names(symptom_labels))) {
  stop(
    "At least one symptom code is not defined in symptom_labels: ",
    paste(setdiff(symptoms, names(symptom_labels)), collapse = ", "),
    call. = FALSE
  )
}

model_labels <- c(
  free_J_free_h = "Fully group-specific",
  equal_J_free_h = "Common interactions",
  free_J_equal_h = "Common activation",
  equal_J_equal_h = "Fully invariant"
)


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

write_csv <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(TABLE_DIR, filename),
    row.names = FALSE,
    na = ""
  )
}

format_fixed <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    formatC(x, digits = digits, format = "f")
  )
}

format_scientific <- function(x, digits = 2) {
  ifelse(
    is.na(x),
    "",
    formatC(x, digits = digits, format = "e")
  )
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < .001,
      "$<.001$",
      sub("^0", "", formatC(x, digits = 3, format = "f"))
    )
  )
}

write_latex_table <- function(
    x,
    filename,
    caption,
    label,
    longtable = FALSE) {

  latex <- knitr::kable(
    x,
    format = "latex",
    booktabs = TRUE,
    longtable = longtable,
    caption = caption,
    label = label,
    row.names = FALSE,
    escape = FALSE,
    linesep = ""
  )

  writeLines(
    latex,
    con = file.path(TABLE_DIR, filename)
  )
}

check_vcov <- function(fit, model_name) {
  if (is.null(fit$vcov)) {
    stop(
      "No covariance matrix was saved for model '",
      model_name,
      "'. Rerun the exact-model script after enabling hessian = TRUE ",
      "and saving solve(hessian) as fit$vcov.",
      call. = FALSE
    )
  }

  if (!is.matrix(fit$vcov)) {
    stop(
      "fit$vcov is not a matrix for model '",
      model_name,
      "'.",
      call. = FALSE
    )
  }
}

get_activation_indices <- function(model_name) {

  if (model_name == "free_J_free_h") {
    # Parameter order:
    # low h, low J, high h, high J
    index_low <- seq_len(p)
    index_high <- p + q + seq_len(p)
  } else if (model_name == "equal_J_free_h") {
    # Parameter order:
    # low h, high h, common J
    index_low <- seq_len(p)
    index_high <- p + seq_len(p)
  } else {
    stop(
      "Activation contrasts can only be extracted from models with ",
      "group-specific activation parameters.",
      call. = FALSE
    )
  }

  list(
    low = index_low,
    high = index_high
  )
}

make_activation_table <- function(model_name) {

  fit <- fits[[model_name]]
  check_vcov(fit, model_name)

  index <- get_activation_indices(model_name)

  h_low <- fit$par[index$low]
  h_high <- fit$par[index$high]

  se_h_low <- sqrt(diag(fit$vcov)[index$low])
  se_h_high <- sqrt(diag(fit$vcov)[index$high])

  contrast_matrix <- matrix(
    0,
    nrow = p,
    ncol = length(fit$par)
  )

  for (i in seq_len(p)) {
    contrast_matrix[i, index$low[i]] <- -1
    contrast_matrix[i, index$high[i]] <- 1
  }

  variance_delta_h <- diag(
    contrast_matrix %*%
      fit$vcov %*%
      t(contrast_matrix)
  )

  se_delta_h <- sqrt(pmax(variance_delta_h, 0))
  delta_h <- h_high - h_low

  data.frame(
    Symptom_code = symptoms,
    Symptom = unname(symptom_labels[symptoms]),
    h_Low = h_low,
    SE_h_Low = se_h_low,
    h_High = h_high,
    SE_h_High = se_h_high,
    Delta_h = delta_h,
    SE_Delta_h = se_delta_h,
    CI_lower = delta_h - 1.96 * se_delta_h,
    CI_upper = delta_h + 1.96 * se_delta_h,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 4. Table A1: optimization and Hessian diagnostics
# ------------------------------------------------------------------------------

appendix_table_A1 <- do.call(
  rbind,
  lapply(
    expected_fit_names,
    function(model_key) {

      fit <- fits[[model_key]]

      data.frame(
        Model_key = model_key,
        Model = unname(model_labels[model_key]),
        Parameters = fit$npar,
        Convergence_code = fit$convergence,
        Maximum_absolute_gradient = fit$max_abs_gradient,
        Mean_absolute_gradient = fit$mean_abs_gradient,
        Minimum_Hessian_eigenvalue =
          fit$min_hessian_eigenvalue,
        Hessian_condition_number =
          fit$hessian_condition_number,
        Positive_definite_Hessian =
          fit$hessian_positive_definite,
        Covariance_matrix_available =
          !is.null(fit$vcov),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_csv(
  appendix_table_A1,
  "table_A1_optimization_diagnostics.csv"
)

table_A1_latex <- appendix_table_A1[
  ,
  c(
    "Model",
    "Parameters",
    "Convergence_code",
    "Maximum_absolute_gradient",
    "Minimum_Hessian_eigenvalue",
    "Hessian_condition_number",
    "Covariance_matrix_available"
  )
]

table_A1_latex$Maximum_absolute_gradient <-
  format_scientific(
    table_A1_latex$Maximum_absolute_gradient,
    digits = 2
  )

table_A1_latex$Minimum_Hessian_eigenvalue <-
  format_fixed(
    table_A1_latex$Minimum_Hessian_eigenvalue,
    digits = 2
  )

table_A1_latex$Hessian_condition_number <-
  format_fixed(
    table_A1_latex$Hessian_condition_number,
    digits = 2
  )

table_A1_latex$Covariance_matrix_available <-
  ifelse(
    table_A1_latex$Covariance_matrix_available,
    "Yes",
    "No"
  )

names(table_A1_latex) <- c(
  "Model",
  "$k$",
  "Code",
  "Max. $|g|$",
  "Min. Hessian eigenvalue",
  "Condition number",
  "Covariance"
)

write_latex_table(
  table_A1_latex,
  "table_A1_optimization_diagnostics.tex",
  paste(
    "Optimization diagnostics for the four exact multigroup Ising models.",
    "A convergence code of zero indicates successful termination."
  ),
  "appendix_optimization"
)


# ------------------------------------------------------------------------------
# 5. Table A2: fully group-specific activation estimates
# ------------------------------------------------------------------------------

appendix_table_A2 <- make_activation_table(
  "free_J_free_h"
)

write_csv(
  appendix_table_A2,
  "table_A2_activation_fully_group_specific.csv"
)

table_A2_latex <- appendix_table_A2[
  ,
  c(
    "Symptom",
    "h_Low",
    "h_High",
    "Delta_h",
    "SE_Delta_h",
    "CI_lower",
    "CI_upper"
  )
]

numeric_columns_A2 <- setdiff(
  names(table_A2_latex),
  "Symptom"
)

table_A2_latex[numeric_columns_A2] <- lapply(
  table_A2_latex[numeric_columns_A2],
  format_fixed,
  digits = 3
)

names(table_A2_latex) <- c(
  "Symptom",
  "$h_{\\mathrm{Low}}$",
  "$h_{\\mathrm{High}}$",
  "$\\Delta h$",
  "SE",
  "95\\% CI lower",
  "95\\% CI upper"
)

write_latex_table(
  table_A2_latex,
  "table_A2_activation_fully_group_specific.tex",
  paste(
    "Symptom-activation estimates from the fully group-specific exact",
    "Ising model. Both interaction and activation parameters were",
    "estimated separately in the Low and High Slow Risk groups."
  ),
  "appendix_activation_free"
)


# ------------------------------------------------------------------------------
# 6. Table A3: common-interaction activation estimates
# ------------------------------------------------------------------------------

appendix_table_A3 <- make_activation_table(
  "equal_J_free_h"
)

write_csv(
  appendix_table_A3,
  "table_A3_activation_common_interactions.csv"
)

table_A3_latex <- appendix_table_A3[
  ,
  c(
    "Symptom",
    "h_Low",
    "h_High",
    "Delta_h",
    "SE_Delta_h",
    "CI_lower",
    "CI_upper"
  )
]

numeric_columns_A3 <- setdiff(
  names(table_A3_latex),
  "Symptom"
)

table_A3_latex[numeric_columns_A3] <- lapply(
  table_A3_latex[numeric_columns_A3],
  format_fixed,
  digits = 3
)

names(table_A3_latex) <- c(
  "Symptom",
  "$h_{\\mathrm{Low}}$",
  "$h_{\\mathrm{High}}$",
  "$\\Delta h$",
  "SE",
  "95\\% CI lower",
  "95\\% CI upper"
)

write_latex_table(
  table_A3_latex,
  "table_A3_activation_common_interactions.tex",
  paste(
    "Symptom-activation estimates from the exact Ising model with",
    "interactions constrained to equality and activation parameters",
    "estimated separately across groups."
  ),
  "appendix_activation_commonJ"
)


# ------------------------------------------------------------------------------
# 7. Figure A1: sensitivity of activation estimates
# ------------------------------------------------------------------------------

figure_data <- rbind(
  transform(
    appendix_table_A2,
    Specification = "Group-specific interactions"
  ),
  transform(
    appendix_table_A3,
    Specification = "Common interactions"
  )
)

symptom_order <- appendix_table_A2$Symptom[
  order(appendix_table_A2$Delta_h)
]

figure_data$Symptom <- factor(
  figure_data$Symptom,
  levels = symptom_order
)

figure_data$Specification <- factor(
  figure_data$Specification,
  levels = c(
    "Group-specific interactions",
    "Common interactions"
  )
)

position <- ggplot2::position_dodge(width = 0.55)

figure_A1 <- ggplot2::ggplot(
  figure_data,
  ggplot2::aes(
    x = Symptom,
    y = Delta_h,
    shape = Specification,
    group = Specification
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    position = position,
    width = 0.12,
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    position = position,
    size = 2.7
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Activation difference (High minus Low Slow Risk)",
    shape = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major.y = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "figure_A1_activation_sensitivity.pdf"
  ),
  plot = figure_A1,
  width = 7.5,
  height = 5.5,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "figure_A1_activation_sensitivity.png"
  ),
  plot = figure_A1,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 400
)


# ------------------------------------------------------------------------------
# 8. Table A4: all 36 NCT edge comparisons
# ------------------------------------------------------------------------------

required_nct_elements <- c("nw1", "nw2", "einv.perm")

missing_nct_elements <- setdiff(
  required_nct_elements,
  names(nct_result)
)

if (length(missing_nct_elements) > 0L) {
  stop(
    "The NCT object is missing: ",
    paste(missing_nct_elements, collapse = ", "),
    call. = FALSE
  )
}

edge_positions <- cbind(
  pair_index[, 1],
  pair_index[, 2]
)

edge_low <- as.numeric(
  nct_result$nw1[edge_positions]
)

edge_high <- as.numeric(
  nct_result$nw2[edge_positions]
)

perm <- nct_result$einv.perm

if (length(dim(perm)) != 3L) {
  stop(
    "nct_result$einv.perm must be a three-dimensional array.",
    call. = FALSE
  )
}

n_permutations <- dim(perm)[3]

observed_abs_difference <- abs(
  edge_high - edge_low
)

p_raw <- vapply(
  seq_len(q),
  function(k) {

    i <- pair_index[k, 1]
    j <- pair_index[k, 2]

    permuted_abs_difference <- abs(
      perm[i, j, seq_len(n_permutations)]
    )

    (
      sum(
        permuted_abs_difference >=
          observed_abs_difference[k]
      ) + 1
    ) / (n_permutations + 1)
  },
  numeric(1)
)

p_holm <- stats::p.adjust(
  p_raw,
  method = "holm"
)

if (length(p_raw) != q || length(p_holm) != q) {
  stop(
    "The number of edgewise p-values does not equal the number of ",
    "unique symptom pairs.",
    call. = FALSE
  )
}

appendix_table_A4 <- data.frame(
  Symptom_1 = unname(
    symptom_labels[symptoms[pair_index[, 1]]]
  ),
  Symptom_2 = unname(
    symptom_labels[symptoms[pair_index[, 2]]]
  ),
  Edge_Low_Slow_Risk = edge_low,
  Edge_High_Slow_Risk = edge_high,
  Difference_High_minus_Low =
    edge_high - edge_low,
  Absolute_difference =
    abs(edge_high - edge_low),
  P_unadjusted = p_raw,
  P_Holm = p_holm,
  stringsAsFactors = FALSE
)

appendix_table_A4 <- appendix_table_A4[
  order(
    appendix_table_A4$P_Holm,
    appendix_table_A4$P_unadjusted
  ),
  ,
  drop = FALSE
]

write_csv(
  appendix_table_A4,
  "table_A4_nct_edge_comparisons.csv"
)

table_A4_latex <- appendix_table_A4

numeric_edge_columns <- c(
  "Edge_Low_Slow_Risk",
  "Edge_High_Slow_Risk",
  "Difference_High_minus_Low",
  "Absolute_difference"
)

table_A4_latex[numeric_edge_columns] <- lapply(
  table_A4_latex[numeric_edge_columns],
  format_fixed,
  digits = 3
)

table_A4_latex$P_unadjusted <-
  format_p(table_A4_latex$P_unadjusted)

table_A4_latex$P_Holm <-
  format_p(table_A4_latex$P_Holm)

names(table_A4_latex) <- c(
  "Symptom 1",
  "Symptom 2",
  "Low",
  "High",
  "$\\Delta J$",
  "$|\\Delta J|$",
  "$p$",
  "$p_{\\mathrm{Holm}}$"
)

write_latex_table(
  table_A4_latex,
  "table_A4_nct_edge_comparisons.tex",
  paste(
    "Edge-specific Network Comparison Test results for the Low and",
    "High Slow Risk groups. The signed difference is High minus Low.",
    "Both unadjusted and Holm-adjusted permutation p-values are shown."
  ),
  "appendix_nct_edges",
  longtable = TRUE
)


# ------------------------------------------------------------------------------
# 9. Verification summary
# ------------------------------------------------------------------------------

cat("\nAppendix outputs created successfully.\n\n")

cat("Exact-model input:\n", normalizePath(EXACT_RESULTS_PATH), "\n\n")
cat("NCT input:\n", normalizePath(NCT_RESULTS_PATH), "\n\n")

cat("Table A1: optimization diagnostics\n")
print(appendix_table_A1)

cat("\nTable A2: fully group-specific activation contrasts\n")
print(
  data.frame(
    Symptom = appendix_table_A2$Symptom,
    round(
      appendix_table_A2[
        ,
        c(
          "Delta_h",
          "SE_Delta_h",
          "CI_lower",
          "CI_upper"
        )
      ],
      3
    ),
    row.names = NULL
  )
)

cat("\nTable A3: common-interaction activation contrasts\n")
print(
  data.frame(
    Symptom = appendix_table_A3$Symptom,
    round(
      appendix_table_A3[
        ,
        c(
          "Delta_h",
          "SE_Delta_h",
          "CI_lower",
          "CI_upper"
        )
      ],
      3
    ),
    row.names = NULL
  )
)

cat("\nFive smallest Holm-adjusted NCT edge p-values\n")
print(
  head(
    appendix_table_A4[
      ,
      c(
        "Symptom_1",
        "Symptom_2",
        "Difference_High_minus_Low",
        "P_unadjusted",
        "P_Holm"
      )
    ],
    5
  )
)

cat(
  "\nFiles written to:\n",
  normalizePath(OUTPUT_DIR),
  "\n"
)

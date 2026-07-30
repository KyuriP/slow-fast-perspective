# ==============================================================================
# R/06_activation_difference_ci.R
#
# Wald confidence intervals for activation differences in the fully free exact
# two-group Ising model:
#
#   free_J_free_h:
#     J_low != J_high
#     h_low != h_high
#
# The exact Fisher information is computed from the 2^9 = 512-state
# distribution. Because the two groups are independent, for each symptom:
#
#   Var(h_high - h_low) = Var(h_high) + Var(h_low).
#
# This script does not refit any model.
#
# Run from the repository root:
#
#   source("R/06_activation_difference_ci.R")
#
# ==============================================================================


# ==============================================================================
# 0. PATHS
# ==============================================================================

RESULT_DIR <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "exact_multigroup_ising"
)

MODEL_FILE <- file.path(
  RESULT_DIR,
  "exact_multigroup_models.rds"
)

GROUP_COUNT_FILE <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "group_counts.csv"
)

OUTPUT_FILE <- file.path(
  RESULT_DIR,
  "activation_differences_free_J_free_h_with_95CI.csv"
)

FIGURE_PDF <- file.path(
  RESULT_DIR,
  "activation_differences_free_J_free_h_95CI.pdf"
)

FIGURE_PNG <- file.path(
  RESULT_DIR,
  "activation_differences_free_J_free_h_95CI.png"
)


# ==============================================================================
# 1. CHECKS AND INPUTS
# ==============================================================================

if (!file.exists(MODEL_FILE)) {
  stop("Could not find model file: ", MODEL_FILE)
}

if (!file.exists(GROUP_COUNT_FILE)) {
  stop("Could not find group-count file: ", GROUP_COUNT_FILE)
}

model_object <- readRDS(MODEL_FILE)
group_counts <- utils::read.csv(
  GROUP_COUNT_FILE,
  stringsAsFactors = FALSE
)

if (!"free_J_free_h" %in% names(model_object$fits)) {
  stop("The RDS file does not contain the free_J_free_h model.")
}

fit <- model_object$fits$free_J_free_h
symptoms <- model_object$symptoms
pair_index <- model_object$pair_index

if (nrow(group_counts) != 2L) {
  stop("group_counts.csv must contain exactly two groups.")
}

n_low <- as.numeric(group_counts$n[1])
n_high <- as.numeric(group_counts$n[2])

p <- length(symptoms)
q <- nrow(pair_index)

expected_parameter_count <- 2L * (p + q)

if (length(fit$par) != expected_parameter_count) {
  stop(
    "Unexpected number of parameters in free_J_free_h: ",
    length(fit$par),
    "; expected ", expected_parameter_count, "."
  )
}


# ==============================================================================
# 2. RECONSTRUCT THE TWO GROUP-SPECIFIC PARAMETER VECTORS
# ==============================================================================

h_low <- fit$par[seq_len(p)]
J_low <- fit$par[p + seq_len(q)]

offset <- p + q
h_high <- fit$par[offset + seq_len(p)]
J_high <- fit$par[offset + p + seq_len(q)]

theta_low <- c(h_low, J_low)
theta_high <- c(h_high, J_high)


# ==============================================================================
# 3. EXACT STATE SPACE AND SUFFICIENT STATISTICS
# ==============================================================================

states <- as.matrix(
  expand.grid(
    replicate(p, c(0, 1), simplify = FALSE)
  )
)

colnames(states) <- symptoms

state_pair_products <-
  states[, pair_index[, 1], drop = FALSE] *
  states[, pair_index[, 2], drop = FALSE]

state_statistics <- cbind(
  states,
  state_pair_products
)


# ==============================================================================
# 4. EXACT FISHER INFORMATION
# ==============================================================================

exact_information <- function(theta, n, state_statistics) {
  eta <- as.numeric(state_statistics %*% theta)

  max_eta <- max(eta)
  weights <- exp(eta - max_eta)
  probabilities <- weights / sum(weights)

  expected_statistics <- as.numeric(
    crossprod(probabilities, state_statistics)
  )

  weighted_statistics <-
    state_statistics * sqrt(probabilities)

  expected_outer_product <- crossprod(weighted_statistics)

  covariance_statistics <-
    expected_outer_product -
    tcrossprod(expected_statistics)

  information <- n * covariance_statistics

  # Symmetrize to remove tiny floating-point asymmetry.
  information <- (information + t(information)) / 2

  information
}

safe_inverse <- function(matrix_x, label) {
  condition_number <- kappa(matrix_x)

  if (!is.finite(condition_number)) {
    stop("Non-finite information-matrix condition number for ", label, ".")
  }

  inverse <- tryCatch(
    solve(matrix_x),
    error = function(e) NULL
  )

  if (is.null(inverse)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(
        "Information matrix was singular for ", label,
        " and MASS is unavailable for a generalized inverse."
      )
    }

    warning(
      "Using a generalized inverse for ", label,
      "; inspect the condition number carefully."
    )

    inverse <- MASS::ginv(matrix_x)
  }

  list(
    inverse = inverse,
    condition_number = condition_number
  )
}

information_low <- exact_information(
  theta = theta_low,
  n = n_low,
  state_statistics = state_statistics
)

information_high <- exact_information(
  theta = theta_high,
  n = n_high,
  state_statistics = state_statistics
)

inverse_low <- safe_inverse(
  information_low,
  "Low Slow Risk"
)

inverse_high <- safe_inverse(
  information_high,
  "High Slow Risk"
)

covariance_low <- inverse_low$inverse
covariance_high <- inverse_high$inverse


# ==============================================================================
# 5. ACTIVATION DIFFERENCE STANDARD ERRORS AND CONFIDENCE INTERVALS
# ==============================================================================

variance_h_low <- diag(covariance_low)[seq_len(p)]
variance_h_high <- diag(covariance_high)[seq_len(p)]

difference <- h_high - h_low

standard_error_difference <- sqrt(
  variance_h_low + variance_h_high
)

z_value <- difference / standard_error_difference

p_value <- 2 * stats::pnorm(
  abs(z_value),
  lower.tail = FALSE
)

critical_value <- stats::qnorm(0.975)

ci_lower <- difference -
  critical_value * standard_error_difference

ci_upper <- difference +
  critical_value * standard_error_difference

symptom_labels <- c(
  anh = "Anhedonia",
  dep = "Depressed mood",
  slp = "Sleep problems",
  ene = "Low energy",
  app = "Appetite change",
  glt = "Guilt",
  con = "Concentration",
  mot = "Psychomotor change",
  sui = "Suicidality"
)

result <- data.frame(
  symptom = symptoms,
  symptom_label = unname(symptom_labels[symptoms]),
  h_low = h_low,
  se_h_low = sqrt(variance_h_low),
  h_high = h_high,
  se_h_high = sqrt(variance_h_high),
  difference_high_minus_low = difference,
  se_difference = standard_error_difference,
  ci_95_lower = ci_lower,
  ci_95_upper = ci_upper,
  z_value = z_value,
  p_value = p_value,
  condition_number_low = inverse_low$condition_number,
  condition_number_high = inverse_high$condition_number,
  stringsAsFactors = FALSE
)

result <- result[
  order(result$difference_high_minus_low),
  ,
  drop = FALSE
]

utils::write.csv(
  result,
  OUTPUT_FILE,
  row.names = FALSE,
  na = ""
)


# ==============================================================================
# 6. KEY ACTIVATION-DIFFERENCE FIGURE
# ==============================================================================

plot_result <- result
plot_result$symptom_label <- factor(
  plot_result$symptom_label,
  levels = plot_result$symptom_label
)

make_plot <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mar = c(4.5, 9.2, 1.2, 1.0),
    las = 1
  )

  y_positions <- seq_len(nrow(plot_result))

  x_range <- range(
    c(
      plot_result$ci_95_lower,
      plot_result$ci_95_upper,
      0
    ),
    finite = TRUE
  )

  padding <- 0.08 * diff(x_range)
  if (!is.finite(padding) || padding == 0) {
    padding <- 0.1
  }

  graphics::plot(
    x = plot_result$difference_high_minus_low,
    y = y_positions,
    xlim = c(
      x_range[1] - padding,
      x_range[2] + padding
    ),
    ylim = c(0.5, nrow(plot_result) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = expression(
      Delta * h[i] == h[i * ", High"] - h[i * ", Low"]
    ),
    pch = 19
  )

  graphics::abline(
    v = 0,
    lty = 2
  )

  graphics::segments(
    x0 = plot_result$ci_95_lower,
    y0 = y_positions,
    x1 = plot_result$ci_95_upper,
    y1 = y_positions
  )

  graphics::points(
    x = plot_result$difference_high_minus_low,
    y = y_positions,
    pch = 19
  )

  graphics::axis(
    side = 2,
    at = y_positions,
    labels = as.character(plot_result$symptom_label),
    las = 1
  )
}

grDevices::pdf(
  FIGURE_PDF,
  width = 7.2,
  height = 5.5,
  useDingbats = FALSE
)

make_plot()
grDevices::dev.off()

grDevices::png(
  FIGURE_PNG,
  width = 1800,
  height = 1350,
  res = 250
)

make_plot()
grDevices::dev.off()


# ==============================================================================
# 7. DIAGNOSTIC OUTPUT
# ==============================================================================

diagnostics <- data.frame(
  quantity = c(
    "n_low",
    "n_high",
    "information_condition_number_low",
    "information_condition_number_high",
    "maximum_absolute_optimization_gradient"
  ),
  value = c(
    n_low,
    n_high,
    inverse_low$condition_number,
    inverse_high$condition_number,
    fit$max_abs_gradient
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  diagnostics,
  file.path(
    RESULT_DIR,
    "activation_ci_diagnostics.csv"
  ),
  row.names = FALSE
)

message(
  "\nActivation-difference confidence intervals complete.\n",
  "CSV: ", normalizePath(OUTPUT_FILE, mustWork = FALSE), "\n",
  "PDF: ", normalizePath(FIGURE_PDF, mustWork = FALSE), "\n",
  "PNG: ", normalizePath(FIGURE_PNG, mustWork = FALSE)
)

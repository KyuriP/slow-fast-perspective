# ==============================================================================
# R/06_validate_activation_covariance.R
#
# Independent numerical validation of activation-difference uncertainty for the
# fully group-specific exact two-group Ising model.
#
# This script:
#   1. reconstructs the exact expected Fisher information separately by group;
#   2. calculates SEs for Delta h = h_high - h_low;
#   3. compares them with SEs derived from the inverse Hessian saved by script 04.
#
# It does not refit models or generate manuscript figures.
#
# Run from the repository root:
#   source("R/06_validate_activation_covariance.R")
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
  "activation_covariance_validation.csv"
)

DIAGNOSTIC_FILE <- file.path(
  RESULT_DIR,
  "activation_covariance_validation_diagnostics.csv"
)

TOLERANCE <- 1e-6

if (!file.exists(MODEL_FILE)) {
  stop("Could not find model file: ", MODEL_FILE, call. = FALSE)
}

if (!file.exists(GROUP_COUNT_FILE)) {
  stop("Could not find group-count file: ", GROUP_COUNT_FILE, call. = FALSE)
}

model_object <- readRDS(MODEL_FILE)
group_counts <- utils::read.csv(
  GROUP_COUNT_FILE,
  stringsAsFactors = FALSE
)

required_model_elements <- c("fits", "symptoms", "pair_index")
missing_model_elements <- setdiff(
  required_model_elements,
  names(model_object)
)

if (length(missing_model_elements) > 0L) {
  stop(
    "The model object is missing: ",
    paste(missing_model_elements, collapse = ", "),
    call. = FALSE
  )
}

if (!"free_J_free_h" %in% names(model_object$fits)) {
  stop("The RDS file does not contain free_J_free_h.", call. = FALSE)
}

fit <- model_object$fits$free_J_free_h

if (is.null(fit$vcov) || !is.matrix(fit$vcov)) {
  stop(
    "The fitted model does not contain a valid covariance matrix. ",
    "Rerun R/04_network_invariance_exact.R with hessian = TRUE.",
    call. = FALSE
  )
}

symptoms <- model_object$symptoms
pair_index <- model_object$pair_index

if (!all(c("n") %in% names(group_counts)) || nrow(group_counts) != 2L) {
  stop(
    "group_counts.csv must contain exactly two rows and an 'n' column.",
    call. = FALSE
  )
}

n_low <- as.numeric(group_counts$n[1])
n_high <- as.numeric(group_counts$n[2])

if (any(!is.finite(c(n_low, n_high))) || any(c(n_low, n_high) <= 0)) {
  stop("Invalid group sample sizes.", call. = FALSE)
}

p <- length(symptoms)
q <- nrow(pair_index)

expected_parameter_count <- 2L * (p + q)

if (length(fit$par) != expected_parameter_count) {
  stop(
    "Unexpected number of parameters: ",
    length(fit$par),
    "; expected ", expected_parameter_count, ".",
    call. = FALSE
  )
}

index_h_low <- seq_len(p)
index_J_low <- p + seq_len(q)

offset <- p + q
index_h_high <- offset + seq_len(p)
index_J_high <- offset + p + seq_len(q)

h_low <- fit$par[index_h_low]
J_low <- fit$par[index_J_low]
h_high <- fit$par[index_h_high]
J_high <- fit$par[index_J_high]

theta_low <- c(h_low, J_low)
theta_high <- c(h_high, J_high)

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

exact_information <- function(theta, n, sufficient_statistics) {
  eta <- as.numeric(sufficient_statistics %*% theta)

  max_eta <- max(eta)
  unnormalized <- exp(eta - max_eta)
  probabilities <- unnormalized / sum(unnormalized)

  expected_statistics <- as.numeric(
    crossprod(probabilities, sufficient_statistics)
  )

  weighted_statistics <-
    sufficient_statistics * sqrt(probabilities)

  expected_outer_product <- crossprod(weighted_statistics)

  covariance_statistics <-
    expected_outer_product -
    tcrossprod(expected_statistics)

  information <- n * covariance_statistics
  (information + t(information)) / 2
}

safe_inverse <- function(x, label) {
  eigenvalues <- eigen(
    x,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  positive_definite <- all(eigenvalues > 0)
  condition_number <- if (positive_definite) {
    max(eigenvalues) / min(eigenvalues)
  } else {
    Inf
  }

  if (!positive_definite) {
    stop(
      "Expected Fisher information is not positive definite for ",
      label,
      ".",
      call. = FALSE
    )
  }

  inverse <- tryCatch(
    solve(x),
    error = function(e) {
      stop(
        "Could not invert expected Fisher information for ",
        label,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  list(
    inverse = inverse,
    condition_number = condition_number,
    min_eigenvalue = min(eigenvalues)
  )
}

information_low <- exact_information(
  theta_low,
  n_low,
  state_statistics
)

information_high <- exact_information(
  theta_high,
  n_high,
  state_statistics
)

inverse_low <- safe_inverse(
  information_low,
  "Low Slow Risk"
)

inverse_high <- safe_inverse(
  information_high,
  "High Slow Risk"
)

variance_h_low_fisher <-
  diag(inverse_low$inverse)[seq_len(p)]

variance_h_high_fisher <-
  diag(inverse_high$inverse)[seq_len(p)]

se_delta_fisher <- sqrt(
  variance_h_low_fisher +
    variance_h_high_fisher
)

contrast_matrix <- matrix(
  0,
  nrow = p,
  ncol = length(fit$par)
)

for (i in seq_len(p)) {
  contrast_matrix[i, index_h_low[i]] <- -1
  contrast_matrix[i, index_h_high[i]] <- 1
}

variance_delta_hessian <- diag(
  contrast_matrix %*%
    fit$vcov %*%
    t(contrast_matrix)
)

se_delta_hessian <- sqrt(
  pmax(variance_delta_hessian, 0)
)

difference <- h_high - h_low

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

result <- data.frame(
  symptom = symptoms,
  symptom_label = unname(symptom_labels[symptoms]),
  h_low = h_low,
  h_high = h_high,
  difference_high_minus_low = difference,
  se_difference_exact_fisher = se_delta_fisher,
  se_difference_saved_hessian = se_delta_hessian,
  absolute_se_difference = abs(
    se_delta_fisher -
      se_delta_hessian
  ),
  stringsAsFactors = FALSE
)

maximum_absolute_se_difference <- max(
  result$absolute_se_difference
)

validation_passed <-
  maximum_absolute_se_difference <= TOLERANCE

utils::write.csv(
  result,
  OUTPUT_FILE,
  row.names = FALSE,
  na = ""
)

diagnostics <- data.frame(
  quantity = c(
    "n_low",
    "n_high",
    "expected_information_condition_number_low",
    "expected_information_condition_number_high",
    "expected_information_min_eigenvalue_low",
    "expected_information_min_eigenvalue_high",
    "maximum_absolute_optimization_gradient",
    "maximum_absolute_se_difference",
    "tolerance",
    "validation_passed"
  ),
  value = c(
    n_low,
    n_high,
    inverse_low$condition_number,
    inverse_high$condition_number,
    inverse_low$min_eigenvalue,
    inverse_high$min_eigenvalue,
    fit$max_abs_gradient,
    maximum_absolute_se_difference,
    TOLERANCE,
    validation_passed
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  diagnostics,
  DIAGNOSTIC_FILE,
  row.names = FALSE,
  na = ""
)

if (!validation_passed) {
  warning(
    "Expected-Fisher and saved-Hessian SEs differ by more than ",
    format(TOLERANCE, scientific = TRUE),
    ". Maximum absolute difference = ",
    format(maximum_absolute_se_difference, scientific = TRUE)
  )
}

message(
  "\nActivation covariance validation complete.\n",
  "Validation passed: ", validation_passed, "\n",
  "Maximum absolute SE difference: ",
  format(maximum_absolute_se_difference, scientific = TRUE), "\n",
  "CSV: ", normalizePath(OUTPUT_FILE, mustWork = FALSE), "\n",
  "Diagnostics: ", normalizePath(DIAGNOSTIC_FILE, mustWork = FALSE)
)

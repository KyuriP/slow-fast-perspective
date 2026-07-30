# ==============================================================================
# R/04_network_invariance_exact.R
#
# Exact multigroup Ising comparison in the manuscript's 0/1 parameterization:
#
#   P(X = x) proportional to
#   exp{ sum_i h_i x_i + sum_{i<j} J_ij x_i x_j }.
#
# This script tests, without a separate beta/inverse-temperature parameter,
# whether Low and High Slow Risk groups differ in:
#
#   1. symptom interactions J,
#   2. activation parameters h,
#   3. both, or
#   4. neither.
#
# Four exact-likelihood models are fitted:
#
#   free_J_free_h   : J_low != J_high, h_low != h_high
#   equal_J_free_h  : J_low == J_high, h_low != h_high
#   free_J_equal_h  : J_low != J_high, h_low == h_high
#   equal_J_equal_h : J_low == J_high, h_low == h_high
#
# Because p = 9, the partition function is evaluated exactly over all
# 2^9 = 512 symptom configurations. No psychonetrics beta parameter is used.
#
# Run from the repository root:
#
#   source("R/04_network_invariance_exact.R")
#
# Outputs:
#
#   results/network_invariance_exact/slow_risk/
#
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

SEED <- 20260729L

# Run the standard Network Comparison Test as a complementary analysis.
# Set FALSE when you already have the final NCT output and only need the exact
# multigroup model comparison.
RUN_NCT <- TRUE

# Start with Low versus High Slow Risk. Set TRUE only after the primary analysis
# has run successfully.
RUN_SECONDARY_COMPARISONS <- FALSE

# NCT settings. Use 1,000 for a working run and 5,000 or more for the final run.
NCT_ITERATIONS <- 1000L
NCT_GAMMA <- 0.25
NCT_EDGE_ADJUSTMENT <- "holm"
NCT_PROGRESS <- TRUE

# Exact-likelihood optimization settings.
OPTIM_MAXIT <- 10000L
OPTIM_RELTOL <- 1e-12
OPTIM_TRACE <- 0L

OUTPUT_ROOT <- file.path("results", "network_invariance_exact")


# ==============================================================================
# 1. PROJECT AND PACKAGE CHECKS
# ==============================================================================

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
    paste(missing_source_files, collapse = ", ")
  )
}

invisible(lapply(source_files, source))

if (!exists("analysis_df_bin", inherits = TRUE)) {
  stop("analysis_df_bin was not created by the preparation scripts.")
}

if (!exists("symptoms", inherits = TRUE)) {
  stop("The symptom vector 'symptoms' was not created by R/00_setup.R.")
}

if (isTRUE(RUN_NCT) &&
    !requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  stop(
    "RUN_NCT is TRUE, but NetworkComparisonTest is not installed.\n",
    "Install it with install.packages('NetworkComparisonTest'), ",
    "or set RUN_NCT <- FALSE."
  )
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)


# ==============================================================================
# 2. OUTPUT HELPERS
# ==============================================================================

write_csv <- function(x, path, row.names = FALSE) {
  utils::write.csv(
    as.data.frame(x),
    file = path,
    row.names = row.names,
    na = ""
  )
  invisible(path)
}

write_matrix_csv <- function(x, path) {
  utils::write.csv(
    as.data.frame(x),
    file = path,
    row.names = TRUE,
    na = ""
  )
  invisible(path)
}

safe_scalar <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  if (length(y) == 0L || !is.finite(y[1])) {
    return(default)
  }
  y[1]
}


# ==============================================================================
# 3. COMPARISON DEFINITIONS
# ==============================================================================

comparison_specs <- list(
  slow_risk = list(
    label = "Low versus High Slow Risk",
    group_var = "slow_group",
    group_levels = c("Low Slow Risk", "High Slow Risk")
  )
)

if (isTRUE(RUN_SECONDARY_COMPARISONS)) {
  comparison_specs$ethnicity <- list(
    label = "Dutch versus Non-Dutch",
    group_var = "dutch_grp",
    group_levels = c("Dutch", "Non-Dutch")
  )

  comparison_specs$age <- list(
    label = "Younger versus Older",
    group_var = "age_grp",
    group_levels = c("Younger", "Older")
  )
}


# ==============================================================================
# 4. DATA PREPARATION
# ==============================================================================

prepare_comparison_data <- function(spec) {
  required_columns <- c(symptoms, spec$group_var)
  missing_columns <- setdiff(required_columns, names(analysis_df_bin))

  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) for ", spec$label, ": ",
      paste(missing_columns, collapse = ", ")
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

  dat <- dat[order(dat[[spec$group_var]]), , drop = FALSE]

  for (v in symptoms) {
    dat[[v]] <- as.integer(dat[[v]])

    bad_values <- setdiff(unique(dat[[v]]), c(0L, 1L))
    if (length(bad_values) > 0L) {
      stop(
        "Symptom ", v, " contains values other than 0 and 1: ",
        paste(bad_values, collapse = ", ")
      )
    }
  }

  group_counts <- table(dat[[spec$group_var]])

  if (length(group_counts) != 2L || any(group_counts == 0L)) {
    stop(
      "Both requested groups must be present for ", spec$label, "."
    )
  }

  for (g in spec$group_levels) {
    group_dat <- dat[
      dat[[spec$group_var]] == g,
      symptoms,
      drop = FALSE
    ]

    constant_variables <- names(which(vapply(
      group_dat,
      function(x) length(unique(x)) < 2L,
      logical(1)
    )))

    if (length(constant_variables) > 0L) {
      stop(
        "Constant symptom(s) in group '", g, "': ",
        paste(constant_variables, collapse = ", ")
      )
    }
  }

  dat
}


# ==============================================================================
# 5. EXACT ISING HELPERS
# ==============================================================================

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

make_state_space <- function(p) {
  as.matrix(
    expand.grid(
      replicate(p, c(0, 1), simplify = FALSE)
    )
  )
}

make_pair_index <- function(p) {
  which(
    upper.tri(matrix(0, p, p)),
    arr.ind = TRUE
  )
}

pair_products <- function(X, pair_index) {
  X[, pair_index[, 1], drop = FALSE] *
    X[, pair_index[, 2], drop = FALSE]
}

make_group_sufficient_statistics <- function(X, pair_index) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  list(
    n = nrow(X),
    observed = c(
      colSums(X),
      colSums(pair_products(X, pair_index))
    ),
    prevalence = colMeans(X)
  )
}

safe_logit_prevalence <- function(prevalence, eps = 0.005) {
  stats::qlogis(
    pmin(
      pmax(prevalence, eps),
      1 - eps
    )
  )
}

group_loglik_and_score <- function(
    h,
    J,
    group_stats,
    state_statistics) {

  theta <- c(h, J)
  eta <- as.numeric(state_statistics %*% theta)

  max_eta <- max(eta)
  exp_eta <- exp(eta - max_eta)
  denominator <- sum(exp_eta)
  probabilities <- exp_eta / denominator

  log_partition <- max_eta + log(denominator)

  expected_statistics <- as.numeric(
    crossprod(probabilities, state_statistics)
  )

  log_likelihood <-
    sum(theta * group_stats$observed) -
    group_stats$n * log_partition

  score <-
    group_stats$observed -
    group_stats$n * expected_statistics

  list(
    log_likelihood = log_likelihood,
    score_h = score[seq_along(h)],
    score_J = score[length(h) + seq_along(J)]
  )
}


# ==============================================================================
# 6. PARAMETER MAPS FOR THE FOUR MODELS
# ==============================================================================

unpack_two_group_parameters <- function(par, model, p, q) {
  if (model == "free_J_free_h") {
    h1 <- par[seq_len(p)]
    J1 <- par[p + seq_len(q)]

    offset <- p + q
    h2 <- par[offset + seq_len(p)]
    J2 <- par[offset + p + seq_len(q)]
  } else if (model == "equal_J_free_h") {
    h1 <- par[seq_len(p)]
    h2 <- par[p + seq_len(p)]
    J_shared <- par[2L * p + seq_len(q)]

    J1 <- J_shared
    J2 <- J_shared
  } else if (model == "free_J_equal_h") {
    h_shared <- par[seq_len(p)]
    J1 <- par[p + seq_len(q)]
    J2 <- par[p + q + seq_len(q)]

    h1 <- h_shared
    h2 <- h_shared
  } else if (model == "equal_J_equal_h") {
    h_shared <- par[seq_len(p)]
    J_shared <- par[p + seq_len(q)]

    h1 <- h_shared
    h2 <- h_shared
    J1 <- J_shared
    J2 <- J_shared
  } else {
    stop("Unknown model: ", model)
  }

  list(
    h1 = h1,
    J1 = J1,
    h2 = h2,
    J2 = J2
  )
}

combine_two_group_scores <- function(
    score_group_1,
    score_group_2,
    model) {

  if (model == "free_J_free_h") {
    score <- c(
      score_group_1$score_h,
      score_group_1$score_J,
      score_group_2$score_h,
      score_group_2$score_J
    )
  } else if (model == "equal_J_free_h") {
    score <- c(
      score_group_1$score_h,
      score_group_2$score_h,
      score_group_1$score_J + score_group_2$score_J
    )
  } else if (model == "free_J_equal_h") {
    score <- c(
      score_group_1$score_h + score_group_2$score_h,
      score_group_1$score_J,
      score_group_2$score_J
    )
  } else if (model == "equal_J_equal_h") {
    score <- c(
      score_group_1$score_h + score_group_2$score_h,
      score_group_1$score_J + score_group_2$score_J
    )
  } else {
    stop("Unknown model: ", model)
  }

  score
}


# ==============================================================================
# 7. EXACT-LIKELIHOOD FITTING
# ==============================================================================

make_objective_functions <- function(
    model,
    group_stats,
    state_statistics,
    p,
    q) {

  objective <- function(par) {
    parameters <- unpack_two_group_parameters(
      par = par,
      model = model,
      p = p,
      q = q
    )

    result_1 <- group_loglik_and_score(
      h = parameters$h1,
      J = parameters$J1,
      group_stats = group_stats[[1]],
      state_statistics = state_statistics
    )

    result_2 <- group_loglik_and_score(
      h = parameters$h2,
      J = parameters$J2,
      group_stats = group_stats[[2]],
      state_statistics = state_statistics
    )

    -(result_1$log_likelihood + result_2$log_likelihood)
  }

  gradient <- function(par) {
    parameters <- unpack_two_group_parameters(
      par = par,
      model = model,
      p = p,
      q = q
    )

    result_1 <- group_loglik_and_score(
      h = parameters$h1,
      J = parameters$J1,
      group_stats = group_stats[[1]],
      state_statistics = state_statistics
    )

    result_2 <- group_loglik_and_score(
      h = parameters$h2,
      J = parameters$J2,
      group_stats = group_stats[[2]],
      state_statistics = state_statistics
    )

    -combine_two_group_scores(
      score_group_1 = result_1,
      score_group_2 = result_2,
      model = model
    )
  }

  list(
    objective = objective,
    gradient = gradient
  )
}

fit_exact_two_group_model <- function(
    model,
    start,
    group_stats,
    state_statistics,
    p,
    q,
    verbose = TRUE) {

  functions <- make_objective_functions(
    model = model,
    group_stats = group_stats,
    state_statistics = state_statistics,
    p = p,
    q = q
  )

  if (isTRUE(verbose)) {
    message("\nFitting exact model: ", model)
  }

  fit <- stats::optim(
    par = start,
    fn = functions$objective,
    gr = functions$gradient,
    method = "BFGS",
    hessian = TRUE,
    control = list(
      maxit = OPTIM_MAXIT,
      reltol = OPTIM_RELTOL,
      trace = OPTIM_TRACE,
      REPORT = 10
    )
  )

  final_gradient <- functions$gradient(fit$par)
  
  # The objective is the negative log-likelihood, so its Hessian is the
  # observed information matrix. Its inverse estimates the covariance matrix.
  
  hessian <- (fit$hessian + t(fit$hessian)) / 2
  
  hessian_eigenvalues <- eigen(
    hessian,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  min_hessian_eigenvalue <- min(hessian_eigenvalues)
  
  hessian_positive_definite <- all(hessian_eigenvalues > 0)
  
  hessian_condition_number <- if (hessian_positive_definite) {
    max(hessian_eigenvalues) / min(hessian_eigenvalues)
  } else {
    Inf
  }
  
  vcov_matrix <- tryCatch(
    solve(hessian),
    error = function(e) {
      warning(
        "Could not invert Hessian for model ",
        model,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (fit$convergence != 0L) {
    warning(
      "Model ", model,
      " returned convergence code ", fit$convergence,
      ": ", fit$message
    )
  }

  list(
    model = model,
    par = fit$par,
    logLik = -fit$value,
    npar = length(fit$par),
    convergence = fit$convergence,
    message = fit$message,
    counts = fit$counts,
    max_abs_gradient = max(abs(final_gradient)),
    mean_abs_gradient = mean(abs(final_gradient)),
    hessian = hessian,
    vcov = vcov_matrix,
    min_hessian_eigenvalue = min_hessian_eigenvalue,
    hessian_positive_definite = hessian_positive_definite,
    hessian_condition_number = hessian_condition_number,
    optim = fit
  )
}


# ==============================================================================
# 8. STARTING VALUES
# ==============================================================================

make_starting_values <- function(group_stats, p, q) {
  h1_initial <- safe_logit_prevalence(group_stats[[1]]$prevalence)
  h2_initial <- safe_logit_prevalence(group_stats[[2]]$prevalence)

  n1 <- group_stats[[1]]$n
  n2 <- group_stats[[2]]$n
  n_total <- n1 + n2

  pooled_prevalence <-
    (
      n1 * group_stats[[1]]$prevalence +
      n2 * group_stats[[2]]$prevalence
    ) / n_total

  h_pooled_initial <- safe_logit_prevalence(pooled_prevalence)

  J_zero <- rep(0, q)

  list(
    free_J_free_h = c(
      h1_initial,
      J_zero,
      h2_initial,
      J_zero
    ),
    equal_J_free_h = c(
      h1_initial,
      h2_initial,
      J_zero
    ),
    free_J_equal_h = c(
      h_pooled_initial,
      J_zero,
      J_zero
    ),
    equal_J_equal_h = c(
      h_pooled_initial,
      J_zero
    )
  )
}

refine_starts_from_fitted_models <- function(
    preliminary_fits,
    p,
    q,
    group_stats) {

  free_parameters <- unpack_two_group_parameters(
    preliminary_fits$free_J_free_h$par,
    model = "free_J_free_h",
    p = p,
    q = q
  )

  pooled_parameters <- unpack_two_group_parameters(
    preliminary_fits$equal_J_equal_h$par,
    model = "equal_J_equal_h",
    p = p,
    q = q
  )

  n1 <- group_stats[[1]]$n
  n2 <- group_stats[[2]]$n
  n_total <- n1 + n2

  weighted_J <-
    (
      n1 * free_parameters$J1 +
      n2 * free_parameters$J2
    ) / n_total

  weighted_h <-
    (
      n1 * free_parameters$h1 +
      n2 * free_parameters$h2
    ) / n_total

  list(
    free_J_free_h = preliminary_fits$free_J_free_h$par,
    equal_J_free_h = c(
      free_parameters$h1,
      free_parameters$h2,
      weighted_J
    ),
    free_J_equal_h = c(
      weighted_h,
      free_parameters$J1,
      free_parameters$J2
    ),
    equal_J_equal_h = c(
      pooled_parameters$h1,
      pooled_parameters$J1
    )
  )
}


# ==============================================================================
# 9. PARAMETER EXPORT
# ==============================================================================

edge_vector_to_matrix <- function(J, pair_index, symptoms) {
  p <- length(symptoms)

  matrix_J <- matrix(0, nrow = p, ncol = p)
  matrix_J[pair_index] <- J
  matrix_J <- matrix_J + t(matrix_J)

  rownames(matrix_J) <- symptoms
  colnames(matrix_J) <- symptoms

  matrix_J
}

export_model_parameters <- function(
    fit,
    spec,
    p,
    q,
    pair_index,
    model_dir) {

  parameters <- unpack_two_group_parameters(
    par = fit$par,
    model = fit$model,
    p = p,
    q = q
  )

  group_labels <- spec$group_levels

  activation_table <- rbind(
    data.frame(
      model = fit$model,
      group = group_labels[1],
      symptom = symptoms,
      h = parameters$h1,
      baseline_probability_when_other_symptoms_zero =
        stats::plogis(parameters$h1),
      stringsAsFactors = FALSE
    ),
    data.frame(
      model = fit$model,
      group = group_labels[2],
      symptom = symptoms,
      h = parameters$h2,
      baseline_probability_when_other_symptoms_zero =
        stats::plogis(parameters$h2),
      stringsAsFactors = FALSE
    )
  )

  edge_table_group_1 <- data.frame(
    model = fit$model,
    group = group_labels[1],
    node_1 = symptoms[pair_index[, 1]],
    node_2 = symptoms[pair_index[, 2]],
    J = parameters$J1,
    stringsAsFactors = FALSE
  )

  edge_table_group_2 <- data.frame(
    model = fit$model,
    group = group_labels[2],
    node_1 = symptoms[pair_index[, 1]],
    node_2 = symptoms[pair_index[, 2]],
    J = parameters$J2,
    stringsAsFactors = FALSE
  )

  edge_table <- rbind(
    edge_table_group_1,
    edge_table_group_2
  )

  activation_difference <- data.frame(
    model = fit$model,
    symptom = symptoms,
    h_group_1 = parameters$h1,
    h_group_2 = parameters$h2,
    h_difference_group2_minus_group1 =
      parameters$h2 - parameters$h1,
    baseline_probability_group_1 =
      stats::plogis(parameters$h1),
    baseline_probability_group_2 =
      stats::plogis(parameters$h2),
    baseline_probability_difference_group2_minus_group1 =
      stats::plogis(parameters$h2) -
      stats::plogis(parameters$h1),
    stringsAsFactors = FALSE
  )

  edge_difference <- data.frame(
    model = fit$model,
    node_1 = symptoms[pair_index[, 1]],
    node_2 = symptoms[pair_index[, 2]],
    J_group_1 = parameters$J1,
    J_group_2 = parameters$J2,
    J_difference_group2_minus_group1 =
      parameters$J2 - parameters$J1,
    absolute_J_difference =
      abs(parameters$J2 - parameters$J1),
    stringsAsFactors = FALSE
  )

  edge_difference <- edge_difference[
    order(edge_difference$absolute_J_difference, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  J_matrix_group_1 <- edge_vector_to_matrix(
    parameters$J1,
    pair_index,
    symptoms
  )

  J_matrix_group_2 <- edge_vector_to_matrix(
    parameters$J2,
    pair_index,
    symptoms
  )

  write_csv(
    activation_table,
    file.path(
      model_dir,
      paste0("activation_parameters_", fit$model, ".csv")
    )
  )

  write_csv(
    edge_table,
    file.path(
      model_dir,
      paste0("interaction_parameters_", fit$model, ".csv")
    )
  )

  write_csv(
    activation_difference,
    file.path(
      model_dir,
      paste0("activation_differences_", fit$model, ".csv")
    )
  )

  write_csv(
    edge_difference,
    file.path(
      model_dir,
      paste0("interaction_differences_", fit$model, ".csv")
    )
  )

  write_matrix_csv(
    J_matrix_group_1,
    file.path(
      model_dir,
      paste0("J_", fit$model, "_group_1.csv")
    )
  )

  write_matrix_csv(
    J_matrix_group_2,
    file.path(
      model_dir,
      paste0("J_", fit$model, "_group_2.csv")
    )
  )

  invisible(
    list(
      activation = activation_table,
      interactions = edge_table,
      activation_difference = activation_difference,
      interaction_difference = edge_difference
    )
  )
}


# ==============================================================================
# 10. MODEL COMPARISON
# ==============================================================================

make_model_comparison_table <- function(fits, n_total) {
  model_order <- c(
    "free_J_free_h",
    "equal_J_free_h",
    "free_J_equal_h",
    "equal_J_equal_h"
  )

  table <- do.call(
    rbind,
    lapply(model_order, function(model_name) {
      fit <- fits[[model_name]]

      data.frame(
        model = model_name,
        npar = fit$npar,
        logLik = fit$logLik,
        deviance = -2 * fit$logLik,
        AIC = 2 * fit$npar - 2 * fit$logLik,
        BIC = log(n_total) * fit$npar - 2 * fit$logLik,
        convergence = fit$convergence,
        max_abs_gradient = fit$max_abs_gradient,
        mean_abs_gradient = fit$mean_abs_gradient,
        stringsAsFactors = FALSE
      )
    })
  )

  free_logLik <- table$logLik[
    table$model == "free_J_free_h"
  ]

  table$logLik_loss_vs_free <-
    free_logLik - table$logLik

  table$logLik_loss_per_participant <-
    table$logLik_loss_vs_free / n_total

  table$delta_AIC <- table$AIC - min(table$AIC)
  table$delta_BIC <- table$BIC - min(table$BIC)

  table
}

likelihood_ratio_test <- function(
    larger_fit,
    smaller_fit,
    comparison_label) {

  statistic <- 2 * (
    larger_fit$logLik -
    smaller_fit$logLik
  )

  degrees_of_freedom <-
    larger_fit$npar -
    smaller_fit$npar

  if (statistic < -1e-5) {
    warning(
      "Restricted model had a higher log likelihood for ",
      comparison_label,
      ". Check convergence."
    )
  }

  statistic <- max(statistic, 0)

  data.frame(
    comparison = comparison_label,
    larger_model = larger_fit$model,
    restricted_model = smaller_fit$model,
    LR_chisq = statistic,
    df = degrees_of_freedom,
    p_value = stats::pchisq(
      statistic,
      df = degrees_of_freedom,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )
}

make_lrt_table <- function(fits) {
  rbind(
    likelihood_ratio_test(
      larger_fit = fits$free_J_free_h,
      smaller_fit = fits$equal_J_free_h,
      comparison_label = "Test equality of J while h remains group-specific"
    ),
    likelihood_ratio_test(
      larger_fit = fits$free_J_free_h,
      smaller_fit = fits$free_J_equal_h,
      comparison_label = "Test equality of h while J remains group-specific"
    ),
    likelihood_ratio_test(
      larger_fit = fits$equal_J_free_h,
      smaller_fit = fits$equal_J_equal_h,
      comparison_label = "Test equality of h given equal J"
    ),
    likelihood_ratio_test(
      larger_fit = fits$free_J_equal_h,
      smaller_fit = fits$equal_J_equal_h,
      comparison_label = "Test equality of J given equal h"
    )
  )
}


# ==============================================================================
# 11. STANDARD NETWORK COMPARISON TEST
# ==============================================================================

run_nct_analysis <- function(dat, spec, comparison_dir) {
  message("\nRunning NCT: ", spec$label)

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

  # Use a direct call. Some NetworkComparisonTest versions can fail when NCT
  # is called through do.call() because the function internally inspects its
  # matched call.
  set.seed(SEED)

  result <- NetworkComparisonTest::NCT(
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
    progressbar = NCT_PROGRESS,
    make.positive.definite = TRUE,
    p.adjust.methods = NCT_EDGE_ADJUSTMENT,
    test.centrality = FALSE,
    verbose = TRUE
  )

  saveRDS(
    result,
    file.path(comparison_dir, "nct_result.rds")
  )

  strengths <- as.numeric(result$glstrinv.sep)
  if (length(strengths) < 2L) {
    strengths <- c(NA_real_, NA_real_)
  }

  summary_table <- data.frame(
    comparison = spec$label,
    group_1 = spec$group_levels[1],
    group_2 = spec$group_levels[2],
    n_group_1 = nrow(group_1),
    n_group_2 = nrow(group_2),
    network_structure_statistic_M =
      safe_scalar(result$nwinv.real),
    network_structure_p =
      safe_scalar(result$nwinv.pval),
    global_strength_group_1 = strengths[1],
    global_strength_group_2 = strengths[2],
    global_strength_difference_S =
      safe_scalar(result$glstrinv.real),
    global_strength_p =
      safe_scalar(result$glstrinv.pval),
    permutations = NCT_ITERATIONS,
    stringsAsFactors = FALSE
  )

  write_csv(
    summary_table,
    file.path(comparison_dir, "nct_summary.csv")
  )

  if (!is.null(result$einv.pvals)) {
    edge_tests <- as.data.frame(result$einv.pvals)
    write_csv(
      edge_tests,
      file.path(comparison_dir, "nct_edge_tests.csv")
    )
  }

  write_matrix_csv(
    result$nw1,
    file.path(comparison_dir, "nct_network_group_1.csv")
  )

  write_matrix_csv(
    result$nw2,
    file.path(comparison_dir, "nct_network_group_2.csv")
  )

  summary_table
}


# ==============================================================================
# 12. RUN ONE COMPARISON
# ==============================================================================

run_exact_comparison <- function(
    comparison_name,
    spec) {

  message(
    "\n",
    paste(rep("=", 78), collapse = ""),
    "\nCOMPARISON: ", spec$label,
    "\n",
    paste(rep("=", 78), collapse = "")
  )

  comparison_dir <- file.path(
    OUTPUT_ROOT,
    comparison_name
  )

  model_dir <- file.path(
    comparison_dir,
    "exact_multigroup_ising"
  )

  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dat <- prepare_comparison_data(spec)

  group_1_data <- as.matrix(
    dat[
      dat[[spec$group_var]] == spec$group_levels[1],
      symptoms,
      drop = FALSE
    ]
  )

  group_2_data <- as.matrix(
    dat[
      dat[[spec$group_var]] == spec$group_levels[2],
      symptoms,
      drop = FALSE
    ]
  )

  storage.mode(group_1_data) <- "numeric"
  storage.mode(group_2_data) <- "numeric"

  group_counts <- data.frame(
    group = spec$group_levels,
    n = c(
      nrow(group_1_data),
      nrow(group_2_data)
    ),
    stringsAsFactors = FALSE
  )

  write_csv(
    group_counts,
    file.path(comparison_dir, "group_counts.csv")
  )

  observed_prevalence <- rbind(
    data.frame(
      group = spec$group_levels[1],
      symptom = symptoms,
      prevalence = colMeans(group_1_data),
      stringsAsFactors = FALSE
    ),
    data.frame(
      group = spec$group_levels[2],
      symptom = symptoms,
      prevalence = colMeans(group_2_data),
      stringsAsFactors = FALSE
    )
  )

  write_csv(
    observed_prevalence,
    file.path(
      comparison_dir,
      "observed_symptom_prevalence.csv"
    )
  )

  p <- length(symptoms)
  pair_index <- make_pair_index(p)
  q <- nrow(pair_index)

  states <- make_state_space(p)
  colnames(states) <- symptoms

  state_pairs <- pair_products(
    states,
    pair_index
  )

  state_statistics <- cbind(
    states,
    state_pairs
  )

  group_stats <- list(
    make_group_sufficient_statistics(
      group_1_data,
      pair_index
    ),
    make_group_sufficient_statistics(
      group_2_data,
      pair_index
    )
  )

  # First pass from neutral, prevalence-based starts.
  initial_starts <- make_starting_values(
    group_stats = group_stats,
    p = p,
    q = q
  )

  preliminary_fits <- lapply(
    names(initial_starts),
    function(model_name) {
      fit_exact_two_group_model(
        model = model_name,
        start = initial_starts[[model_name]],
        group_stats = group_stats,
        state_statistics = state_statistics,
        p = p,
        q = q,
        verbose = TRUE
      )
    }
  )

  names(preliminary_fits) <- names(initial_starts)

  # Second pass uses estimates from the free and pooled models as stronger
  # starting values. Because the Ising log likelihood is concave in its natural
  # parameters, this should converge to the same optimum but improves numerical
  # reliability.
  refined_starts <- refine_starts_from_fitted_models(
    preliminary_fits = preliminary_fits,
    p = p,
    q = q,
    group_stats = group_stats
  )

  fits <- lapply(
    names(refined_starts),
    function(model_name) {
      fit_exact_two_group_model(
        model = model_name,
        start = refined_starts[[model_name]],
        group_stats = group_stats,
        state_statistics = state_statistics,
        p = p,
        q = q,
        verbose = TRUE
      )
    }
  )

  names(fits) <- names(refined_starts)

  n_total <- nrow(group_1_data) + nrow(group_2_data)

  model_comparison <- make_model_comparison_table(
    fits = fits,
    n_total = n_total
  )

  lrt_table <- make_lrt_table(fits)

  write_csv(
    model_comparison,
    file.path(
      model_dir,
      "model_comparison_exact.csv"
    )
  )

  write_csv(
    lrt_table,
    file.path(
      model_dir,
      "nested_likelihood_ratio_tests_exact.csv"
    )
  )

  exported_parameters <- lapply(
    fits,
    function(fit) {
      export_model_parameters(
        fit = fit,
        spec = spec,
        p = p,
        q = q,
        pair_index = pair_index,
        model_dir = model_dir
      )
    }
  )

  saveRDS(
    list(
      specification = spec,
      symptoms = symptoms,
      pair_index = pair_index,
      fits = fits,
      model_comparison = model_comparison,
      likelihood_ratio_tests = lrt_table,
      exported_parameters = exported_parameters
    ),
    file.path(
      model_dir,
      "exact_multigroup_models.rds"
    )
  )

  convergence_summary <- model_comparison[
    ,
    c(
      "model",
      "convergence",
      "max_abs_gradient",
      "mean_abs_gradient"
    ),
    drop = FALSE
  ]

  write_csv(
    convergence_summary,
    file.path(
      model_dir,
      "optimization_diagnostics.csv"
    )
  )

  interpretation_guide <- c(
    paste0("Comparison: ", spec$label),
    "",
    "MODEL DEFINITIONS",
    "free_J_free_h: group-specific interactions and activation parameters.",
    "equal_J_free_h: equal interactions, group-specific activation parameters.",
    "free_J_equal_h: group-specific interactions, equal activation parameters.",
    "equal_J_equal_h: equal interactions and equal activation parameters.",
    "",
    "NESTED TESTS",
    "free_J_free_h vs equal_J_free_h tests exact equality of the 36 J edges.",
    "free_J_free_h vs free_J_equal_h tests exact equality of the 9 h values.",
    "equal_J_free_h vs equal_J_equal_h tests equality of h given equal J.",
    "free_J_equal_h vs equal_J_equal_h tests equality of J given equal h.",
    "",
    "NON-NESTED COMPARISON",
    "equal_J_free_h and free_J_equal_h are not nested.",
    "Compare them using AIC and BIC, not a likelihood-ratio test.",
    "",
    "LARGE-SAMPLE CAUTION",
    "With approximately 24,000 observations, small departures from exact",
    "equality can produce very small p-values. Interpret likelihood-ratio tests",
    "alongside AIC, BIC, per-participant log-likelihood loss, and parameter",
    "magnitudes.",
    "",
    "BASELINE PROBABILITIES",
    "plogis(h_i) is the fitted probability of symptom i when all other symptoms",
    "are set to zero. It is not the observed marginal prevalence.",
    "",
    "THIS SCRIPT USES NO BETA PARAMETER.",
    "Therefore, equality of J and h is tested directly in the same 0/1",
    "parameterization used in the manuscript."
  )

  writeLines(
    interpretation_guide,
    file.path(
      model_dir,
      "INTERPRETATION_GUIDE.txt"
    )
  )

  nct_summary <- NULL

  if (isTRUE(RUN_NCT)) {
    nct_summary <- run_nct_analysis(
      dat = dat,
      spec = spec,
      comparison_dir = comparison_dir
    )
  }

  list(
    comparison = spec,
    group_counts = group_counts,
    model_comparison = model_comparison,
    likelihood_ratio_tests = lrt_table,
    nct_summary = nct_summary,
    fits = fits
  )
}


# ==============================================================================
# 13. RUN REQUESTED COMPARISONS
# ==============================================================================

all_results <- list()

for (comparison_name in names(comparison_specs)) {
  all_results[[comparison_name]] <- run_exact_comparison(
    comparison_name = comparison_name,
    spec = comparison_specs[[comparison_name]]
  )
}

saveRDS(
  all_results,
  file.path(
    OUTPUT_ROOT,
    "all_exact_invariance_results.rds"
  )
)


# ==============================================================================
# 14. REPRODUCIBILITY INFORMATION
# ==============================================================================

settings <- data.frame(
  setting = c(
    "seed",
    "run_nct",
    "run_secondary_comparisons",
    "nct_iterations",
    "nct_gamma",
    "nct_edge_adjustment",
    "optim_maxit",
    "optim_reltol"
  ),
  value = c(
    as.character(SEED),
    as.character(RUN_NCT),
    as.character(RUN_SECONDARY_COMPARISONS),
    as.character(NCT_ITERATIONS),
    as.character(NCT_GAMMA),
    NCT_EDGE_ADJUSTMENT,
    as.character(OPTIM_MAXIT),
    format(OPTIM_RELTOL, scientific = TRUE)
  ),
  stringsAsFactors = FALSE
)

write_csv(
  settings,
  file.path(
    OUTPUT_ROOT,
    "run_settings.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_ROOT,
    "sessionInfo.txt"
  )
)

message(
  "\nExact invariance analysis complete.\nResults written to: ",
  normalizePath(
    OUTPUT_ROOT,
    mustWork = FALSE
  )
)

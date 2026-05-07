# ───────────────────────────────────────────────────────────────────────────────
# JOINT ISING MODEL WITH ONE GLOBAL SHARED J
#
# Model:
#   P(Y_n | X_n) ∝ exp(
#     sum_i h_i(X_n) Y_ni + sum_{i<j} J_ij Y_ni Y_nj
#   )
#
# where:
#   h_i(X_n) = alpha_i
#            + gamma_slow_i * HighSlow_n
#            + gamma_eth_i  * NonDutch_n
#            + gamma_age_i  * Older_n
#
# One J is shared across everyone.
# Baseline activation is allowed to vary with group indicators.
# All parameters are estimated jointly.
# ───────────────────────────────────────────────────────────────────────────────

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

make_state_space <- function(p) {
  as.matrix(expand.grid(replicate(p, c(0, 1), simplify = FALSE)))
}

make_pair_index <- function(p) {
  which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
}

pair_products <- function(X, pair_index) {
  out <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(pair_index))
  for (k in seq_len(nrow(pair_index))) {
    out[, k] <- X[, pair_index[k, 1]] * X[, pair_index[k, 2]]
  }
  out
}

fit_ising_globalJ_covH <- function(df_bin, symptoms, covars,
                                   maxit = 1000, verbose = TRUE) {
  
  dat <- df_bin %>%
    dplyr::select(dplyr::all_of(c(symptoms, covars))) %>%
    tidyr::drop_na()
  
  X <- as.matrix(dat[, symptoms, drop = FALSE])
  storage.mode(X) <- "numeric"
  
  Z <- as.matrix(dat[, covars, drop = FALSE])
  storage.mode(Z) <- "numeric"
  
  p <- length(symptoms)
  q <- length(covars)
  
  pair_index <- make_pair_index(p)
  n_edges <- nrow(pair_index)
  
  states <- make_state_space(p)
  colnames(states) <- symptoms
  state_pairs <- pair_products(states, pair_index)
  
  # Group rows by unique covariate profile.
  # With 3 binary covariates, there are at most 8 profiles.
  profile_id <- apply(Z, 1, paste, collapse = "_")
  profiles <- unique(profile_id)
  
  profile_stats <- lapply(profiles, function(pid) {
    idx <- which(profile_id == pid)
    Xc <- X[idx, , drop = FALSE]
    Zc <- Z[idx[1], , drop = FALSE]
    
    list(
      n = nrow(Xc),
      z = as.numeric(Zc),
      S = colSums(Xc),
      T = colSums(pair_products(Xc, pair_index))
    )
  })
  
  # Parameter vector:
  #   alpha: p baseline external fields
  #   Gamma: p x q covariate effects on external fields
  #   J: n_edges shared pairwise interactions
  #
  # par = c(alpha, as.vector(Gamma), J_edges)
  
  eps <- 0.005
  alpha_init <- qlogis(pmin(pmax(colMeans(X), eps), 1 - eps))
  Gamma_init <- rep(0, p * q)
  J_init <- rep(0, n_edges)
  
  par_init <- c(alpha_init, Gamma_init, J_init)
  
  neg_loglik <- function(par) {
    alpha <- par[1:p]
    
    Gamma_vec <- par[(p + 1):(p + p * q)]
    Gamma <- matrix(Gamma_vec, nrow = p, ncol = q)
    
    J <- par[(p + p * q + 1):(p + p * q + n_edges)]
    
    ll <- 0
    
    for (st in profile_stats) {
      h <- alpha + as.numeric(Gamma %*% st$z)
      
      eta <- as.numeric(states %*% h + state_pairs %*% J)
      logZ <- logsumexp(eta)
      
      ll <- ll +
        sum(h * st$S) +
        sum(J * st$T) -
        st$n * logZ
    }
    
    -ll
  }
  
  fit <- optim(
    par = par_init,
    fn = neg_loglik,
    method = "BFGS",
    control = list(
      maxit = maxit,
      trace = ifelse(verbose, 1, 0),
      REPORT = 10
    )
  )
  
  if (fit$convergence != 0) {
    warning("optim() did not fully converge. convergence code = ", fit$convergence)
  }
  
  par_hat <- fit$par
  
  alpha_hat <- par_hat[1:p]
  names(alpha_hat) <- symptoms
  
  Gamma_vec_hat <- par_hat[(p + 1):(p + p * q)]
  Gamma_hat <- matrix(Gamma_vec_hat, nrow = p, ncol = q)
  rownames(Gamma_hat) <- symptoms
  colnames(Gamma_hat) <- covars
  
  J_edges_hat <- par_hat[(p + p * q + 1):(p + p * q + n_edges)]
  
  J_hat <- matrix(0, p, p)
  J_hat[pair_index] <- J_edges_hat
  J_hat <- J_hat + t(J_hat)
  rownames(J_hat) <- colnames(J_hat) <- symptoms
  
  list(
    alpha = alpha_hat,
    Gamma = Gamma_hat,
    J = J_hat,
    covars = covars,
    logLik = -fit$value,
    npar = length(par_hat),
    nobs = nrow(dat),
    AIC = 2 * length(par_hat) - 2 * (-fit$value),
    BIC = log(nrow(dat)) * length(par_hat) - 2 * (-fit$value),
    optim = fit
  )
}



# ───────────────────────────────────────────────────────────────────────────────
# Prepare covariates for external-field shifts
# Reference profile:
#   Low Slow Risk, Dutch, Younger
# ───────────────────────────────────────────────────────────────────────────────

analysis_df_bin <- analysis_df_bin %>%
  mutate(
    high_slow = as.integer(slow_group == "High Slow Risk"),
    non_dutch = as.integer(dutch_grp == "Non-Dutch"),
    older     = as.integer(age_grp == "Older")
  )

covars <- c("high_slow", "non_dutch", "older")

joint_global <- fit_ising_globalJ_covH(
  df_bin = analysis_df_bin,
  symptoms = symptoms,
  covars = covars,
  maxit = 1000,
  verbose = TRUE
)





# ───────────────────────────────────────────────────────────────────────────────
# Extract adjusted baseline activation probabilities for figure
# Other covariates are held at their sample means.
# ───────────────────────────────────────────────────────────────────────────────

get_baseline_prob <- function(fit, z_profile) {
  h <- fit$alpha + as.numeric(fit$Gamma %*% z_profile)
  tibble::tibble(
    symptom = names(h),
    h = as.numeric(h),
    baseline_prob = plogis(h)
  )
}

z_mean <- colMeans(analysis_df_bin[, covars], na.rm = TRUE)

profiles <- tibble::tibble(
  facet = c(
    "Slow-Risk Group", "Slow-Risk Group",
    "Ethnicity", "Ethnicity",
    "Age Group", "Age Group"
  ),
  group = c(
    "Low Slow Risk", "High Slow Risk",
    "Dutch", "Non-Dutch",
    "Younger", "Older"
  ),
  high_slow = c(0, 1, z_mean["high_slow"], z_mean["high_slow"], z_mean["high_slow"], z_mean["high_slow"]),
  non_dutch = c(z_mean["non_dutch"], z_mean["non_dutch"], 0, 1, z_mean["non_dutch"], z_mean["non_dutch"]),
  older     = c(z_mean["older"], z_mean["older"], z_mean["older"], z_mean["older"], 0, 1)
)

h_long_joint <- purrr::pmap_dfr(
  profiles,
  function(facet, group, high_slow, non_dutch, older) {
    z_profile <- c(
      high_slow = high_slow,
      non_dutch = non_dutch,
      older = older
    )
    
    get_baseline_prob(joint_global, z_profile) %>%
      mutate(
        facet = facet,
        group = group
      )
  }
) %>%
  mutate(
    facet = factor(facet, levels = c("Slow-Risk Group", "Ethnicity", "Age Group")),
    group = factor(group, levels = names(group_colors))
  )







veil_alpha <- 0.6

veil_df <- data.frame(
  facet = factor(
    c("Ethnicity", "Age Group"),
    levels = c("Slow-Risk Group", "Ethnicity", "Age Group")
  )
)

symptom_labels <- c(
  sui = "Suicidality",
  slp = "Sleep problems",
  mot = "Psychomotor change",
  glt = "Guilt",
  ene = "Low energy",
  dep = "Depressed mood",
  con = "Concentration",
  app = "Appetite change",
  anh = "Anhedonia"
)

# Order by average baseline activation in the Slow-Risk comparison
symptom_order_fig3 <- h_long_joint %>%
  filter(facet == "Slow-Risk Group") %>%
  group_by(symptom) %>%
  summarise(mean_prob = mean(baseline_prob), .groups = "drop") %>%
  arrange(mean_prob) %>%
  pull(symptom)

h_long_joint <- h_long_joint %>%
  mutate(symptom = factor(symptom, levels = symptom_order_fig3))

p_baseline_joint <- ggplot(
  h_long_joint,
  aes(x = symptom, y = baseline_prob, color = group, group = group)
) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.5) +
  geom_rect(
    data = veil_df,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "white",
    alpha = veil_alpha,
    inherit.aes = FALSE
  ) +
  coord_flip() +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_x_discrete(labels = symptom_labels) +
  labs(
    title = "Baseline symptom activation by group",
    x = NULL,
    y = "Baseline activation probability",
    color = NULL
  ) +
  theme_paper +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

fig_threshold_validation <- p_baseline_joint

fig_threshold_validation

# ggsave(
#   "figs/fig3_baseline_activation_joint_globalJ.pdf",
#   fig_threshold_validation,
#   width = 8.2,
#   height = 5.4,
#   units = "in"
# )
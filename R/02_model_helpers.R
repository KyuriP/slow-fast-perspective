to_binary01 <- function(v, thr = 1L) {
  if (is.factor(v)) v <- as.character(v)
  v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v), NA_real_, ifelse(v >= thr, 1, 0))
}

circle_layout <- function(p) {
  th <- seq(0, 2*pi, length.out = p + 1L)[-(p + 1L)]
  cbind(cos(th), sin(th))
}

estimate_h_with_offset <- function(df_bin, symptoms, J, G_col = "G") {
  H0  <- setNames(numeric(length(symptoms)), symptoms)
  Gam <- setNames(numeric(length(symptoms)), symptoms)
  
  J <- J[symptoms, symptoms, drop = FALSE]
  
  for (i in seq_along(symptoms)) {
    yi   <- df_bin[[symptoms[i]]]
    nbrs <- setdiff(symptoms, symptoms[i])
    
    eta <- as.numeric(as.matrix(df_bin[, nbrs, drop = FALSE]) %*% J[symptoms[i], nbrs])
    
    mf <- data.frame(y = yi, G = df_bin[[G_col]], off = eta)
    mf <- mf[complete.cases(mf), , drop = FALSE]
    if (nrow(mf) == 0 || length(unique(mf$y)) < 2) next
    
    fit <- glm(y ~ G + offset(off), data = mf, family = binomial())
    co  <- coef(fit)
    
    H0[i]  <- unname(co["(Intercept)"])
    Gam[i] <- if (!is.na(co["G"])) unname(co["G"]) else 0
  }
  
  list(h_level1 = H0, h_level2 = H0 + Gam, gamma = Gam)
}

fit_thresholds_twolevel <- function(df_bin, group_var, level1, level2, symptoms, J) {
  df_g <- df_bin %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(G = as.integer(.data[[group_var]] == level2))
  
  estimate_h_with_offset(df_g, symptoms, J, G_col = "G")
}

simulate_prev <- function(J, h_vec, nsim = 40000, symptoms) {
  J <- J[symptoms, symptoms, drop = FALSE]
  X <- IsingSampler::IsingSampler(nsim, graph = J, thresholds = h_vec, method = "MH")
  setNames(colMeans(X), symptoms)
}

estimate_J <- function(X_df) {
  bootnet::estimateNetwork(X_df, default = "IsingFit")$graph
}

make_binary <- function(x) {
  as.integer(x > 0)
}
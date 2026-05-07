# ───────────────────────────────────────────────────────────────────────────────
# 03_fit_shared_model.R
# Binarize symptoms, estimate shared J, and fit group-specific thresholds
# ───────────────────────────────────────────────────────────────────────────────

thr <- 1L

analysis_df_bin <- analysis_df
analysis_df_bin[symptoms] <- lapply(analysis_df_bin[symptoms], to_binary01, thr = thr)

distinct_counts <- vapply(
  analysis_df_bin[symptoms],
  function(v) length(unique(v[!is.na(v)])),
  integer(1)
)

stopifnot(all(distinct_counts <= 2))

dat_all <- analysis_df_bin %>%
  select(all_of(symptoms)) %>%
  drop_na()

J_shared <- estimate_J(dat_all)
J_shared <- J_shared[symptoms, symptoms, drop = FALSE]

h_eth <- fit_thresholds_twolevel(
  analysis_df_bin,
  "dutch_grp",
  "Dutch",
  "Non-Dutch",
  symptoms,
  J_shared
)

h_age <- fit_thresholds_twolevel(
  analysis_df_bin,
  "age_grp",
  "Younger",
  "Older",
  symptoms,
  J_shared
)

h_slow <- fit_thresholds_twolevel(
  analysis_df_bin,
  "slow_group",
  "Low Slow Risk",
  "High Slow Risk",
  symptoms,
  J_shared
)
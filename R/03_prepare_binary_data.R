# ───────────────────────────────────────────────────────────────────────────────
# 03_prepare_binary_data.R
# Binarize PHQ-9 symptoms for the joint Ising model
# ───────────────────────────────────────────────────────────────────────────────

thr <- 1L

analysis_df_bin <- analysis_df

analysis_df_bin[symptoms] <- lapply(
  analysis_df_bin[symptoms],
  to_binary01,
  thr = thr
)

distinct_counts <- vapply(
  analysis_df_bin[symptoms],
  function(v) length(unique(v[!is.na(v)])),
  integer(1)
)

stopifnot(all(distinct_counts <= 2))
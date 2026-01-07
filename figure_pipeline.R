# ───────────────────────────────────────────────────────────────────────────────
# FIGURE PIPELINE
# Purpose:
#   1) Build "slow layer" composite (Slow Risk Load; SRL)
#   2) Estimate a shared Ising coupling matrix J (pooled data)
#   3) Estimate group-specific thresholds h_g under fixed J (offset GLMs)
#   4) Generate paper figures + summary tables
#
# Notes:
#   - The HELIUS dataset is not included. This script expects: data/helius.sav
# ───────────────────────────────────────────────────────────────────────────────

# necessary packages
suppressPackageStartupMessages({
  # data wrangling + plotting
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(ggrepel)
  
  # network / Ising tooling
  library(bootnet)       # estimateNetwork(..., default="IsingFit")
  library(IsingSampler)  # simulate from Ising(J, h)
  library(qgraph)        # network visualisations
  
  # I/O + summaries used later in the script
  library(haven)         # read_sav
  library(skimr)         # skim()
  library(broom)         # tidy(glm)
})

# set the seed
set.seed(1)


# ───────────────────────────────────────────────────────────────────────────────
# 1) LOAD DATA
# ───────────────────────────────────────────────────────────────────────────────
# HELIUS data are in SPSS format; read_sav preserves labelled vectors.
helius <- haven::read_sav("data/helius.sav")

# ───────────────────────────────────────────────────────────────────────────────
# 2) BUILD SLOW-LAYER INPUTS (candidate “structural burden” indicators)
# ───────────────────────────────────────────────────────────────────────────────
# We treat the slow layer as a relatively slowly varying contextual burden proxy.
# Construction is intentionally transparent (z-score inputs, domain means, composite mean).
slow_vars <- helius %>%
  select(
    gender      = H1_geslacht,
    age         = H1_lft,
    inc_dif     = H1_InkHhMoeite,
    work_sit    = H1_WerkSit,
    edu_cat     = H1_Opleid,
    arbeid      = H1_Arbeidsparticipatie,
    occ_level   = H1_BeroepsNiveau,
    pcs         = H1_PCS12,
    activity    = H1_Squash_totmwk,
    smoking     = H1_PackYears,
    alcohol     = H1_AlcoholConsumption,
    health_lit  = H1_SBSQ_meanscore,   # SBSQ = health literacy
    soc_support = H1_SSQT,             # SSQT = social support
    discrim     = H1_Discr_sumscore,
    stress_work = H1_StressWerk,
    stress_home = H1_StressHuish,
    lang_diff   = H1_MoeiteNLtaal,
    cult_dist   = H1_CultDistMeanScore6_corrected,
    ethn        = H1_etniciteit
  ) %>%
  # Step 2.1: HELIUS missing codes → NA
  # (-9 and -1 are “special missing” codes in the dataset)
  mutate(
    across(everything(), ~ na_if(.x, -9)),
    across(everything(), ~ na_if(.x, -1))
  ) %>%
  # Step 2.2: harmonise some fields + cohort-specific recodes
  mutate(
    # culture distance: ensure numeric, discard invalid values
    cult_dist   = as.numeric(cult_dist),
    cult_dist   = ifelse(cult_dist < 0, NA, cult_dist),
    
    # stress items are intended to be 1–4; everything else becomes NA
    stress_home = ifelse(stress_home %in% 1:4, stress_home, NA),
    stress_work = ifelse(stress_work %in% 1:4, stress_work, NA),
    
    # ethnicity code: 6 treated as missing here
    ethn  = as.numeric(ethn),
    ethn  = ifelse(ethn == 6, NA, ethn),
    
    # Dutch participants (ethn==1): language difficulty / cultural distance are conceptually
    # “not applicable”. Here, missing is coded as 0 = “no difficulty”.
    cult_dist = ifelse(ethn == 1 & is.na(cult_dist), 0, cult_dist),
    lang_diff = ifelse(ethn == 1 & is.na(lang_diff), 0, lang_diff)
  ) %>%
  # Step 2.3: z-score standardize all slow-layer indicators that enter the composite.
  # This places indicators on comparable scales before domain aggregation.
  mutate(
    across(
      c(inc_dif, work_sit, edu_cat, arbeid, occ_level,
        pcs, activity, smoking, alcohol,
        health_lit, soc_support, discrim,
        stress_work, stress_home, lang_diff, cult_dist),
      \(x) as.numeric(scale(x))
    )
  ) %>%
  # Step 2.4: reverse "protective" indicators so higher = higher risk.
  # This ensures all contributions point in the same conceptual direction.
  mutate(
    edu_rev      = -edu_cat,
    occ_rev      = -occ_level,
    pcs_rev      = -pcs,
    activity_rev = -activity,
    lit_rev      = -health_lit,
    support_rev  = -soc_support
  )

# ───────────────────────────────────────────────────────────────────────────────
# 3) SYMPTOMS (PHQ-9 items + sum score)
# ───────────────────────────────────────────────────────────────────────────────
# Note: items are ordinal 0–3 in HELIUS. Later, they are dichotomized for Ising analyses.
symptoms_df <- helius %>%
  select(
    anh    = H1_WlbvRecent1,
    dep    = H1_WlbvRecent2,
    slp    = H1_WlbvRecent3,
    ene    = H1_WlbvRecent4,
    app    = H1_WlbvRecent5,
    glt    = H1_WlbvRecent6,
    con    = H1_WlbvRecent7,
    mot    = H1_WlbvRecent_89,
    sui    = H1_WlbvRecent10,
    PHQsum = H1_PHQ9_sumscore
  ) %>%
  mutate(
    across(everything(), ~ na_if(.x, -9)),
    across(everything(), ~ na_if(.x, -1))
  )

# ───────────────────────────────────────────────────────────────────────────────
# 4) CONSTRUCT Slow Risk Load (3-domain composite)
# ───────────────────────────────────────────────────────────────────────────────
# Domain scores are rowMeans of z-scored indicators; composite is mean of domains.
# Rationale:
#   - simple, transparent, and yields one parsimonious slow-layer dimension.
slow_vars <- slow_vars %>%
  mutate(
    # Socioeconomic/resources domain
    SES_risk = rowMeans(cbind(
      inc_dif, work_sit, arbeid, edu_rev, occ_rev, lit_rev, support_rev
    ), na.rm = TRUE),
    
    # Structural / psychosocial stress domain
    structural_stress = rowMeans(cbind(
      discrim, lang_diff, cult_dist, stress_work, stress_home
    ), na.rm = TRUE),
    
    # Health / lifestyle domain
    health_risk = rowMeans(cbind(
      pcs_rev, activity_rev, smoking, alcohol
    ), na.rm = TRUE),
    
    # Composite slow-layer burden (the "slow layer" proxy; not yet z-scored)
    slow_risk_load = rowMeans(cbind(
      SES_risk, structural_stress, health_risk
    ), na.rm = TRUE)
  )

# ───────────────────────────────────────────────────────────────────────────────
# 5) BUILD ANALYSIS DATAFRAME (MAIN WORKING TABLE)
# ───────────────────────────────────────────────────────────────────────────────
analysis_df <- bind_cols(slow_vars, symptoms_df) %>%
  filter(!is.na(slow_risk_load), !is.na(PHQsum)) %>%   # ensure key variables exist
  mutate(
    # Ethnicity grouping used throughout: Dutch (ethn==1) vs Non-Dutch (everything else not missing)
    dutch_grp = ifelse(ethn == 1, "Dutch", "Non-Dutch"),
    
    # Gender coding for covariate adjustment
    gender_grp = dplyr::case_when(
      gender == 1 ~ "Male",
      gender == 2 ~ "Female",
      TRUE        ~ NA_character_
    ),
    
    # Age split: median split (Older vs Younger)
    age_grp = ifelse(
      age >= median(age, na.rm = TRUE),
      "Older", "Younger"
    ),
    
    # Slow-risk split: median split (High vs Low)
    slow_group = ifelse(
      slow_risk_load > median(slow_risk_load, na.rm = TRUE),
      "High Slow Risk", "Low Slow Risk"
    )
  ) %>%
  mutate(
    # Z-score SRL for any figure/text that uses “(z)” or “per 1 SD” interpretations.
    slow_risk_z = as.numeric(scale(slow_risk_load))
  )

# Quick sanity check on constructed variables (optional)
skimr::skim(analysis_df[, c("slow_risk_load", "slow_risk_z", "PHQsum",
                            "SES_risk", "structural_stress", "health_risk")])

symptom_vars <- c("anh","dep","slp","ene","app","glt","con","mot","sui")

# ───────────────────────────────────────────────────────────────────────────────
# 5) GLOBAL SETTINGS (palette, theme, symptom order)
# ───────────────────────────────────────────────────────────────────────────────
group_colors <- c(
  # Slow Risk
  "Low Slow Risk"  = "#1B9E77",
  "High Slow Risk" = "#D95F02",
  # Ethnicity
  "Dutch"     = "#7570B3",
  "Non-Dutch" = "#E7298A",
  # Age
  "Younger" = "#66A61E",
  "Older"   = "#E6AB02"
)

theme_paper <- theme_minimal(base_size = 13, base_family = "Palatino") +
  theme(panel.grid.minor = element_blank())

# IMPORTANT: From here on, `symptoms` refers to a character vector of symptom names.
symptoms <- symptom_vars

# Preferred order in plots (top-to-bottom after coord_flip)
symptom_order <- c("sui","slp","mot","glt","ene","dep","con","app","anh")
symptoms_plot <- if (all(symptom_order %in% symptoms)) symptom_order else symptoms

# ───────────────────────────────────────────────────────────────────────────────
# 6) HELPERS (binarization, Ising threshold fitting, simulation, network estimation)
# ───────────────────────────────────────────────────────────────────────────────

# Convert ordinal 0–3 → binary 0/1 (thr=1 means "any endorsement" counts as 1)
to_binary01 <- function(v, thr = 1L) {
  if (is.factor(v)) v <- as.character(v)
  v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v), NA_real_, ifelse(v >= thr, 1, 0))
}

# Fixed circular layout for qgraph panels (keeps node positions stable across plots)
circle_layout <- function(p) {
  th <- seq(0, 2*pi, length.out = p + 1L)[-(p + 1L)]
  cbind(cos(th), sin(th))
}

# Estimate group-specific thresholds h_g under a FIXED coupling matrix J via offset-GLMs.
#
# Background (Ising conditionals for binary symptoms Y_i ∈ {0,1}):
#   logit P(Y_i = 1 | Y_-i) = h_i + Σ_{j≠i} J_{ij} Y_j
#
# Constraint imposed here:
#   - J is shared across groups (estimated once from pooled data)
#   - thresholds differ by group:
#       h_i(level2) = h_i(level1) + γ_i
#
# Implementation detail:
#   - For node i, compute η_i = Σ_{j≠i} J_{ij} Y_j and include it as an OFFSET (fixed coefficient = 1)
#   - Fit y ~ G + offset(η_i), where G = 1{group == level2}
#   - intercept = h_i(level1); coef(G) = γ_i; therefore h_i(level2) = intercept + coef(G)
estimate_h_with_offset <- function(df_bin, symptoms, J, G_col = "G") {
  H0  <- setNames(numeric(length(symptoms)), symptoms)  # thresholds for G=0 (level1)
  Gam <- setNames(numeric(length(symptoms)), symptoms)  # threshold shift for G=1 (level2 vs level1)
  
  # Ensure J indexing works and matches symptom order
  J <- J[symptoms, symptoms, drop = FALSE]
  
  for (i in seq_along(symptoms)) {
    yi   <- df_bin[[symptoms[i]]]
    nbrs <- setdiff(symptoms, symptoms[i])
    
    # Offset term = fixed network contribution from all other symptoms
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

# Wrapper for any two-level grouping variable:
#   - encodes group membership into G (0/1)
#   - returns thresholds for level1 (G=0) and level2 (G=1)
fit_thresholds_twolevel <- function(df_bin, group_var, level1, level2, symptoms, J) {
  df_g <- df_bin %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(G = as.integer(.data[[group_var]] == level2))
  
  estimate_h_with_offset(df_g, symptoms, J, G_col = "G")
}

# Simulate from an Ising model with given J and thresholds h, then compute prevalences.
# Note: This checks marginal base rates, not conditional associations (edges).
simulate_prev <- function(J, h_vec, nsim = 40000, symptoms) {
  J <- J[symptoms, symptoms, drop = FALSE]
  X <- IsingSampler::IsingSampler(nsim, graph = J, thresholds = h_vec, method = "MH")
  setNames(colMeans(X), symptoms)
}

# Re-estimate an Ising network using the same estimator everywhere (important for comparability).
# NOTE: IsingFit is penalized; tiny Sim−Emp discrepancies can reflect sparsity/penalization.
estimate_J <- function(X_df) {
  bootnet::estimateNetwork(X_df, default = "IsingFit")$graph
}


# ───────────────────────────────────────────────────────────────────────────────
# 7) ENSURE GROUP VARIABLES + FACTOR LEVELS (for consistent plotting)
# ───────────────────────────────────────────────────────────────────────────────
analysis_df <- analysis_df %>%
  mutate(
    slow_group = factor(slow_group, levels = c("Low Slow Risk","High Slow Risk")),
    dutch_grp  = factor(as.character(dutch_grp), levels = c("Dutch","Non-Dutch")),
    age_grp    = factor(as.character(age_grp),   levels = c("Younger","Older"))
  )


# ───────────────────────────────────────────────────────────────────────────────
# 8) FIG 2: Slow Risk Load distributions + PHQ-9 means (with 99% CI)
# ───────────────────────────────────────────────────────────────────────────────

# Simple x-axis labels (print only; underlying values stay the same)
xlab_group <- function(x) dplyr::recode(
  x,
  "Low Slow Risk"  = "Low Risk",
  "High Slow Risk" = "High Risk",
  .default = x
)

# Long table for faceting (SRL is plotted as z here, consistent with labels/captions)
slow_facet_df <- analysis_df %>%
  select(slow_risk_z, slow_group, dutch_grp, age_grp) %>%
  pivot_longer(c(slow_group, dutch_grp, age_grp),
               names_to = "group_type", values_to = "group") %>%
  filter(!is.na(group)) %>%
  mutate(
    group_type = factor(group_type,
                        levels = c("slow_group","dutch_grp","age_grp"),
                        labels = c("Slow-Risk Group","Ethnicity","Age Group")),
    group = factor(as.character(group), levels = names(group_colors))
  )

# Top panel: boxplots (FREE x-scale per facet so only relevant labels show)
p_top <- ggplot(slow_facet_df, aes(x = group, y = slow_risk_z, fill = group)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.20) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Slow Risk Load (z)") +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(strip.text = element_text(face = "bold", size = 13),
        plot.margin = margin(b = 2))

# Bottom panel: PHQ mean + 99% CI
phq_long <- analysis_df %>%
  select(PHQsum, slow_group, dutch_grp, age_grp) %>%
  pivot_longer(c(slow_group, dutch_grp, age_grp),
               names_to = "group_type", values_to = "group") %>%
  filter(!is.na(PHQsum), !is.na(group)) %>%
  mutate(
    group_type = factor(group_type,
                        levels = c("slow_group","dutch_grp","age_grp"),
                        labels = c("Slow-Risk Group","Ethnicity","Age Group")),
    group = factor(as.character(group), levels = names(group_colors))
  )

# 99% CI uses z_{0.995} ≈ 2.5758 (two-sided)
z99 <- qnorm(0.995)

phq_sum <- phq_long %>%
  group_by(group_type, group) %>%
  summarise(
    n  = dplyr::n(),
    m  = mean(PHQsum, na.rm = TRUE),
    se = sd(PHQsum,  na.rm = TRUE) / sqrt(n),
    lo = m - z99 * se,
    hi = m + z99 * se,
    .groups = "drop"
  )

p_bottom <- ggplot(phq_sum, aes(x = group, y = m, color = group)) +
  geom_point(size = 2.8, shape=1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.8) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "PHQ-9 mean (99% CI)") +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_color_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank(),
    plot.margin  = margin(t = 2)
  )

fig2 <- p_top / p_bottom + plot_layout(heights = c(1, 1))

# ggsave("figs/fig2_slow_top_phq_bottom.pdf", fig2, width = 7, height = 9, units = "in")


# ───────────────────────────────────────────────────────────────────────────────
# 9) BINARIZE SYMPTOMS + ESTIMATE SHARED J (pooled fast-layer coupling structure)
# ───────────────────────────────────────────────────────────────────────────────

thr <- 1L  # 1 = any endorsement (>=1 becomes 1)

analysis_df_bin <- analysis_df
analysis_df_bin[symptoms] <- lapply(analysis_df_bin[symptoms], to_binary01, thr = thr)

# Sanity check: each symptom must be binary (0/1) once missing is removed
distinct_counts <- vapply(
  analysis_df_bin[symptoms],
  function(v) length(unique(v[!is.na(v)])),
  integer(1)
)
stopifnot(all(distinct_counts <= 2))

# Pooled symptom matrix (complete cases) used to estimate the shared coupling matrix J
dat_all <- analysis_df_bin %>% select(all_of(symptoms)) %>% drop_na()

# Shared J estimated once, using eLasso IsingFit in bootnet
J_shared <- estimate_J(dat_all)
J_shared <- J_shared[symptoms, symptoms, drop = FALSE]  # enforce same ordering


# ───────────────────────────────────────────────────────────────────────────────
# 10) FIG 3: (a) group thresholds h_g under fixed J, (b) observed vs simulated prevalence
# ───────────────────────────────────────────────────────────────────────────────

# Step 3.1: estimate thresholds for each grouping (two-level comparisons)
h_eth  <- fit_thresholds_twolevel(analysis_df_bin, "dutch_grp",  "Dutch","Non-Dutch",              symptoms, J_shared)
h_age  <- fit_thresholds_twolevel(analysis_df_bin, "age_grp",    "Younger","Older",               symptoms, J_shared)
h_slow <- fit_thresholds_twolevel(analysis_df_bin, "slow_group", "Low Slow Risk","High Slow Risk", symptoms, J_shared)

# helper to stack thresholds into a tidy long format for ggplot
as_long_h <- function(h_obj, g1, g2, facet_label) {
  tibble(
    symptom = names(h_obj$h_level1),
    !!g1 := as.numeric(h_obj$h_level1),
    !!g2 := as.numeric(h_obj$h_level2)
  ) %>%
    pivot_longer(cols = c(!!g1, !!g2), names_to = "group", values_to = "h") %>%
    mutate(facet = facet_label)
}

# Combine all groupings into one thresholds table for faceting
h_long <- bind_rows(
  as_long_h(h_slow, "Low Slow Risk","High Slow Risk", "Slow-Risk Group"),
  as_long_h(h_eth,  "Dutch","Non-Dutch",              "Ethnicity"),
  as_long_h(h_age,  "Younger","Older",                "Age Group")
) %>%
  mutate(
    facet   = factor(facet, levels = c("Slow-Risk Group","Ethnicity","Age Group")),
    symptom = factor(symptom, levels = symptoms_plot),
    group   = factor(group, levels = names(group_colors))
  )

# Panel (a): thresholds
p_h_facets <- ggplot(h_long,
                     aes(x = fct_rev(symptom), y = h, color = group, group = group)) +
  geom_point(size = 2.1) +
  geom_line(linewidth = 0.7) +
  coord_flip() +
  facet_wrap(~ facet, nrow = 1) +
  labs(title = "Ising thresholds (h) by group and grouping", x = NULL, y = NULL, color = NULL) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_paper +
  theme(strip.text = element_text(face = "bold", size = 13),
        legend.position = "bottom",
        legend.text = element_text(size = 12))

# Panel (b): observed vs simulated prevalences, by group and grouping
# This is a marginal/base-rate validation: points near y=x indicate good recovery of prevalence profiles.
one_group_prev <- function(df_bin, group_var, level1, level2, h_obj, symptoms, J) {
  emp1 <- df_bin %>% filter(.data[[group_var]] == level1) %>%
    select(all_of(symptoms)) %>% drop_na()
  emp2 <- df_bin %>% filter(.data[[group_var]] == level2) %>%
    select(all_of(symptoms)) %>% drop_na()
  
  tibble(
    group_type = group_var,
    group      = rep(c(level1, level2), each = length(symptoms)),
    symptom    = rep(symptoms, times = 2),
    prev_emp   = c(colMeans(emp1), colMeans(emp2)),
    prev_sim   = c(simulate_prev(J, h_obj$h_level1, symptoms = symptoms),
                   simulate_prev(J, h_obj$h_level2, symptoms = symptoms))
  )
}

cmp_prev_all <- bind_rows(
  one_group_prev(analysis_df_bin, "slow_group", "Low Slow Risk","High Slow Risk", h_slow, symptoms, J_shared),
  one_group_prev(analysis_df_bin, "dutch_grp",  "Dutch","Non-Dutch",             h_eth,  symptoms, J_shared),
  one_group_prev(analysis_df_bin, "age_grp",    "Younger","Older",               h_age,  symptoms, J_shared)
) %>%
  mutate(
    group_type = recode(group_type,
                        slow_group = "Slow-Risk Group",
                        dutch_grp  = "Ethnicity",
                        age_grp    = "Age Group"),
    group_type = factor(group_type, levels = c("Slow-Risk Group","Ethnicity","Age Group")),
    group      = factor(group, levels = names(group_colors))
  )

p_prev_all <- ggplot(cmp_prev_all, aes(prev_emp, prev_sim, color = group)) +
  geom_abline(linetype = 2, color = "grey40") +
  geom_point(size = 1, alpha = 0.55) +
  coord_equal() +
  facet_wrap(~ group_type, nrow = 1, scales = "fixed") +
  labs(title = "Observed vs. simulated prevalence", x = NULL, y = NULL, color = NULL) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  theme_paper +
  theme(strip.text = element_blank(),
        legend.position = "bottom")

# Optional: label each point by symptom
p_prev_all_lab <- p_prev_all +
  ggrepel::geom_text_repel(aes(label = symptom),
                           size = 3, show.legend = FALSE, max.overlaps = 100)

fig3 <- ggpubr::ggarrange(
  p_h_facets, p_prev_all_lab,
  legend = "bottom", common.legend = TRUE,
  nrow = 2, labels = c("(a)", "(b)"),
  align = "hv",
  font.label = list(family = "Palatino")
)

# ggsave("figs/fig3_thresholds_prevalence.pdf", fig3, width = 8, height = 8, units = "in")


# ───────────────────────────────────────────────────────────────────────────────
# 11) FIG 4: NETWORK RECOVERY MAPS (Sim − Emp), 6 panels
# ───────────────────────────────────────────────────────────────────────────────
# Edge-level diagnostic:
#   Even if marginal prevalences match, conditional associations (edges) might differ.
#   Here we test whether (shared J + group-specific h) can reproduce each group's
#   estimated edge profile when we simulate under the model and re-estimate networks.
#
# Procedure per group:
#   1) Fit h_g under fixed J_shared
#   2) Simulate data under (J_shared, h_g)
#   3) Re-estimate Ising networks from empirical and simulated data using the SAME estimator
#   4) Plot Δ = J_hat_sim − J_hat_emp

estimate_panel_delta <- function(df_bin, group_var, level1, level2, symptoms, J, n_sim = 50000) {
  # thresholds for the two groups under fixed J
  h_obj <- fit_thresholds_twolevel(df_bin, group_var, level1, level2, symptoms, J)
  
  # empirical datasets (complete cases on symptom set)
  emp1 <- df_bin %>% filter(.data[[group_var]] == level1) %>%
    select(all_of(symptoms)) %>% drop_na()
  emp2 <- df_bin %>% filter(.data[[group_var]] == level2) %>%
    select(all_of(symptoms)) %>% drop_na()
  
  # simulated datasets under the model (shared J, group thresholds)
  sim1 <- IsingSampler::IsingSampler(n_sim, graph = J, thresholds = h_obj$h_level1, method = "MH")
  sim2 <- IsingSampler::IsingSampler(n_sim, graph = J, thresholds = h_obj$h_level2, method = "MH")
  colnames(sim1) <- symptoms
  colnames(sim2) <- symptoms
  
  # re-estimate networks from empirical and simulated data with the SAME estimator
  G_emp1 <- estimate_J(emp1)
  G_emp2 <- estimate_J(emp2)
  G_sim1 <- estimate_J(as.data.frame(sim1))
  G_sim2 <- estimate_J(as.data.frame(sim2))
  
  list(
    level1 = level1, D1 = (G_sim1 - G_emp1),
    level2 = level2, D2 = (G_sim2 - G_emp2)
  )
}

# Compute deltas for the three groupings (each yields 2 panels)
d_slow <- estimate_panel_delta(analysis_df_bin, "slow_group", "Low Slow Risk","High Slow Risk", symptoms, J_shared)
d_eth  <- estimate_panel_delta(analysis_df_bin, "dutch_grp",  "Dutch","Non-Dutch",             symptoms, J_shared)
d_age  <- estimate_panel_delta(analysis_df_bin, "age_grp",    "Younger","Older",               symptoms, J_shared)

plot_six_deltas <- function(d_slow, d_eth, d_age, symptoms, trim = 0) {
  # fixed node coordinates (stable visual comparison across panels)
  lay <- circle_layout(length(symptoms))
  
  # optional: hide tiny differences below 'trim' to reduce clutter
  trim_mat <- function(M) { M[abs(M) < trim] <- 0; M }
  if (trim > 0) {
    d_slow$D1 <- trim_mat(d_slow$D1); d_slow$D2 <- trim_mat(d_slow$D2)
    d_eth$D1  <- trim_mat(d_eth$D1);  d_eth$D2  <- trim_mat(d_eth$D2)
    d_age$D1  <- trim_mat(d_age$D1);  d_age$D2  <- trim_mat(d_age$D2)
  }
  
  # common edge scale across all six panels
  global_max <- max(abs(c(d_slow$D1, d_slow$D2, d_eth$D1, d_eth$D2, d_age$D1, d_age$D2)), na.rm = TRUE)
  
  plot_one <- function(M, title_txt) {
    # Edge colors by sign of Δ: red = simulated edge stronger than empirical (Δ > 0),
    #                            blue = simulated edge weaker than empirical (Δ < 0).
    edge_col <- ifelse(M > 0, "firebrick3", "steelblue3")
    edge_col[is.na(edge_col)] <- NA
    
    qgraph::qgraph(
      M,
      layout = lay,
      labels = symptoms,
      edge.color = edge_col,
      color = "white",
      vsize = 13,
      esize = 12,
      label.cex = 1.2,
      minimum = 0,
      maximum = global_max + 0.2,
      title = NULL
    )
    
    title(
      main = title_txt,
      adj  = 0.5,
      cex.main = 2.8,
      font.main = 2,
      family = "Palatino"
    )
  }
  
  # base plotting for a 2x3 grid of qgraphs
  oldpar <- par(mfcol = c(2,3), mar = c(1,1,3,1), oma = c(0,0,4,0), family = "Palatino")
  on.exit(par(oldpar), add = TRUE)
  
  plot_one(d_slow$D1, paste0("Slow-Risk: ", d_slow$level1))
  plot_one(d_slow$D2, paste0("Slow-Risk: ", d_slow$level2))
  plot_one(d_eth$D1,  paste0("Ethnicity: ", d_eth$level1))
  plot_one(d_eth$D2,  paste0("Ethnicity: ", d_eth$level2))
  plot_one(d_age$D1,  paste0("Age: ", d_age$level1))
  plot_one(d_age$D2,  paste0("Age: ", d_age$level2))
  
  # common title
  mtext("Network recovery (Sim - Emp)", side = 3, outer = TRUE, line = 1, cex = 2.2, font = 2)
}

# pdf("figs/fig4_network_recovery.pdf", width = 16, height = 12, family = "Palatino")
plot_six_deltas(d_slow, d_eth, d_age, symptoms = symptoms, trim = 0)
# dev.off()

# ───────────────────────────────────────────────────────────────────────────────
# 12) FIG 5: (a) Symptom sensitivity (log-odds) vs SRL
#        (b) Strength centrality vs Spearman association with SRL
# ───────────────────────────────────────────────────────────────────────────────

# Define binary endorsement consistent with the dichotomization logic:
# here "1" means any endorsement (>0 on the original 0–3 scale)
make_binary <- function(x) as.integer(x > 0)

# Symptom-wise logistic regressions:
#   make_binary(symptom) ~ SRL(z) + age + gender + ethnicity
# Note: using slow_risk_z ensures the coefficient is interpretable as “per 1 SD SRL”.
symptom_reg <- purrr::map_dfr(symptom_vars, function(sym) {
  f <- as.formula(paste0("make_binary(", sym, ") ~ slow_risk_z + age + gender_grp + dutch_grp"))
  m <- glm(f, data = analysis_df, family = binomial())
  broom::tidy(m) %>%
    filter(term == "slow_risk_z") %>%
    mutate(symptom = sym)
}) %>%
  mutate(
    OR        = exp(estimate),
    conf.low  = exp(estimate - 1.96*std.error),
    conf.high = exp(estimate + 1.96*std.error)
  ) %>%
  arrange(desc(abs(estimate)))

# Panel A plot uses the log-odds (estimate) with Wald 95% CI
sens_df <- symptom_reg %>%
  transmute(
    symptom,
    logOR      = estimate,
    logOR_low  = estimate - 1.96 * std.error,
    logOR_high = estimate + 1.96 * std.error
  ) %>%
  arrange(logOR) %>%
  mutate(symptom = factor(symptom, levels = symptom))

pA <- ggplot(sens_df, aes(x = symptom, y = logOR)) +
  geom_errorbar(aes(ymin = logOR_low, ymax = logOR_high),
                width = 0.2, color = "steelblue4") +
  geom_point(size = 2.8, color = "steelblue") +
  coord_flip() +
  labs(
    title = "(a) Symptom sensitivity to Slow Risk Load",
    x = NULL,
    y = "Sensitivity (log odds per 1 SD Slow Risk Load)"
  ) +
  theme_paper

# Panel B: strength centrality from the shared J, and Spearman association with SRL.
# Strength = sum of absolute edge weights for each node (weighted degree).
W <- J_shared
strength <- rowSums(abs(W), na.rm = TRUE)
names(strength) <- rownames(W)

# Spearman correlation: robust to nonlinearity and ordinal scaling of PHQ items.
slow_corr_spearman <- vapply(symptoms, function(v) {
  cor(analysis_df[[v]], analysis_df$slow_risk_z,
      method = "spearman", use = "pairwise.complete.obs")
}, numeric(1))

plot_tbl <- tibble(
  symptom   = symptoms,
  strength  = as.numeric(strength[symptoms]),
  corr_slow = as.numeric(slow_corr_spearman)
)

pB <- ggplot(plot_tbl, aes(x = strength, y = corr_slow, label = symptom)) +
  geom_point(size = 2.8, color = "gold2") +
  ggrepel::geom_text_repel(size = 3.6, max.overlaps = Inf) +
  labs(
    title = "(b) Central symptoms covary more strongly with Slow Risk Load",
    x = "Strength (shared symptom network)",
    y = expression("Spearman " * rho * " with SRL (z)")
  ) +
  theme_paper

fig5 <- (pA / pB) + plot_layout(heights = c(1, 1))
# ggsave("figs/fig5_slow_to_fast.pdf", fig5, width = 8, height = 8, units = "in")





## ───────────────────────────────────────────────────────────────────────────────
# 13) TABLE 1: Sample characteristics + slow-layer indices + PHQ outcomes
# Uses analysis_df exactly as constructed in your pipeline.
# Output: a data.frame you can feed into knitr::kable(..., format="latex")
# ───────────────────────────────────────────────────────────────────────────────

# ---- formatting helpers ----
fmt_mean_sd <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), sd(x))
}

fmt_median_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  q <- quantile(x, probs = c(.25, .5, .75), names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), q[2], q[1], q[3])
}

fmt_n_pct <- function(x, level = NULL) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(NA_character_)
  if (is.null(level)) stop("Provide `level=` for fmt_n_pct().")
  k <- sum(as.character(x) == level)
  sprintf("%d (%.1f%%)", k, 100 * k / n)
}

pct_missing <- function(x) {
  sprintf("%.1f%%", 100 * mean(is.na(x)))
}

# ---- binary endorsement consistent with Fig 5 panel (a): x > 0 ----
endorse01 <- function(x) as.integer(x > 0)

# ---- Ensure z-scored versions exist for domains + SRL (for correct “(z)” labels) ----
analysis_df <- analysis_df %>%
  mutate(
    slow_risk_z = as.numeric(scale(slow_risk_load)),
    SES_z       = as.numeric(scale(SES_risk)),
    stress_z    = as.numeric(scale(structural_stress)),
    health_z    = as.numeric(scale(health_risk))
  )

# ---- build Table 1 ----
tab1 <- bind_rows(
  tibble(
    Variable = "N (analytic sample)",
    Summary  = as.character(nrow(analysis_df)),
    Missing  = "—"
  ),
  tibble(
    Variable = "Age (years), mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$age, digits = 1),
    Missing  = pct_missing(analysis_df$age)
  ),
  tibble(
    Variable = "Gender: Female, n (%)",
    Summary  = fmt_n_pct(analysis_df$gender_grp, level = "Female"),
    Missing  = pct_missing(analysis_df$gender_grp)
  ),
  tibble(
    Variable = "Ethnicity: Non-Dutch, n (%)",
    Summary  = fmt_n_pct(analysis_df$dutch_grp, level = "Non-Dutch"),
    Missing  = pct_missing(analysis_df$dutch_grp)
  ),
  tibble(
    Variable = "Slow Risk Load (z), mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$slow_risk_z, digits = 2),
    Missing  = pct_missing(analysis_df$slow_risk_z)
  ),
  tibble(
    Variable = "SES risk domain (z), mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$SES_z, digits = 2),
    Missing  = pct_missing(analysis_df$SES_risk)  # missingness identical to SES_z
  ),
  tibble(
    Variable = "Structural stress domain (z), mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$stress_z, digits = 2),
    Missing  = pct_missing(analysis_df$structural_stress)
  ),
  tibble(
    Variable = "Health risk domain (z), mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$health_z, digits = 2),
    Missing  = pct_missing(analysis_df$health_risk)
  ),
  tibble(
    Variable = "PHQ-9 sum score, mean (SD)",
    Summary  = fmt_mean_sd(analysis_df$PHQsum, digits = 2),
    Missing  = pct_missing(analysis_df$PHQsum)
  )
)


# Add symptom endorsement rates (any endorsement > 0), consistent with Fig 5 panel (a)
symptom_rows <- lapply(symptom_vars, function(v) {
  x <- analysis_df[[v]]
  x01 <- endorse01(x)
  n <- sum(!is.na(x01))
  k <- sum(x01 == 1, na.rm = TRUE)
  tibble(
    Variable = paste0("PHQ-9 item: ", v, " endorsed (>0), n (%)"),
    Summary  = sprintf("%d (%.1f%%)", k, 100 * k / n),
    Missing  = pct_missing(x)
  )
}) %>% bind_rows()

tab1 <- bind_rows(tab1, symptom_rows)
tab1



# ───────────────────────────────────────────────────────────────────────────────
# 14) FIG 6: Overall distributions of slow-layer indices (densities) + PHQ sum
# NOTE: use the z-scored domain scores (SES_z / stress_z / health_z) so the x-axis
# label "Standardized score (z)" is accurate for every panel.
# ───────────────────────────────────────────────────────────────────────────────

slow_dist_long <- analysis_df %>%
  select(slow_risk_z, SES_z, stress_z, health_z) %>%
  pivot_longer(everything(), names_to = "measure", values_to = "value") %>%
  mutate(
    measure = factor(
      measure,
      levels = c("slow_risk_z","SES_z","stress_z","health_z"),
      labels = c("Slow Risk Load","SES risk","Structural stress","Health risk")
    )
  )

p_slow_histdens <- ggplot(slow_dist_long, aes(x = value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    boundary = 0,
    closed = "left",
    fill = "grey80",
    color = "grey60",
    linewidth = 0.1,
    na.rm = TRUE
  ) +
  geom_density(linewidth = 0.2, na.rm = TRUE) +
  facet_wrap(~ measure, ncol = 2, scales = "free_y") +
  labs(x = "Standardized score (z)", y = "Density") +
  theme_paper

p_phq_hist <- ggplot(analysis_df, aes(x = PHQsum)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, boundary = -0.5, closed = "left",
                 fill = "grey75", color = "grey40", linewidth = 0.1, na.rm = TRUE) +
  geom_density(linewidth = 0.2, na.rm = TRUE, adjust = 1.1) +
  scale_x_continuous(breaks = seq(0, 27, by = 3)) +
  labs(x = "PHQ-9 sum score", y = "Density") +
  theme_paper

supp_dist2 <- p_slow_histdens / p_phq_hist +
  plot_layout(heights = c(2, 1), axis_titles = "collect")

# ggsave("figs/fig6_distribution_hist_density.pdf", supp_dist2, width = 7, height = 7.5)

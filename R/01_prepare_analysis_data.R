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
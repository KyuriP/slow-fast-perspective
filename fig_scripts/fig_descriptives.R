# ───────────────────────────────────────────────────────────────────────────────
# FIGURE: Descriptives
# Slow Risk Load distributions + PHQ-9 means, with secondary comparisons veiled
# ───────────────────────────────────────────────────────────────────────────────

veil_alpha <- 0.6

# Simple x-axis labels
xlab_group <- function(x) dplyr::recode(
  x,
  "Low Slow Risk"  = "Low Risk",
  "High Slow Risk" = "High Risk",
  .default = x
)

# ───────────────────────────────────────────────────────────────────────────────
# Plotting data
# ───────────────────────────────────────────────────────────────────────────────

# Long table for Slow Risk Load faceting
slow_facet_df <- analysis_df %>%
  select(slow_risk_z, slow_group, dutch_grp, age_grp) %>%
  pivot_longer(
    c(slow_group, dutch_grp, age_grp),
    names_to = "group_type",
    values_to = "group"
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    group_type = factor(
      group_type,
      levels = c("slow_group", "dutch_grp", "age_grp"),
      labels = c("Slow-Risk Group", "Ethnicity", "Age Group")
    ),
    group = factor(as.character(group), levels = names(group_colors))
  )

# Long table for PHQ-9 means
phq_long <- analysis_df %>%
  select(PHQsum, slow_group, dutch_grp, age_grp) %>%
  pivot_longer(
    c(slow_group, dutch_grp, age_grp),
    names_to = "group_type",
    values_to = "group"
  ) %>%
  filter(!is.na(PHQsum), !is.na(group)) %>%
  mutate(
    group_type = factor(
      group_type,
      levels = c("slow_group", "dutch_grp", "age_grp"),
      labels = c("Slow-Risk Group", "Ethnicity", "Age Group")
    ),
    group = factor(as.character(group), levels = names(group_colors))
  )

# 99% CI
ci_level <- 0.99
z_crit <- qnorm((1 + ci_level) / 2)

phq_sum <- phq_long %>%
  group_by(group_type, group) %>%
  summarise(
    n  = dplyr::n(),
    m  = mean(PHQsum, na.rm = TRUE),
    se = sd(PHQsum, na.rm = TRUE) / sqrt(n),
    lo = m - z_crit * se,
    hi = m + z_crit * se,
    .groups = "drop"
  )

# Facets to visually de-emphasize
# xmin/xmax/ymin/ymax are explicit to avoid geom_rect length warnings.
veil_df <- data.frame(
  group_type = factor(
    c("Ethnicity", "Age Group"),
    levels = c("Slow-Risk Group", "Ethnicity", "Age Group")
  ),
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)

# ───────────────────────────────────────────────────────────────────────────────
# Top panel: PHQ-9 mean + 99% CI
# ───────────────────────────────────────────────────────────────────────────────

p_phq_top <- ggplot(phq_sum, aes(x = group, y = m, color = group)) +
  geom_point(size = 2.8, shape = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.5) +
  geom_rect(
    data = veil_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    alpha = veil_alpha,
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = paste0(
      "Average depressive\nsymptom score (",
      round(ci_level * 100),
      "% CI)"
    )
  ) +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_color_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 2)
  )

# ───────────────────────────────────────────────────────────────────────────────
# Bottom panel: Slow Risk Load distributions
# ───────────────────────────────────────────────────────────────────────────────

p_srl_bottom <- ggplot(
  slow_facet_df,
  aes(x = group, y = slow_risk_z, fill = group)
) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.20) +
  geom_rect(
    data = veil_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    alpha = veil_alpha,
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "Slow Risk Load (standardized)"
  ) +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 2)
  )

# ───────────────────────────────────────────────────────────────────────────────
# Combine
# ───────────────────────────────────────────────────────────────────────────────

fig_descriptives <- p_phq_top / p_srl_bottom +
  plot_layout(heights = c(1, 1))

# ggsave("figs/fig2_refined.pdf", fig_descriptives, width = 7, height = 9, units = "in")

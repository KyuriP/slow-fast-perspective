# ==============================================================================
# Supplementary Figure S1: secondary ethnicity and age descriptives
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

ensure_analysis_data()

group_long <- analysis_df %>%
  select(
    PHQsum,
    slow_risk_z,
    dutch_grp,
    age_grp
  ) %>%
  pivot_longer(
    cols = c(dutch_grp, age_grp),
    names_to = "comparison",
    values_to = "group"
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("dutch_grp", "age_grp"),
      labels = c("Ethnicity", "Age group")
    ),
    group = factor(
      as.character(group),
      levels = names(GROUP_COLORS)
    )
  )

phq_summary <- group_long %>%
  filter(!is.na(PHQsum)) %>%
  group_by(comparison, group) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(PHQsum),
    se = stats::sd(PHQsum) / sqrt(n),
    lower = mean - stats::qnorm(0.975) * se,
    upper = mean + stats::qnorm(0.975) * se,
    .groups = "drop"
  )

p_top <- ggplot(
  phq_summary,
  aes(x = group, y = mean, color = group)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(size = 3.1) +
  facet_wrap(~ comparison, scales = "free_x") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "A  Depressive symptom level",
    x = NULL,
    y = "Mean PHQ-9 total (95% CI)"
  ) +
  theme_revised +
  theme(axis.ticks.x = element_blank())

p_bottom <- ggplot(
  group_long,
  aes(x = group, y = slow_risk_z, fill = group)
) +
  geom_boxplot(
    alpha = 0.82,
    outlier.alpha = 0.10,
    linewidth = 0.55
  ) +
  facet_wrap(~ comparison, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "B  Slow Risk Load",
    x = NULL,
    y = "Slow Risk Load (standardized)"
  ) +
  theme_revised +
  theme(axis.ticks.x = element_blank())

fig_S1 <- p_top / p_bottom +
  plot_layout(heights = c(0.85, 1.15))

save_revised_figure(
  plot = fig_S1,
  filename = "supp_figure_S1_secondary_descriptives",
  width = 7.5,
  height = 7.0
)

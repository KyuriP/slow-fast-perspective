# ==============================================================================
# Main Figure 2: Slow Risk Load and observed depressive symptom expression
#
# Panel A: continuous Slow Risk Load distribution and median split
# Panel B: mean PHQ-9 total by Low/High Slow Risk
# Panel C: observed PHQ-9 symptom prevalence by Low/High Slow Risk
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

ensure_analysis_data()

# ------------------------------------------------------------------------------
# Panel A: Slow Risk distribution
# ------------------------------------------------------------------------------

median_slow_z <- stats::median(analysis_df$slow_risk_z, na.rm = TRUE)

slow_distribution <- analysis_df %>%
  filter(!is.na(slow_risk_z), !is.na(slow_group)) %>%
  mutate(
    slow_group = factor(
      slow_group,
      levels = c("Low Slow Risk", "High Slow Risk")
    )
  )

x_positions <- slow_distribution %>%
  group_by(slow_group) %>%
  summarise(
    x = stats::median(slow_risk_z, na.rm = TRUE),
    .groups = "drop"
  )

pA <- ggplot(
  slow_distribution,
  aes(x = slow_risk_z, fill = slow_group)
) +
  geom_histogram(
    bins = 55,
    alpha = 0.80,
    position = "identity",
    color = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = median_slow_z,
    linetype = 2,
    linewidth = 0.75,
    color = DARK_GREY
  ) +
  geom_text(
    data = x_positions,
    aes(
      x = x,
      y = Inf,
      label = ifelse(
        slow_group == "Low Slow Risk",
        "Low Slow Risk",
        "High Slow Risk"
      ),
      color = slow_group
    ),
    inherit.aes = FALSE,
    vjust = 1.45,
    fontface = "bold",
    family = FONT_FAMILY,
    size = 3.8
  ) +
  annotate(
    "text",
    x = median_slow_z,
    y = Inf,
    label = "Median split",
    vjust = 3.05,
    hjust = -0.08,
    family = FONT_FAMILY,
    size = 3.3,
    color = DARK_GREY
  ) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "A  Slow Risk Load distribution",
    x = "Slow Risk Load (standardized)",
    y = "Participants"
  ) +
  theme_revised +
  theme(
    panel.grid.major.x = element_blank()
  )

# ------------------------------------------------------------------------------
# Panel B: PHQ-9 total
# ------------------------------------------------------------------------------

phq_summary <- analysis_df %>%
  filter(!is.na(PHQsum), !is.na(slow_group)) %>%
  group_by(slow_group) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(PHQsum),
    sd = stats::sd(PHQsum),
    se = sd / sqrt(n),
    lower = mean - stats::qnorm(0.975) * se,
    upper = mean + stats::qnorm(0.975) * se,
    .groups = "drop"
  ) %>%
  mutate(
    slow_group = factor(
      slow_group,
      levels = c("Low Slow Risk", "High Slow Risk")
    ),
    group_label = paste0(
      ifelse(
        slow_group == "Low Slow Risk",
        "Low Slow Risk",
        "High Slow Risk"
      ),
      "\n",
      "n = ",
      scales::comma(n)
    )
  )

pB <- ggplot(
  phq_summary,
  aes(x = slow_group, y = mean, color = slow_group)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_point(size = 3.6) +
  geom_text(
    aes(label = sprintf("%.2f", mean)),
    nudge_y = 0.42,
    family = FONT_FAMILY,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = setNames(phq_summary$group_label, phq_summary$slow_group)) +
  scale_y_continuous(
    limits = c(0, max(phq_summary$upper) + 1.0),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "B  Depressive symptom level",
    x = NULL,
    y = "Mean PHQ-9 total (95% CI)"
  ) +
  theme_revised +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ------------------------------------------------------------------------------
# Panel C: symptom prevalence
# ------------------------------------------------------------------------------

prevalence_table <- purrr::map_dfr(SYMPTOMS, function(symptom_name) {
  analysis_df_bin %>%
    filter(!is.na(slow_group), !is.na(.data[[symptom_name]])) %>%
    group_by(slow_group) %>%
    summarise(
      symptom = symptom_name,
      active = sum(.data[[symptom_name]] == 1),
      n = dplyr::n(),
      prevalence = active / n,
      .groups = "drop"
    ) %>%
    mutate(
      lower = purrr::map2_dbl(
        active,
        n,
        ~ unname(binomial_ci(.x, .y)[["lower"]])
      ),
      upper = purrr::map2_dbl(
        active,
        n,
        ~ unname(binomial_ci(.x, .y)[["upper"]])
      )
    )
}) %>%
  mutate(
    symptom_label = unname(SYMPTOM_LABELS[symptom]),
    symptom_label = factor(
      symptom_label,
      levels = rev(unname(SYMPTOM_LABELS[SYMPTOMS]))
    ),
    slow_group = factor(
      slow_group,
      levels = c("Low Slow Risk", "High Slow Risk")
    )
  )

write.csv(
  prevalence_table,
  "results/figure_data/figure_02_symptom_prevalence.csv",
  row.names = FALSE
)

pC <- ggplot(
  prevalence_table,
  aes(
    x = symptom_label,
    y = prevalence,
    color = slow_group,
    shape = slow_group
  )
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.18,
    linewidth = 0.55,
    position = position_dodge(width = 0.48)
  ) +
  geom_point(
    size = 2.55,
    position = position_dodge(width = 0.48)
  ) +
  coord_flip() +
  scale_color_manual(values = GROUP_COLORS) +
  scale_shape_manual(values = c(
    "Low Slow Risk" = 16,
    "High Slow Risk" = 17
  )) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, max(prevalence_table$upper) + 0.035)
  ) +
  labs(
    title = "C  Observed symptom prevalence",
    x = NULL,
    y = "Participants reporting the symptom",
    color = NULL,
    shape = NULL
  ) +
  theme_revised +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.major.y = element_blank()
  )

top_row <- pA + pB +
  plot_layout(widths = c(1.45, 0.85))

fig_02 <- top_row / pC +
  plot_layout(heights = c(0.82, 1.28))

save_revised_figure(
  plot = fig_02,
  filename = "figure_02_slowrisk_descriptives",
  width = 8.2,
  height = 8.5
)

# ==============================================================================
# fix_figure_02_only.R
#
# Final repair of Main Figure 2.
#
# Run from the repository root:
#   source("fix_figure_02_only.R")
# ==============================================================================

source("fig_scripts/revised/00_revised_figure_setup.R")
ensure_analysis_data()

median_slow_z <- stats::median(
  analysis_df$slow_risk_z,
  na.rm = TRUE
)

slow_distribution <- analysis_df %>%
  dplyr::filter(
    !is.na(slow_risk_z),
    !is.na(slow_group)
  ) %>%
  dplyr::mutate(
    slow_group = factor(
      slow_group,
      levels = c(
        "Low Slow Risk",
        "High Slow Risk"
      )
    )
  )

# ------------------------------------------------------------------------------
# A. Slow Risk Load distribution
# ------------------------------------------------------------------------------

pA <- ggplot(
  slow_distribution,
  aes(
    x = slow_risk_z,
    fill = slow_group
  )
) +
  geom_histogram(
    bins = 55,
    alpha = 0.82,
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
  annotate(
    "label",
    x = median_slow_z + 0.12,
    y = Inf,
    label = "Median",
    hjust = 0,
    vjust = 1.35,
    family = FONT_FAMILY,
    size = 3.15,
    label.size = 0,
    fill = scales::alpha("white", 0.85),
    color = DARK_GREY
  ) +
  scale_fill_manual(
    values = GROUP_COLORS,
    guide = "none"
  ) +
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
# B. PHQ-9 total
# ------------------------------------------------------------------------------

phq_summary <- analysis_df %>%
  dplyr::filter(
    !is.na(PHQsum),
    !is.na(slow_group)
  ) %>%
  dplyr::group_by(slow_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(PHQsum),
    sd = stats::sd(PHQsum),
    se = sd / sqrt(n),
    lower = mean - stats::qnorm(0.975) * se,
    upper = mean + stats::qnorm(0.975) * se,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    slow_group = factor(
      slow_group,
      levels = c(
        "Low Slow Risk",
        "High Slow Risk"
      )
    )
  )

pB <- ggplot(
  phq_summary,
  aes(
    y = slow_group,
    x = mean,
    color = slow_group
  )
) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point(size = 3.7) +
  geom_text(
    aes(label = sprintf("%.2f", mean)),
    nudge_x = 0.32,
    hjust = 0,
    family = FONT_FAMILY,
    size = 3.6,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 8.2),
    breaks = c(0, 2, 4, 6, 8),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "B  PHQ-9 total score",
    x = "Mean PHQ-9 total (95% CI)",
    y = NULL
  ) +
  theme_revised +
  theme(
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

# ------------------------------------------------------------------------------
# C. Observed symptom prevalence
# ------------------------------------------------------------------------------

prevalence_table <- purrr::map_dfr(
  SYMPTOMS,
  function(symptom_name) {
    analysis_df_bin %>%
      dplyr::filter(
        !is.na(slow_group),
        !is.na(.data[[symptom_name]])
      ) %>%
      dplyr::group_by(slow_group) %>%
      dplyr::summarise(
        symptom = symptom_name,
        active = sum(.data[[symptom_name]] == 1),
        n = dplyr::n(),
        prevalence = active / n,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
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
  }
) %>%
  dplyr::mutate(
    symptom_label = unname(
      SYMPTOM_LABELS[symptom]
    ),
    symptom_label = factor(
      symptom_label,
      levels = rev(
        unname(SYMPTOM_LABELS[SYMPTOMS])
      )
    ),
    slow_group = factor(
      slow_group,
      levels = c(
        "Low Slow Risk",
        "High Slow Risk"
      )
    )
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
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.18,
    linewidth = 0.55,
    position = position_dodge(width = 0.48)
  ) +
  geom_point(
    size = 2.55,
    position = position_dodge(width = 0.48)
  ) +
  coord_flip() +
  scale_color_manual(
    values = GROUP_COLORS
  ) +
  scale_shape_manual(
    values = c(
      "Low Slow Risk" = 16,
      "High Slow Risk" = 17
    )
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(
      0,
      max(prevalence_table$upper) + 0.035
    )
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
  patchwork::plot_layout(
    widths = c(1.18, 0.82)
  )

fig_02 <- top_row / pC +
  patchwork::plot_layout(
    heights = c(0.82, 1.28)
  )

save_revised_figure(
  plot = fig_02,
  filename = "figure_02_slowrisk_descriptives",
  width = 8.6,
  height = 8.5
)

message(
  "\nFinal Figure 2 saved to figs/revised/.\n"
)


# ==============================================================================
# Main Figure 3: activation differences from the fully group-specific model
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

input_file <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "exact_multigroup_ising",
  "activation_differences_free_J_free_h_with_95CI.csv"
)

require_result_file(input_file)

activation <- read.csv(
  input_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    symptom_label = ifelse(
      is.na(symptom_label) | symptom_label == "",
      unname(SYMPTOM_LABELS[symptom]),
      symptom_label
    ),
    symptom_label = ifelse(
      symptom == "con",
      "Concentration problems",
      symptom_label
    ),
    direction = ifelse(
      difference_high_minus_low >= 0,
      "Higher in High Slow Risk",
      "Higher in Low Slow Risk"
    )
  ) %>%
  arrange(difference_high_minus_low) %>%
  mutate(
    symptom_label = factor(
      symptom_label,
      levels = symptom_label
    ),
    direction = factor(
      direction,
      levels = c(
        "Higher in Low Slow Risk",
        "Higher in High Slow Risk"
      )
    )
  )

write.csv(
  activation,
  "results/figure_data/figure_03_activation_differences.csv",
  row.names = FALSE
)

x_min <- min(activation$ci_95_lower) - 0.04
x_max <- max(activation$ci_95_upper) + 0.04

fig_03 <- ggplot(
  activation,
  aes(
    x = difference_high_minus_low,
    y = symptom_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.65,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = ci_95_lower,
      xend = ci_95_upper,
      yend = symptom_label
    ),
    linewidth = 0.9,
    color = "grey35"
  ) +
  geom_point(
    aes(fill = direction),
    shape = 21,
    size = 3.5,
    stroke = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Higher in Low Slow Risk" = LOW_COLOR,
      "Higher in High Slow Risk" = HIGH_COLOR
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq(-0.4, 0.8, by = 0.2)
  ) +
  coord_cartesian(
    xlim = c(x_min, x_max)
  ) +
  labs(
    title = "Symptom activation differences",
    subtitle = paste(
      "\u2190 Higher in Low Slow Risk",
      "     Higher in High Slow Risk \u2192"
    ),
    x = "Activation difference (High minus Low Slow Risk)",
    y = NULL
  ) +
  theme_revised +
  theme(
    plot.subtitle = element_text(
      hjust = 0.5,
      color = DARK_GREY,
      margin = margin(b = 6)
    ),
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.45
    ),
    axis.text.y = element_text(size = 11)
  )

save_revised_figure(
  plot = fig_03,
  filename = "figure_03_activation_differences",
  width = 7.5,
  height = 5.8
)

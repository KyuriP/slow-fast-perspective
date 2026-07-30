# ==============================================================================
# Main Figure 4: symptom-specific associations with continuous Slow Risk Load
#
# Replaces the old two-panel symptom-map figure. Only the coefficient plot is
# retained in the main manuscript.
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

ensure_analysis_data()

analysis_bin <- analysis_df %>%
  mutate(
    across(
      all_of(SYMPTOMS),
      ~ as.integer(.x > 0),
      .names = "{.col}_bin"
    )
  )

fit_one_symptom <- function(symptom_name, adjust_for_other_symptoms = FALSE) {
  response <- paste0(symptom_name, "_bin")

  predictors <- c(
    "slow_risk_z",
    "age",
    "gender_grp",
    "dutch_grp"
  )

  if (isTRUE(adjust_for_other_symptoms)) {
    predictors <- c(
      predictors,
      paste0(setdiff(SYMPTOMS, symptom_name), "_bin")
    )
  }

  formula <- reformulate(
    termlabels = predictors,
    response = response
  )

  model <- glm(
    formula,
    data = analysis_bin,
    family = binomial()
  )

  broom::tidy(model) %>%
    filter(term == "slow_risk_z") %>%
    transmute(
      symptom = symptom_name,
      model = ifelse(
        adjust_for_other_symptoms,
        "Symptom-adjusted",
        "Covariate-adjusted"
      ),
      estimate = estimate,
      std_error = std.error,
      lower = estimate - stats::qnorm(0.975) * std.error,
      upper = estimate + stats::qnorm(0.975) * std.error,
      p_value = p.value
    )
}

coef_table <- bind_rows(
  purrr::map_dfr(
    SYMPTOMS,
    fit_one_symptom,
    adjust_for_other_symptoms = FALSE
  ),
  purrr::map_dfr(
    SYMPTOMS,
    fit_one_symptom,
    adjust_for_other_symptoms = TRUE
  )
)

symptom_order <- coef_table %>%
  filter(model == "Symptom-adjusted") %>%
  arrange(estimate) %>%
  pull(symptom)

coef_table <- coef_table %>%
  mutate(
    symptom_label = unname(SYMPTOM_LABELS[symptom]),
    symptom_label = factor(
      symptom_label,
      levels = unname(SYMPTOM_LABELS[symptom_order])
    ),
    model = factor(
      model,
      levels = c("Covariate-adjusted", "Symptom-adjusted")
    )
  )

write.csv(
  coef_table,
  "results/figure_data/figure_04_slowrisk_coefficients.csv",
  row.names = FALSE
)

fig_04 <- ggplot(
  coef_table,
  aes(
    x = estimate,
    y = symptom_label,
    color = model,
    shape = model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.55,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    width = 0.18,
    linewidth = 0.65,
    position = position_dodge(width = 0.52)
  ) +
  geom_point(
    size = 3.0,
    position = position_dodge(width = 0.52)
  ) +
  scale_color_manual(
    values = c(
      "Covariate-adjusted" = GREY_COLOR,
      "Symptom-adjusted" = BLUE_COLOR
    )
  ) +
  scale_shape_manual(
    values = c(
      "Covariate-adjusted" = 17,
      "Symptom-adjusted" = 16
    )
  ) +
  labs(
    title = "Associations with continuous Slow Risk Load vary across symptoms",
    x = "Association with symptom presence\n(log-odds per 1 SD higher Slow Risk Load)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_revised +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    axis.text.y = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.45)
  )

save_revised_figure(
  plot = fig_04,
  filename = "figure_04_slowrisk_coefficients",
  width = 7.6,
  height = 6.1
)

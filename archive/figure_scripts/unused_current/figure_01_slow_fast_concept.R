# ==============================================================================
# Main Figure 1: revised slow-fast conceptual framework
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

nodes <- tibble::tribble(
  ~name,            ~x,   ~y,
  "Sleep",          7.25, 5.05,
  "Energy",         8.40, 4.45,
  "Mood",           8.15, 3.25,
  "Activity",       6.75, 3.00,
  "Concentration",  6.25, 4.15
)

edges <- tibble::tribble(
  ~from,            ~to,
  "Sleep",          "Energy",
  "Energy",         "Mood",
  "Mood",           "Activity",
  "Activity",       "Concentration",
  "Concentration",  "Sleep",
  "Energy",         "Concentration",
  "Mood",           "Concentration"
) %>%
  left_join(nodes %>% select(from = name, x_from = x, y_from = y), by = "from") %>%
  left_join(nodes %>% select(to = name, x_to = x, y_to = y), by = "to")

slow_items <- tibble::tribble(
  ~label,                  ~x,   ~y,
  "Economic resources",   2.25, 5.25,
  "Work and housing",     2.25, 4.55,
  "Psychosocial stress",  2.25, 3.85,
  "Physical health",      2.25, 3.15,
  "Lifestyle conditions", 2.25, 2.45
)

fig_01 <- ggplot() +
  annotate(
    "rect",
    xmin = 0.45, xmax = 4.15,
    ymin = 1.75, ymax = 6.35,
    fill = "#E7F4EF",
    color = LOW_COLOR,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = 2.30, y = 5.95,
    label = "Slower contextual conditions",
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.7
  ) +
  geom_label(
    data = slow_items,
    aes(x = x, y = y, label = label),
    family = FONT_FAMILY,
    size = 3.45,
    label.size = 0.25,
    label.padding = grid::unit(0.16, "lines"),
    fill = "white",
    color = DARK_GREY
  ) +

  annotate(
    "rect",
    xmin = 5.85, xmax = 9.55,
    ymin = 1.75, ymax = 6.35,
    fill = "#EDF4F9",
    color = BLUE_COLOR,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = 7.70, y = 5.95,
    label = "Fast symptom layer",
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.7
  ) +
  geom_segment(
    data = edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    color = "grey58",
    linewidth = 0.9
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 8.5,
    stroke = 0.7,
    color = DARK_GREY,
    fill = "white"
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = name),
    family = FONT_FAMILY,
    size = 2.85,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 7.70, y = 2.25,
    label = "Symptoms are linked by pairwise interactions",
    family = FONT_FAMILY,
    size = 3.35,
    color = DARK_GREY
  ) +

  # Slow-to-fast modulation.
  annotate(
    "segment",
    x = 4.25, y = 4.55,
    xend = 5.72, yend = 4.55,
    linewidth = 1.0,
    color = HIGH_COLOR,
    arrow = grid::arrow(
      length = grid::unit(0.22, "cm"),
      type = "closed"
    )
  ) +
  annotate(
    "text",
    x = 4.98, y = 5.02,
    label = "may shift\nsymptom activation",
    family = FONT_FAMILY,
    size = 3.35,
    color = HIGH_COLOR,
    lineheight = 0.95
  ) +

  # Fast-to-slow feedback.
  annotate(
    "curve",
    x = 5.72, y = 2.15,
    xend = 4.25, yend = 2.15,
    curvature = -0.35,
    linewidth = 0.8,
    linetype = 2,
    color = DARK_GREY,
    arrow = grid::arrow(
      length = grid::unit(0.18, "cm"),
      type = "closed"
    )
  ) +
  annotate(
    "text",
    x = 4.98, y = 1.55,
    label = "symptoms may feed back over time",
    family = FONT_FAMILY,
    size = 3.10,
    color = DARK_GREY
  ) +

  annotate(
    "rect",
    xmin = 1.05, xmax = 8.95,
    ymin = 0.25, ymax = 1.05,
    fill = "white",
    color = "grey55",
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = 5.00, y = 0.65,
    label = paste(
      "Groups may differ in symptom activation,",
      "symptom interactions, or both."
    ),
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.05
  ) +
  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0, 6.8),
    clip = "off"
  ) +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    plot.margin = margin(8, 8, 8, 8)
  )

save_revised_figure(
  plot = fig_01,
  filename = "figure_01_slow_fast_concept",
  width = 8.3,
  height = 5.6
)

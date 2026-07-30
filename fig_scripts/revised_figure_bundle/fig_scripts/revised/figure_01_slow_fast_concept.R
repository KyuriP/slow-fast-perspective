# ==============================================================================
# Main Figure 1: revised slow-fast conceptual framework
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

# Network nodes in the fast layer.
nodes <- tibble::tribble(
  ~name,         ~x,   ~y,
  "Sleep",       7.00, 5.10,
  "Energy",      8.40, 4.55,
  "Mood",        8.15, 3.05,
  "Activity",    6.55, 2.75,
  "Concentration", 6.10, 4.05
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
  ~label,                      ~x,   ~y,
  "Economic resources",       2.45, 5.70,
  "Work and housing",         2.45, 5.05,
  "Psychosocial stress",      2.45, 4.40,
  "Physical health",          2.45, 3.75,
  "Lifestyle conditions",     2.45, 3.10
)

fig_01 <- ggplot() +
  # Slow-layer box
  annotate(
    "rect",
    xmin = 0.45, xmax = 4.45,
    ymin = 2.40, ymax = 6.55,
    fill = "#E7F4EF",
    color = LOW_COLOR,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = 2.45, y = 6.12,
    label = "Slower contextual conditions",
    family = FONT_FAMILY,
    fontface = "bold",
    size = 5.1
  ) +
  geom_label(
    data = slow_items,
    aes(x = x, y = y, label = label),
    family = FONT_FAMILY,
    size = 3.7,
    label.size = 0.25,
    fill = "white",
    color = DARK_GREY
  ) +

  # Fast-layer box
  annotate(
    "rect",
    xmin = 5.15, xmax = 9.70,
    ymin = 1.45, ymax = 6.55,
    fill = "#EDF4F9",
    color = BLUE_COLOR,
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = 7.43, y = 6.12,
    label = "Faster-changing symptom system",
    family = FONT_FAMILY,
    fontface = "bold",
    size = 5.1
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
    size = 9.0,
    stroke = 0.7,
    color = DARK_GREY,
    fill = "white"
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = name),
    family = FONT_FAMILY,
    size = 3.0,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 7.40, y = 1.82,
    label = "Edges represent symptom interactions  (J\u1d62\u2c7c)",
    family = FONT_FAMILY,
    size = 3.5,
    color = DARK_GREY
  ) +

  # Slow-to-fast arrow
  annotate(
    "segment",
    x = 4.48, y = 4.85,
    xend = 5.05, yend = 4.85,
    linewidth = 1.0,
    color = HIGH_COLOR,
    arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed")
  ) +
  annotate(
    "text",
    x = 4.77, y = 5.35,
    label = "may shift\nsymptom activation  (h\u1d62)",
    family = FONT_FAMILY,
    size = 3.35,
    color = HIGH_COLOR,
    lineheight = 0.95
  ) +

  # Fast-to-slow feedback arrow
  annotate(
    "curve",
    x = 5.20, y = 2.15,
    xend = 4.38, yend = 2.78,
    curvature = -0.35,
    linewidth = 0.8,
    linetype = 2,
    color = DARK_GREY,
    arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")
  ) +
  annotate(
    "text",
    x = 4.95, y = 1.72,
    label = "symptoms may feed back over time",
    family = FONT_FAMILY,
    size = 3.2,
    color = DARK_GREY
  ) +

  # Bottom synthesis statement
  annotate(
    "rect",
    xmin = 1.15, xmax = 9.00,
    ymin = 0.20, ymax = 1.00,
    fill = "white",
    color = "grey55",
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = 5.08, y = 0.60,
    label = "Group differences may occur in symptom activation, symptom interactions, or both.",
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.25
  ) +
  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0, 7),
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
  height = 5.8
)

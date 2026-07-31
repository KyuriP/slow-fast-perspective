# ───────────────────────────────────────────────────────────────────────────────
# FIGURE: Slow Risk Load symptom map
#
#   (a) Symptom-specific links to Slow Risk Load
#   (b) Symptom-adjusted Slow Risk Load associations in the shared symptom network
# ───────────────────────────────────────────────────────────────────────────────

# Required packages
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(broom)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggrepel)
library(grid)
library(patchwork)
library(scales)

# Assumes these objects already exist:
#   analysis_df
#   symptom_vars
#   joint_global       # from the joint shared-network Ising model
#   theme_paper


# ───────────────────────────────────────────────────────────────────────────────
# 0) Setup
# ───────────────────────────────────────────────────────────────────────────────

make_binary <- function(x) as.integer(x > 0)

# Symptom-label dictionary used throughout the figure
symptom_names <- c(
  dep = "depressed mood",
  mot = "psychomotor change",
  glt = "guilt",
  sui = "suicidality",
  app = "appetite change",
  slp = "sleep problems",
  con = "concentration",
  anh = "anhedonia",
  ene = "low energy"
)

# Colors
col_covariate_adjusted <- "grey65"
col_symptom_adjusted   <- "#2C7FB8"

# Font
pal_family <- theme_paper$text$family

if (is.null(pal_family) || identical(pal_family, "")) {
  pal_family <- "Palatino"
}

# Text controls
title_size       <- 19
axis_title_size  <- 17
axis_text_size   <- 16
legend_text_size <- 16

node_label_size <- 5.4
legend_title_sz <- 14
legend_text_sz  <- 13
abbr_title_sz   <- 4.8
abbr_text_sz    <- 4.8
abbr_dot_size   <- 4.8


# ───────────────────────────────────────────────────────────────────────────────
# 1) Prepare binary symptom data
# ───────────────────────────────────────────────────────────────────────────────

analysis_bin <- analysis_df %>%
  mutate(
    across(
      all_of(symptom_vars),
      ~ make_binary(.x),
      .names = "{.col}_bin"
    )
  )

# Check naming consistency with the shared Ising network
if (!identical(sort(symptom_vars), sort(rownames(joint_global$J)))) {
  warning("symptom_vars and rownames(joint_global$J) differ; check naming consistency.")
}


# ───────────────────────────────────────────────────────────────────────────────
# 2) Covariate-adjusted logistic regressions
#    symptom_j ~ Slow Risk Load + age + gender + ethnicity
# ───────────────────────────────────────────────────────────────────────────────

symptom_reg_covariate <- purrr::map_dfr(symptom_vars, function(sym) {
  response <- paste0(sym, "_bin")
  
  f <- reformulate(
    termlabels = c(
      "slow_risk_z",
      "age",
      "gender_grp",
      "dutch_grp"
    ),
    response = response
  )
  
  m <- glm(f, data = analysis_bin, family = binomial())
  
  broom::tidy(m) %>%
    dplyr::filter(term == "slow_risk_z") %>%
    dplyr::mutate(
      symptom = sym,
      model = "Covariate-adjusted"
    )
}) %>%
  mutate(
    logOR      = estimate,
    logOR_low  = estimate - 1.96 * std.error,
    logOR_high = estimate + 1.96 * std.error
  )


# ───────────────────────────────────────────────────────────────────────────────
# 3) Symptom-adjusted logistic regressions
#    symptom_j ~ Slow Risk Load + age + gender + ethnicity + all other symptoms
# ───────────────────────────────────────────────────────────────────────────────

symptom_reg_adjusted <- purrr::map_dfr(symptom_vars, function(sym) {
  response <- paste0(sym, "_bin")
  other_terms <- paste0(setdiff(symptom_vars, sym), "_bin")
  
  f <- reformulate(
    termlabels = c(
      "slow_risk_z",
      "age",
      "gender_grp",
      "dutch_grp",
      other_terms
    ),
    response = response
  )
  
  m <- glm(f, data = analysis_bin, family = binomial())
  
  broom::tidy(m) %>%
    dplyr::filter(term == "slow_risk_z") %>%
    dplyr::mutate(
      symptom = sym,
      model = "Symptom-adjusted"
    )
}) %>%
  mutate(
    logOR      = estimate,
    logOR_low  = estimate - 1.96 * std.error,
    logOR_high = estimate + 1.96 * std.error
  )


# ───────────────────────────────────────────────────────────────────────────────
# 4) Combine model results for panel (a)
# ───────────────────────────────────────────────────────────────────────────────

compare_srl <- bind_rows(
  symptom_reg_covariate %>%
    select(symptom, model, logOR, logOR_low, logOR_high),
  symptom_reg_adjusted %>%
    select(symptom, model, logOR, logOR_low, logOR_high)
)

# Order symptoms by symptom-adjusted association, strongest at top
symptom_order <- symptom_reg_adjusted %>%
  arrange(desc(logOR)) %>%
  pull(symptom)

compare_srl <- compare_srl %>%
  mutate(
    model = factor(
      model,
      levels = c("Covariate-adjusted", "Symptom-adjusted")
    ),
    symptom = factor(
      symptom,
      levels = rev(symptom_order)
    )
  )


# ───────────────────────────────────────────────────────────────────────────────
# 5) Panel A: symptom-specific Slow Risk Load associations
# ───────────────────────────────────────────────────────────────────────────────

pA <- ggplot(
  compare_srl,
  aes(x = symptom, y = logOR, shape = model, color = model)
) +
  geom_errorbar(
    aes(ymin = logOR_low, ymax = logOR_high),
    width = 0.18,
    position = position_dodge(width = 0.55),
    linewidth = 0.65
  ) +
  geom_point(
    size = 2.9,
    position = position_dodge(width = 0.55)
  ) +
  coord_flip() +
  scale_x_discrete(labels = symptom_names) +
  scale_color_manual(
    values = c(
      "Covariate-adjusted" = col_covariate_adjusted,
      "Symptom-adjusted"   = col_symptom_adjusted
    )
  ) +
  scale_shape_manual(
    values = c(
      "Covariate-adjusted" = 17,
      "Symptom-adjusted"   = 16
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3.4,
        linewidth = 1.0
      ),
      byrow = TRUE
    )
  ) +
  labs(
    title = "(a) Symptom-specific links to Slow Risk Load",
    x = NULL,
    y = "Association between Slow Risk Load and symptom activation\n(log-odds change per 1 SD Slow Risk Load)",
    color = NULL,
    shape = NULL
  ) +
  theme_paper +
  theme(
    text = element_text(size = 18, family = pal_family),
    plot.title = element_text(
      size = title_size,
      face = "bold",
      hjust = 0,
      family = pal_family
    ),
    axis.title.x = element_text(size = axis_title_size, family = pal_family),
    axis.text = element_text(size = axis_text_size, family = pal_family),
    legend.text = element_text(size = legend_text_size, family = pal_family),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm"),
    legend.spacing.x = unit(0.55, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )


# ───────────────────────────────────────────────────────────────────────────────
# 6) Panel B: shared network with symptom-adjusted associations
# ───────────────────────────────────────────────────────────────────────────────

# Shared network from jointly estimated Ising model
W_plot <- as.matrix(joint_global$J)
diag(W_plot) <- 0

edge_plot_tbl <- which(upper.tri(W_plot), arr.ind = TRUE) %>%
  as.data.frame() %>%
  transmute(
    from = rownames(W_plot)[row],
    to = colnames(W_plot)[col],
    abs_weight = abs(W_plot[cbind(row, col)])
  ) %>%
  filter(abs_weight > 0)

node_plot_tbl <- symptom_reg_adjusted %>%
  transmute(
    name = symptom,
    symptom_adjusted_logOR = logOR,
    full_name = recode(symptom, !!!symptom_names)
  )

g <- igraph::graph_from_data_frame(
  d = edge_plot_tbl,
  vertices = node_plot_tbl,
  directed = FALSE
)

igraph::E(g)$weight <- edge_plot_tbl$abs_weight

cond_min <- min(node_plot_tbl$symptom_adjusted_logOR, na.rm = TRUE)
cond_max <- max(node_plot_tbl$symptom_adjusted_logOR, na.rm = TRUE)

blue_fun <- scales::col_numeric(
  palette = c("#EAF3FB", "#A8C8E8", "#5B95C8", "#08519C"),
  domain = c(cond_min, cond_max)
)

set.seed(123)


# ───────────────────────────────────────────────────────────────────────────────
# 6.1) Main network plot
# ───────────────────────────────────────────────────────────────────────────────

p_network_main <- ggraph(g, layout = "fr", weights = igraph::E(g)$weight) +
  geom_edge_link(
    aes(edge_width = abs_weight),
    color = "grey80",
    alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = symptom_adjusted_logOR),
    shape = 21,
    size = 7.2,
    color = "black",
    stroke = 0.5
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = node_label_size,
    family = pal_family
  ) +
  scale_edge_width(range = c(0.35, 2.6)) +
  scale_fill_gradientn(
    colours = c("#EAF3FB", "#A8C8E8", "#5B95C8", "#08519C"),
    values = scales::rescale(c(cond_min, 0.27, 0.34, cond_max)),
    limits = c(cond_min, cond_max),
    breaks = c(0.25, 0.30, 0.35, 0.40),
    name = "Symptom-adjusted\nlink to Slow Risk Load"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      barheight = unit(8.8, "cm"),
      barwidth = unit(0.62, "cm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  labs(title = NULL) +
  theme_paper +
  theme(
    text = element_text(family = pal_family, size = 18),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = legend_title_sz, family = pal_family),
    legend.text = element_text(size = legend_text_sz, family = pal_family),
    plot.margin = margin(0, 0, 0, 0)
  )


# ───────────────────────────────────────────────────────────────────────────────
# 6.2) Abbreviation key with colored dots
# ───────────────────────────────────────────────────────────────────────────────

abbr_df <- node_plot_tbl %>%
  mutate(
    label = paste0(name, " = ", full_name),
    fill_col = blue_fun(symptom_adjusted_logOR)
  ) %>%
  arrange(desc(symptom_adjusted_logOR)) %>%
  mutate(y = rev(seq_len(n())))

p_abbrev_key <- ggplot(abbr_df, aes(x = 0, y = y)) +
  geom_point(
    aes(fill = symptom_adjusted_logOR),
    shape = 21,
    size = abbr_dot_size,
    color = "black",
    stroke = 0.4
  ) +
  geom_text(
    aes(x = 0.18, label = label),
    hjust = 0,
    family = pal_family,
    size = abbr_text_sz
  ) +
  scale_fill_gradientn(
    colours = c("#EAF3FB", "#A8C8E8", "#5B95C8", "#08519C"),
    values = scales::rescale(c(cond_min, 0.27, 0.34, cond_max)),
    limits = c(cond_min, cond_max),
    guide = "none"
  ) +
  annotate(
    "text",
    x = 0,
    y = max(abbr_df$y) + 1.2,
    label = "Abbreviations",
    hjust = 0,
    family = pal_family,
    size = abbr_title_sz,
    fontface = "bold"
  ) +
  xlim(0, 2.45) +
  ylim(0.5, max(abbr_df$y) + 1.7) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )


# ───────────────────────────────────────────────────────────────────────────────
# 6.3) Final panel B: network + colorbar + abbreviation key
# ───────────────────────────────────────────────────────────────────────────────

pB_title <- ggplot() +
  labs(
    title = "(b) Symptom-adjusted Slow Risk Load links in the shared network"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = title_size,
      face = "bold",
      hjust = 0,
      family = pal_family
    ),
    plot.margin = margin(0, 0, 2, 0)
  )

pB_body <- p_network_main + p_abbrev_key +
  plot_layout(widths = c(4.6, 3.0))

pB_final <- pB_title / pB_body +
  plot_layout(heights = c(0.10, 1))


# ───────────────────────────────────────────────────────────────────────────────
# Final figure
# ───────────────────────────────────────────────────────────────────────────────

fig_slow_risk_symptom_map <- (pA / pB_final) +
  plot_layout(heights = c(1.05, 1))

fig_slow_risk_symptom_map

# ggsave(
#   filename = "figs/fig5_refined.pdf",
#   plot = fig_slow_risk_symptom_map,
#   width = 12,
#   height = 12,
#   units = "in"
# )
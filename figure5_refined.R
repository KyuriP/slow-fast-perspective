# ───────────────────────────────────────────────────────────────────────────────
# FIGURE 5 (final): 
#   (a) Marginal and conditional associations with slow-risk burden
#   (b) Shared symptom network with highlighted symptom subset
# ───────────────────────────────────────────────────────────────────────────────

# Required packages:
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(tibble)
# library(broom)
# library(ggplot2)
# library(ggrepel)
# library(patchwork)
# library(igraph)
# library(ggraph)

# Assumes these objects already exist:
#   analysis_df
#   symptom_vars
#   symptom_reg          # earlier marginal models
#   J_shared
#   theme_paper

# symptom_reg should be the earlier marginal regressions:
#   make_binary(symptom) ~ slow_risk_z + age + gender_grp + dutch_grp


# ───────────────────────────────────────────────────────────────────────────────
# 0) Setup
# ───────────────────────────────────────────────────────────────────────────────


make_binary <- function(x) as.integer(x > 0)

# Use one consistent symptom-label dictionary throughout the figure
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

# Order for panel (a), strongest conditional association at top
symptom_code_order <- c("dep", "mot", "glt", "sui", "app", "slp", "con", "anh", "ene")

# Colors
col_covariate_adjusted <- "grey65"
col_symptom_adjusted   <- "#2C7FB8"

col_marginal    <- "grey65"
col_conditional <- "#2C7FB8"   # blue
col_other_nodes <- "grey80"
col_cluster     <- "#2C7FB8"
col_other_edges <- "grey75"
col_cluster_edg <- "#2C7FB8"


# ───────────────────────────────────────────────────────────────────────────────
# 1) Conditional / nodewise logistic regressions
#    symptom_j ~ SRL + age + gender + ethnicity + all other symptoms
# ───────────────────────────────────────────────────────────────────────────────

analysis_bin <- analysis_df %>%
  mutate(
    across(
      all_of(symptom_vars),
      ~ make_binary(.x),
      .names = "{.col}_bin"
    )
  )

# Basic checks
if (!identical(sort(symptom_vars), sort(rownames(J_shared)))) {
  warning("symptom_vars and rownames(J_shared) differ; check naming consistency.")
}

symptom_reg_conditional <- purrr::map_dfr(symptom_vars, function(sym) {
  
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
      model = "Conditional"
    )
}) %>%
  mutate(
    logOR      = estimate,
    logOR_low  = estimate - 1.96 * std.error,
    logOR_high = estimate + 1.96 * std.error
  )


# ───────────────────────────────────────────────────────────────────────────────
# 2) Combine with marginal estimates
# ───────────────────────────────────────────────────────────────────────────────

symptom_reg_marginal2 <- symptom_reg %>%
  mutate(
    model      = "Marginal",
    logOR      = estimate,
    logOR_low  = estimate - 1.96 * std.error,
    logOR_high = estimate + 1.96 * std.error
  )

compare_srl <- bind_rows(
  symptom_reg_marginal2 %>%
    select(symptom, model, logOR, logOR_low, logOR_high),
  symptom_reg_conditional %>%
    select(symptom, model, logOR, logOR_low, logOR_high)
)

# Order by conditional coefficients, strongest at top
symptom_order <- symptom_reg_conditional %>%
  arrange(desc(logOR)) %>%
  pull(symptom)

compare_srl <- compare_srl %>%
  mutate(
    model = factor(
      model,
      levels = c("Marginal", "Conditional"),
      labels = c("Covariate-adjusted", "Symptom-adjusted")
    ),
    symptom = factor(symptom, levels = rev(symptom_order))
  )


# ───────────────────────────────────────────────────────────────────────────────
# 3) Panel A
# ───────────────────────────────────────────────────────────────────────────────

pA <- ggplot(
  compare_srl,
  aes(x = symptom, y = logOR, shape = model, color = model)
) +
  geom_errorbar(
    aes(ymin = logOR_low, ymax = logOR_high),
    width = 0.18,
    position = position_dodge(width = 0.55),
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.6,
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
    labs(
      title = "(a) Symptom-specific associations with Slow Risk Load",
      x = NULL,
      y = "Association with symptom activation\n(log odds per 1 SD Slow Risk Load)",
    color = NULL,
    shape = NULL
  ) +
  theme_paper +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 19, face = "bold", hjust = 0),
    axis.title.x = element_text(size = 17),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

# ───────────────────────────────────────────────────────────────────────────────
# Figure 5 text-size + label refinements
# Apply this after pA is created, then replace the network panel code below.
# ───────────────────────────────────────────────────────────────────────────────

# ---- Panel (a): update labels and scale text up ----
pA <- pA +
  labs(
    title = "(a) Log odds per 1 SD Slow Risk Load",
    x = NULL,
    y = "Slow-layer associations with symptom activation",
    color = NULL,
    shape = NULL
  ) +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 19, face = "bold", hjust = 0),
    axis.title.x = element_text(size = 17),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )


# ───────────────────────────────────────────────────────────────────────────────
# Refined network panel for Figure 5
# Node fill = conditional Slow Risk Load association from panel (a)
# Edge width = shared-network coupling strength
# ───────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tibble)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggrepel)
library(grid)
library(patchwork)
library(scales)

# Font
pal_family <- theme_paper$text$family
if (is.null(pal_family) || identical(pal_family, "")) {
  pal_family <- "Palatino"
}



# Larger text controls
title_size      <- 19
node_label_size <- 5.4
legend_title_sz <- 14
legend_text_sz  <- 13
abbr_title_sz   <- 4.8
abbr_text_sz    <- 4.8
abbr_dot_size   <- 4.8

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

# Shared network from jointly estimated model
W_plot <- as.matrix(joint_global$J)
diag(W_plot) <- 0

edge_plot_tbl <- which(upper.tri(W_plot), arr.ind = TRUE) %>%
  as.data.frame() %>%
  transmute(
    from = rownames(W_plot)[row],
    to   = colnames(W_plot)[col],
    abs_weight = abs(W_plot[cbind(row, col)])
  ) %>%
  filter(abs_weight > 0)

node_plot_tbl <- symptom_reg_conditional %>%
  transmute(
    name = symptom,
    conditional_logOR = logOR,
    full_name = recode(symptom, !!!symptom_names)
  )

g <- igraph::graph_from_data_frame(
  d = edge_plot_tbl,
  vertices = node_plot_tbl,
  directed = FALSE
)

igraph::E(g)$weight <- edge_plot_tbl$abs_weight

cond_min <- min(node_plot_tbl$conditional_logOR, na.rm = TRUE)
cond_max <- max(node_plot_tbl$conditional_logOR, na.rm = TRUE)

blue_fun <- scales::col_numeric(
  palette = c("#EAF3FB", "#A8C8E8", "#5B95C8", "#08519C"),
  domain = c(cond_min, cond_max)
)

set.seed(123)

# --- Main network plot ---
p_network_main <- ggraph(g, layout = "fr", weights = igraph::E(g)$weight) +
  geom_edge_link(
    aes(edge_width = abs_weight),
    color = "grey80",
    alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = conditional_logOR),
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
    name = "Conditional association\nwith symptom activation"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      barheight = unit(8.8, "cm"),
      barwidth  = unit(0.62, "cm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  labs(
    title = "(b) Conditional Slow Risk Load associations in the shared network"
  ) +
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
    plot.title = element_text(size = title_size, face = "bold", hjust = 0, family = pal_family),
    plot.margin = margin(0, 0, 0, 0)
  )


# --- Abbreviation key with colored dots ---
abbr_df <- node_plot_tbl %>%
  mutate(
    label = paste0(name, " = ", full_name),
    fill_col = blue_fun(conditional_logOR)
  ) %>%
  arrange(desc(conditional_logOR)) %>%
  mutate(y = rev(seq_len(n())))

p_abbrev_key <- ggplot(abbr_df, aes(x = 0, y = y)) +
  geom_point(
    aes(fill = conditional_logOR),
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
  xlim(0, 2.15) +
  ylim(0.5, max(abbr_df$y) + 1.7) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

# --- Final panel (b): network + abbreviation key ---
pB_final <- p_network_main + p_abbrev_key +
  plot_layout(widths = c(4.6, 2.9))

pB_final


# ───────────────────────────────────────────────────────────────────────────────
# Final Figure 5
# ───────────────────────────────────────────────────────────────────────────────

fig5_final <- (pA / pB_final)
fig5_final

# ggsave(
#   filename = "figs/fig5_refined.pdf",
#   plot = fig5_final,
#   width = 12,
#   height = 12,
#   units = "in"
# )

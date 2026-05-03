# ───────────────────────────────────────────────────────────────────────────────
# FIGURE 2 (REVISED v3)
#
# Approach: keep the original colors entirely, then layer a semi-transparent
# white/gray rectangle over the Ethnicity and Age Group facets only.
# This "veils" the secondary comparisons without touching any of the original
# color or data logic.
#
# Tune `veil_alpha` (0 = invisible, 1 = fully opaque white) to taste.
# 0.45 gives a clear but not total washout.
# ───────────────────────────────────────────────────────────────────────────────

veil_alpha <- 0.6   # ← adjust this one number to control how strong the veil is

# Data frame that tells geom_rect which facets to veil
# (one row per facet that should be grayed out)
veil_df <- data.frame(
  group_type = factor(
    c("Ethnicity", "Age Group"),
    levels = c("Slow-Risk Group", "Ethnicity", "Age Group")
  )
)

# ── re-use ALL original objects (group_colors, phq_sum, slow_facet_df, etc.) ──
# Only the two panel-building blocks change.

# ── TOP PANEL: PHQ-9 mean (99% CI) — original code + veil layer ────────────────
p_phq_top <- ggplot(phq_sum, aes(x = group, y = m, color = group)) +
  geom_point(size = 2.8, shape = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.5) +
  # ── veil: white rectangle over secondary facets ──
  geom_rect(
    data         = veil_df,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill         = "white",
    alpha        = veil_alpha,
    inherit.aes  = FALSE   # don't map color/group from the main data
  ) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  
  labs(x = NULL,     y = paste0("Average depressive\nsymptom score (", round(ci_level * 100), "% CI)")) +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_color_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(
    strip.text   = element_text(face = "bold", size = 13),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(b = 2)
  )

# ── BOTTOM PANEL: SRL boxplots — original code + veil layer ────────────────────
p_srl_bottom <- ggplot(slow_facet_df,
                       aes(x = group, y = slow_risk_z, fill = group)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.20) +
  # ── veil: same white rectangle ──
  geom_rect(
    data        = veil_df,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill        = "white",
    alpha       = veil_alpha,
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ group_type, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Slow Risk Load (standardized)") +
  scale_x_discrete(labels = xlab_group, drop = TRUE) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  theme_paper +
  theme(
    strip.text  = element_blank(),
    plot.margin = margin(t = 2)
  )

# ── COMBINE ────────────────────────────────────────────────────────────────────
fig2_v3 <- p_phq_top / p_srl_bottom +
  plot_layout(heights = c(1, 1))

fig2_v3

# ggsave("figs/fig2_refined.pdf", fig2_v3, width = 7, height = 9, units = "in")

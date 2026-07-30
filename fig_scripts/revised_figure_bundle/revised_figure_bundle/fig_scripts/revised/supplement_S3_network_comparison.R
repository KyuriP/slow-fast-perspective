# ==============================================================================
# Supplementary Figure S3: Low and High Slow Risk NCT networks
#
# Both panels use the same node layout and the same edge-width scale.
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

low_file <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_10000",
  "nct_network_low_risk_10000.csv"
)

high_file <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_10000",
  "nct_network_high_risk_10000.csv"
)

summary_file <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "nct_10000",
  "nct_summary_10000.csv"
)

require_result_file(low_file)
require_result_file(high_file)
require_result_file(summary_file)

read_network_matrix <- function(path) {
  x <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE
  )

  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "numeric"

  matrix_x
}

network_low <- read_network_matrix(low_file)
network_high <- read_network_matrix(high_file)
nct_summary <- read.csv(summary_file, stringsAsFactors = FALSE)

if (!identical(dim(network_low), dim(network_high))) {
  stop("The Low and High network matrices do not have the same dimensions.")
}

if (is.null(rownames(network_low))) {
  rownames(network_low) <- colnames(network_low) <- SYMPTOMS
}

if (is.null(rownames(network_high))) {
  rownames(network_high) <- colnames(network_high) <- SYMPTOMS
}

node_labels <- c(
  anh = "Anh",
  dep = "Dep",
  slp = "Slp",
  ene = "Ene",
  app = "App",
  glt = "Glt",
  con = "Con",
  mot = "Mot",
  sui = "Sui"
)[rownames(network_low)]

layout_shared <- qgraph::averageLayout(
  network_low,
  network_high
)

maximum_edge <- max(
  abs(c(network_low, network_high)),
  na.rm = TRUE
)

pdf_path <- file.path(
  "figs",
  "revised",
  "supp_figure_S3_network_comparison.pdf"
)

png_path <- file.path(
  "figs",
  "revised",
  "supp_figure_S3_network_comparison.png"
)

draw_networks <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(2.0, 1.0, 4.0, 1.0),
    family = FONT_FAMILY
  )

  qgraph::qgraph(
    network_low,
    layout = layout_shared,
    labels = node_labels,
    label.cex = 1.15,
    color = "white",
    border.color = LOW_COLOR,
    border.width = 2,
    vsize = 7.6,
    esize = 8,
    maximum = maximum_edge,
    minimum = 0,
    cut = 0,
    details = FALSE,
    legend = FALSE,
    title = paste0(
      "Low Slow Risk\nGlobal strength = ",
      sprintf("%.2f", nct_summary$global_strength_group_1)
    )
  )

  qgraph::qgraph(
    network_high,
    layout = layout_shared,
    labels = node_labels,
    label.cex = 1.15,
    color = "white",
    border.color = HIGH_COLOR,
    border.width = 2,
    vsize = 7.6,
    esize = 8,
    maximum = maximum_edge,
    minimum = 0,
    cut = 0,
    details = FALSE,
    legend = FALSE,
    title = paste0(
      "High Slow Risk\nGlobal strength = ",
      sprintf("%.2f", nct_summary$global_strength_group_2)
    )
  )
}

grDevices::cairo_pdf(
  filename = pdf_path,
  width = 10.5,
  height = 5.5
)

draw_networks()
grDevices::dev.off()

grDevices::png(
  filename = png_path,
  width = 3360,
  height = 1760,
  res = 320,
  bg = "white"
)

draw_networks()
grDevices::dev.off()

message("Saved: ", pdf_path)
message("Saved: ", png_path)

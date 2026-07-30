# ==============================================================================
# Main Table 2: exact multigroup Ising model comparison
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_revised_figure_setup.R")
}

input_file <- file.path(
  "results",
  "network_invariance_exact",
  "slow_risk",
  "exact_multigroup_ising",
  "model_comparison_exact.csv"
)

require_result_file(input_file)

model_names <- c(
  free_J_free_h = "Group-specific $J$, group-specific $h$",
  equal_J_free_h = "Common $J$, group-specific $h$",
  free_J_equal_h = "Group-specific $J$, common $h$",
  equal_J_equal_h = "Common $J$, common $h$"
)

model_table <- read.csv(
  input_file,
  stringsAsFactors = FALSE
) %>%
  mutate(
    model_label = unname(model_names[model])
  ) %>%
  select(
    model,
    model_label,
    npar,
    logLik,
    AIC,
    delta_AIC,
    BIC,
    delta_BIC
  )

write.csv(
  model_table,
  "tables/revised/table_02_model_comparison.csv",
  row.names = FALSE
)

format_number <- function(x, digits = 2) {
  formatC(
    x,
    format = "f",
    digits = digits,
    big.mark = ","
  )
}

rows <- apply(model_table, 1, function(row) {
  paste0(
    row[["model_label"]],
    " & ",
    row[["npar"]],
    " & ",
    format_number(as.numeric(row[["logLik"]])),
    " & ",
    format_number(as.numeric(row[["AIC"]])),
    " & ",
    format_number(as.numeric(row[["delta_AIC"]])),
    " & ",
    format_number(as.numeric(row[["BIC"]])),
    " & ",
    format_number(as.numeric(row[["delta_BIC"]])),
    " \\\\"
  )
})

latex_table <- c(
  "\\begin{table}[htbp]",
  "\\caption{Exact multigroup Ising model comparison for Low and High Slow Risk groups.}",
  "\\label{tab:ising_model_comparison}",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabular}{@{}lrrrrrr@{}}",
  "\\toprule",
  "\\textbf{Model} & \\textbf{$k$} & \\textbf{LogLik} & \\textbf{AIC} &",
  "\\textbf{$\\Delta$AIC} & \\textbf{BIC} & \\textbf{$\\Delta$BIC} \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\vspace{3pt}",
  "\\footnotesize",
  "\\emph{Note.} $J$ denotes pairwise symptom-interaction parameters and",
  "$h$ denotes symptom-activation parameters. Lower AIC and BIC values",
  "indicate better relative support under the respective criterion.",
  "\\end{table}"
)

writeLines(
  latex_table,
  "tables/revised/table_02_model_comparison.tex"
)

message("Saved: tables/revised/table_02_model_comparison.csv")
message("Saved: tables/revised/table_02_model_comparison.tex")

# Revised figure bundle for `slow-fast-perspective`

This bundle replaces the old figure pipeline that centered the joint shared-\(J\)
model.

## 1. Copy the files into the repository

Copy the folder:

```text
fig_scripts/revised/
```

into the existing repository so the final path is:

```text
slow-fast-perspective/
  fig_scripts/
    revised/
```

The scripts expect the current repository structure:

```text
R/00_setup.R
R/01_prepare_analysis_data.R
R/02_model_helpers.R
R/03_prepare_binary_data.R
data/helius.sav
```

They also expect the finalized analysis outputs:

```text
results/network_invariance_exact/slow_risk/
  group_counts.csv

results/network_invariance_exact/slow_risk/exact_multigroup_ising/
  exact_multigroup_models.rds
  model_comparison_exact.csv
  activation_differences_free_J_free_h_with_95CI.csv

results/network_invariance_exact/slow_risk/nct_10000/
  nct_summary_10000.csv
  nct_network_low_risk_10000.csv
  nct_network_high_risk_10000.csv
```

## 2. Run everything

From the repository root:

```r
source("fig_scripts/revised/make_revised_figures.R")
```

## 3. Main outputs

```text
figs/revised/figure_01_slow_fast_concept.pdf
figs/revised/figure_02_slowrisk_descriptives.pdf
figs/revised/figure_03_activation_differences.pdf
figs/revised/figure_04_slowrisk_coefficients.pdf
tables/revised/table_02_model_comparison.tex
```

PNG versions are also generated.

## 4. Supplementary outputs

```text
figs/revised/supp_figure_S1_secondary_descriptives.pdf
figs/revised/supp_figure_S2_activation_sensitivity.pdf
figs/revised/supp_figure_S3_network_comparison.pdf
```

## 5. Manuscript insertion code

### Figure 1

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]
{figs/revised/figure_01_slow_fast_concept.pdf}
\caption{\textbf{Slow--fast perspective on contextual conditions and depressive
symptoms.} Contextual conditions generally change more slowly than depressive
symptoms and may shift symptom activation. Symptoms are connected through
pairwise interactions and may also feed back into contextual conditions over
time. Groups may differ in symptom activation, symptom interactions, or both.}
\label{fig:slowfast_concept}
\end{figure}
```

### Figure 2

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]
{figs/revised/figure_02_slowrisk_descriptives.pdf}
\caption{\textbf{Slow Risk Load and depressive symptom expression.}
\textbf{(A)} Distribution of standardized Slow Risk Load and the median used
to define Low and High Slow Risk groups.
\textbf{(B)} Mean PHQ-9 total scores with 95\% confidence intervals.
\textbf{(C)} Observed prevalence of each binary PHQ-9 symptom with 95\%
confidence intervals.}
\label{fig:slowrisk_descriptives}
\end{figure}
```

### Table 2

Insert:

```latex
\input{tables/revised/table_02_model_comparison.tex}
```

### Figure 3

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.82\textwidth]
{figs/revised/figure_03_activation_differences.pdf}
\caption{\textbf{Differences in symptom activation between High and Low Slow
Risk groups under the fully group-specific Ising model.} Points show
$\Delta h_i=h_{i,\mathrm{High}}-h_{i,\mathrm{Low}}$ and horizontal lines show
95\% Wald confidence intervals. Positive values indicate higher activation in
the High Slow Risk group. Both interactions and activation parameters were
estimated separately across groups.}
\label{fig:activation_differences}
\end{figure}
```

### Figure 4

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.82\textwidth]
{figs/revised/figure_04_slowrisk_coefficients.pdf}
\caption{\textbf{Symptom-specific associations with continuous Slow Risk
Load.} Covariate-adjusted models include age, gender, and ethnicity.
Symptom-adjusted models additionally include the other eight depressive
symptoms. Points show log-odds coefficients per 1 SD higher Slow Risk Load and
horizontal lines show 95\% confidence intervals.}
\label{fig:slowrisk_coefficients}
\end{figure}
```

## 6. What the new pipeline replaces

Retire these old main-text outputs:

```text
figs/fig2_refined.pdf
figs/fig3_baseline_activation_joint_globalJ.pdf
figs/fig5_refined.pdf
```

The old shared-\(J\) activation result is retained only as a sensitivity analysis
in Supplementary Figure S2. The old force-directed symptom-map panel is not
recreated because it encouraged a causal interpretation that the regressions do
not support.

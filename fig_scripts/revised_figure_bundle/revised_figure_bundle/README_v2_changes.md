# Figure bundle v2

Replace the existing `fig_scripts/revised/` folder with the folder in this
bundle, then run:

```r
source("fig_scripts/revised/make_revised_figures.R")
```

Changes from v1:

- Figure 1: removed overlapping annotations and simplified the conceptual flow.
- Figure 2: fixed the overlapping histogram labels and clipped Panel B title;
  Panel B is now a horizontal mean-and-CI display.
- Figure 3: fixed the missing suicidality confidence interval by removing the
  destructive scale limit; uses the full "Concentration problems" label.
- Figure 4: shortened the title and replaced technical legend labels with
  "Covariates only" and "Covariates + other symptoms."
- Supplement S2: shortened title.
- Supplement S3: uses readable symptom abbreviations and removes automatic
  "Maximum" annotations.

- Table 2: now uses `tabularx`, a wrapping first column, and the heading `Log likelihood`.

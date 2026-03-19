# Optimal Sample Size Calculation in Cost-Effectiveness Longitudinal Cluster Randomized Trials

[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/haowangbiostat/Cost-effectiveness_LCRT/blob/main/LICENSE)
[![Shiny App](https://img.shields.io/badge/Shiny-Interactive%20App-blue)](https://f07k8s-hao-wang.shinyapps.io/Cost-effectiveness_LCRT/)

This repository contains the R code to reproduce the results presented in [Optimal Sample Size Calculation in Cost-Effectiveness Longitudinal Cluster Randomized Trials](TBD).

### Overview

An illustration of three L-CRT design variants with $J = 6$ periods:

<p align="center">
  <img src="figures/designs.pdf" width="70%">
</p>

Purple cells indicate intervention periods and white cells indicate control periods. Panel (a) displays a PA-LCRT, panel (b) shows a CRXO trial, and panel (c) presents an SW-CRT.

### Quickstart

To reproduce the results, please download this repo on a machine with R, run each R script in the [`codes`](codes) directory without modification, and then the results are saved in [`figures`](figures). To generate tables in the main paper, please run the scripts in [`tables`](tables). All the R scripts can be run standalone. To run the R scripts, you do not need to set any pathnames; everything is relative.

Required R packages: dplyr, ggh4x, ggplot2, nloptr, patchwork, tidyr, xtable.

### Compute Optimal Designs

Run our [`Shiny App`](https://f07k8s-hao-wang.shinyapps.io/Cost-effectiveness_LCRT/) to find optimal sample sizes for cost-effectiveness L-CRTs:

A step-by-step tutorial with worked examples is provided in Web Appendix D of the paper.

### Generate Tables and Figures

#### LODs for CRXO Trials and PA-LCRTs (Table 2)

* Run [`table_2_crxo_pa_lod.R`](codes/table_2_crxo_pa_lod.R)
  + LODs under varying ICC and design parameters for CRXO trials and PA-LCRTs with $J \in \{2, 4, 6\}$ periods
  + generate [`table_2.tex`](tables/table_2.tex)

#### LODs With Varying CAC (Figure 2)

* Run [`figure_2_lod_J4.R`](codes/figure_2_lod_J4.R)
  + LODs for CRXO trials, PA-LCRTs, and SW-CRTs with varying cluster autocorrelation (CAC) when $J = 4$
  + generate [`lod_J4.pdf`](figures/lod_J4.pdf)

#### LODs With Varying Standardized Ceiling Ratio (Web Figure 1)

* Run [`web_figure_1_lod_J4_varying_lambda_r.R`](codes/web_figure_1_lod_J4_varying_lambda_r.R)
  + LODs for CRXO trials, PA-LCRTs, and SW-CRTs with $\lambda r \in \{0.1, 1, 2\}$ when $J = 4$
  + generate [`lod_J4_varying_lambda_r.pdf`](figures/lod_J4_varying_lambda_r.pdf)

#### MMDs for CRXO Trials and PA-LCRTs (Table 3)

* Run [`table_3_crxo_pa_mmd.R`](codes/table_3_crxo_pa_mmd.R)
  + MMDs under varying parameter space specifications for CRXO trials and PA-LCRTs with $J \in \{2, 4, 6\}$ periods
  + generate [`table_3.tex`](tables/table_3.tex)

#### MMDs With Varying Maximum Between-Period Effect ICC (Figure 3)

* Run [`figure_3_mmd_J4.R`](codes/figure_3_mmd_J4.R)
  + MMDs for CRXO trials, PA-LCRTs, and SW-CRTs with varying $\rho_{1,\max}^E$ when $J = 4$
  + generate [`mmd_J4.pdf`](figures/mmd_J4.pdf)

#### LODs and MMDs for SW-CRTs (Tables 4–5)

* Run [`table_4_swcrt_lod.R`](codes/table_4_swcrt_lod.R) and [`table_5_swcrt_mmd.R`](codes/table_5_swcrt_mmd.R)
  + LODs and MMDs for SW-CRTs with $Q \in \{3, 5, 7\}$ treatment sequences
  + generate [`table_4.tex`](tables/table_4.tex) and [`table_5.tex`](tables/table_5.tex)

### Helper Functions

The following scripts in [`codes`](codes) contain helper functions used throughout the analysis. These are automatically sourced by the main scripts — you do not need to run them separately:

* `utils_PA.R`: variance calculation and LOD/MMD optimization for parallel-arm LCRTs
* `utils_CRXO.R`: variance calculation and LOD/MMD optimization for cluster randomized crossover trials
* `utils_SWCRT.R`: variance calculation and LOD/MMD optimization for stepped-wedge CRTs

### Shiny App

The [`shiny_app`](shiny_app) folder contains the source code for the interactive R Shiny application deployed at [https://f07k8s-hao-wang.shinyapps.io/Cost-effectiveness_LCRT/](https://f07k8s-hao-wang.shinyapps.io/Cost-effectiveness_LCRT/). The application automates LOD and MMD calculations for all three L-CRT design variants given user-specified parameters. To run the application locally:

```r
shiny::runApp("shiny_app")
```

### Reference

```
@article{Wang2025costeffectiveness,
  title = {Optimal Sample Size Calculation in Cost-Effectiveness Longitudinal Cluster Randomized Trials},
  journal = {Submitted to Statistics in Medicine},
  author = {Wang, Hao and Liu, Jingxia and Tong, Jiaqi and Cameron, Drew B. and Spiegelman, Donna and Li, Fan},
  year = {2025}
}
```

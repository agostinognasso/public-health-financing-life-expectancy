# Replication material

**Public health financing and life expectancy: Within-country evidence
and implications for planning**

Agostino Gnasso, Girolamo Franchetti, Carmela Iorio, Massimiliano
Politano. University of Naples Federico II.

This folder contains the data extracts, the R code and the output
behind the tables, figures and estimates of the article and of its
Supplementary Material. One command regenerates all of it from the
extracts, with no manual step and no network access.

## Contents

    README.md      this file
    LICENSE.txt    MIT for the code, CC BY 4.0 for data and output
    data/          the six World Bank extracts the analysis starts from
    R/             34 numbered scripts, the shared setup 00_common.R
                   and the driver run_all.R
    output/        estimates and intermediate files, written by the scripts
    figures/       the figures of the article and of the Supplementary
                   Material, under their printed numbers, as PDF and PNG

## Requirements

R 4.6.0 and the CRAN packages dplyr, tidyr, plm, lme4, lmerTest,
rstanarm, posterior, bayesplot, performance, pls, compositions,
zCompositions, changepoint, tsfeatures, WDI, countrycode, MASS,
ggplot2, patchwork, forcats and scales. The version of each package
used is recorded in `output/sessionInfo.txt`, the session information
of the last full run (macOS 26, Apple silicon, R 4.6.0 with OpenBLAS). The full run
takes about 15 minutes on that machine, most of it for the
Bayesian model of script 24.

## Running it

    cd R
    Rscript run_all.R

`run_all.R` runs the scripts in order, each in a fresh R process, and
writes one log per script to a folder `log/`, which it creates. A
script can also be run on its own from `R/`, for example
`Rscript 16_longrun.R`, once the scripts
it reads from have run; the order constraints are listed in
`run_all.R`. `00_common.R` holds the paths, the construction of the
estimation panel (`build_panel()`) and the helpers for long-run
multipliers and Fieller sets, and is sourced by the other scripts.

Random numbers are seeded: `set.seed(42)` in scripts 03, 08 and 23,
`seed = 42` for the Stan sampler in script 24, `set.seed(20260921)`
in script 26. On the same platform and package versions a full run
reproduces the deposited CSV files exactly. On another platform the
deterministic estimators agree to rounding error; the Markov chains of
script 24 and the bootstrap draws of scripts 03, 08, 23 and 26 may
differ in their last digits, because they depend on the compiler and
on the linear algebra library.

Some logs carry warnings that do not stop the run: divergent
transitions in the Stan fits of script 24, which the script counts and
re-examines at adapt_delta = 0.99; PLS fits in some bootstrap draws of
scripts 08 and 23 that do not converge within 100 iterations; warnings
from `pgmm` on the deeper instrument sets of script 31; and the
deprecation of `geom_errorbarh()` in ggplot2 4.0 (script 28).

The number of replications of the Kiviet-Bruno parametric bootstrap
(script 26) can be lowered for a quick check, for example
`LSDVC_BREP=50 Rscript 26_lsdvc.R`; the article uses 500.

## Data

All series come from World Bank Open Data (<https://data.worldbank.org>)
through the R package WDI. The health financing shares originate in
the WHO Global Health Expenditure Database
(<https://apps.who.int/nha/database>) and are republished in the World
Development Indicators. Life expectancy, the driver panel and the two
shock series are from the World Development Indicators release of
8 April 2026 (the `lastupdated` attribute of the extracts); the
level-spending extract was downloaded on 19 September 2026 and the
total-spending extract on 4 November 2026. The financing shares and
current health expenditure per head in US dollars of the 4 November
download are identical, on every observation of the estimation panel,
to those of the 8 April extract, so the health expenditure series come
from the same release of the Global Health Expenditure Database.

| File | Content |
|---|---|
| `df_le_long.rds` | life expectancy at birth, SP.DYN.LE00.IN, 144 countries, 1960-2024 |
| `df_panel_drivers.rds` | the driver panel, 144 countries, 2000-2024 (variables below) |
| `HIV_pct_shock.rds` | HIV prevalence, % of population aged 15-49, SH.DYN.AIDS.ZS, 1995-2024 |
| `HOMICIDE_shock.rds` | intentional homicides per 100,000, VC.IHR.PSRC.P5, 1995-2024 |
| `wdi_level_spending.rds` | public health spending as % of GDP (SH.XPD.GHED.GD.ZS) and per head in current US$ (SH.XPD.GHED.PC.CD), and current health expenditure as % of GDP (SH.XPD.CHEX.GD.ZS), 2000-2024 |
| `wdi_total_spending.rds` | current health expenditure per head, PPP, current international $ (SH.XPD.CHEX.PP.CD, `CHE_PPP`), and in current US$ (SH.XPD.CHEX.PC.CD, `CHE_USD`), 2000-2024; observed for all 2350 country-years of the estimation panel |

Variables of `df_panel_drivers.rds` used by the scripts:

| Variable | Definition | Indicator |
|---|---|---|
| `log_GNI` | log(1 + GNI per head, constant 2015 US$) | NY.GNP.PCAP.KD (`GNI_REAL`) |
| `PUB` | domestic general government health expenditure, % of CHE | SH.XPD.GHED.CH.ZS |
| `PRIV` | domestic private health expenditure, % of CHE | SH.XPD.PVTD.CH.ZS |
| `OOP_share` | out-of-pocket expenditure, % of CHE | SH.XPD.OOPC.CH.ZS |
| `PRIV_NONOOP` | max(PRIV - OOP_share, 0) | |
| `HEALTH_PC` | current health expenditure per head, current US$ | SH.XPD.CHEX.PC.CD |
| `EDU_GROSS` | secondary school enrolment, % gross | SE.SEC.ENRR |
| `IMM_DPT` | DPT immunisation, % of children aged 12-23 months | SH.IMM.IDPT |
| `SANITATION` | people using at least basic sanitation services, % | SH.STA.BASS.ZS |
| `region` | World Bank region | |
| `GNI_label` | World Bank income group | |
| `is_covid_year` | 1 in 2020, 0 otherwise | |

CHE is current health expenditure. The external share is computed in
script 01 as EXT = max(100 - PUB - PRIV, 0). The remaining columns of
the extract are not used. The extracts keep the 144 countries, after
removal of regional and income-group aggregates with the countrycode
dictionary, that have at least 30 non-missing annual values of life
expectancy and at least 15 of income (Section 3.1; script 25 re-applies
the life expectancy rule on 2000-2024). The estimation panel built by
`build_panel()` keeps the country-years with life expectancy and every
driver observed: 2350 observations on 138 countries, 2000-2024.

Scripts 11, 22 and 37 download their series through WDI only if the
extract is missing from `data/`. Deleting an extract therefore
refreshes it from the current World Bank release, which will change
the results.

## Output

Each script writes its results to `output/`, in files whose name
starts with the script number from script 16 on. Two fitted model
objects are kept, `panel_lmer_model.rds` (script 04) and
`panel_ilr_model.rds` (script 08), because script 32 reads them; no
other script saves a fitted model, so a full run writes exactly the
files deposited here. Script 33 prints the LaTeX bodies of Tables 1
and 3 to the console; Supplementary Table S1 is written to
`output/35_hiv_table.tex`, the body of Table 4 to
`output/37_total_spending_table.tex` and the body of Table 5 to
`output/39_planning_table.tex`. `output/34_headline.csv` collects the
numbers the article quotes most often, and `output/34_samples.csv` the
observations and countries behind each analysis.

## What produces what

| Article element | Script | Main output |
|---|---|---|
| Table 1, descriptive statistics | 33; 29, 19 | console output of 33, `29_table1_full.csv`, `19_within_between_sd.csv` |
| Table 2, the estimator ladder | 28, from 20, 21, 26 | `28_ladder_impact.csv`, `28_ladder_longrun.csv` |
| Table 3, System GMM and long-run multipliers | 16, 33 | `16_sysgmm_shortrun_ci.csv`, `16_sysgmm_longrun.csv` |
| Table 4 and Section 6.1, with and without total health expenditure per head | 37 | `37_total_spending.csv`, `37_total_spending_diag.csv`, `37_total_spending_table.tex`; country-years of each sample in `37_total_spending_rows.csv` |
| Table 5 and Section 7.1, a planned ten-point shift in the public share | 39, from 16, 21, 26, 35, 37, 38 | `39_planning.csv`, `39_planning_table.tex` |
| Figure 1, mixed-model coefficients | 28 | `figures/Figure1`, `28_mixed_coefficients.csv` |
| Figure 2, the ladder on two estimands | 28, 26 | `figures/Figure2` |
| Section 3.1, sample and regions | 34, 25 | `34_samples.csv` |
| Section 3.2, from ILR1 to percentage points | 17 | `17_ilr_interpretive_map.csv` |
| Section 5.1, within-country t statistics, ICC | 29 | `29_within_tstats.csv` |
| Sections 5.2 and 5.3, long run, MDE, equivalence | 16, 19 | `19_mde.csv`, `19_mde_longrun.csv`, `19_tost.csv` |
| Sections 5.1 and 5.2, distributed lag, exogeneity classification | 21 | `21_exogeneity.csv`, `21_dlag_cumulative.csv` |
| Section 5.3, Kiviet-Bruno correction | 26 | `26_lsdvc.csv`, `26_lsdvc_bootstrap.csv` |
| Section 5.2, what the year effects take out | 36 | `36_fwl.csv`, `36_fwl_periods.csv` |
| Sections 5.2 and 6.2, the AR(3) rejection | 31 | `31_ar3_*.csv` |
| Section 5.4, the public share by income group in 2000 | 38 | `38_income.csv`, `38_income_groups.csv` |
| Section 6.2, the pandemic years | 20 | `20_ladder.csv`, `20_dlag_windows.csv` |
| Section 6.2, zero treatment | 18 | `18_delta_sensitivity.csv` |
| Section 6.1, measure of public spending | 22, 21 | `22_level_ladder.csv` |
| Section 6.2, missingness and gaps | 25 | `25_*.csv` |
| Section 6.2, HIV in the sample and in the model | 35 | `35_hiv*.csv` |
| Section 6.2, the implausible outcome value | 30 | `30_outlier_influence.csv` |
| Section 6.2, dynamic panel diagnostics | 15, 31 | `gmm_robustness_diagnostics.csv` |
| Supplementary Section S1.1, shape features | 02 | `shape_features_post.rds`, `shape_pca_*.csv` |
| Supplementary Section S1.2, predictive design | 03 | `pls/` |
| Supplementary Section S1.3, robustness of the design | 09 | `sensitivity_window.csv`, `influence_pls_loo.csv` |
| Supplementary Section S1.4, level against shape | 03, 08, 23 | `pls/VIP_*.csv`, `23_vip_bootstrap.csv` |
| Supplementary Section S2, Bayesian model and regional transfer | 24, 23 | `24_*.csv`, `23_loro_decomposed.csv` |
| Supplementary Section S4, Table S1 | 35 | `35_hiv_table.tex` |
| Supplementary Section S5, clipping and zeros in the simplex | 18 | `18_zero_audit.csv` |
| Supplementary Section S5, correlations and information criteria | 32 | `32_simplex_ic.csv` |
| Supplementary Section S6 and Table S2, cross-sectional R2 with and without HIV | 27 | `27_table2.csv` |
| Supplementary Figures S1-S4 | 05 | `figures/FigureS1` to `FigureS4` |
| Supplementary Figures S5, S6 | 10 | `figures/FigureS5`, `FigureS6` |
| Supplementary Figures S7-S12 | 28 | `figures/FigureS7` to `FigureS12` |

Scripts 01, 04 and 11 build inputs for the others: the panel with the
external share and the cross-sectional targets (01), the raw-share
mixed model that script 32 compares with the ILR model (04), and the
HIV and homicide series (11).

## The Kiviet-Bruno correction

`26_lsdvc.R` implements the bias-corrected LSDV estimator of Kiviet
(1995), extended to unbalanced panels by Bruno (2005), which Stata
provides as `xtlsdvc` and for which we found no CRAN package. The
approximation formulas are written out in the script header and coded
directly, and the script checks them by simulation before applying
them: on a simulated panel of the same size (N = 138, T = 25,
gamma = 0.87) the correction brings the LSDV bias in the persistence
parameter from -0.072 to -0.00004, and from -0.076 to +0.0001
with 32 per cent of cells missing. The uncorrected LSDV agrees with
`plm` to 2e-15.

On the data it gives alpha = 0.8932 with a parametric bootstrap
standard error of 0.0155 (Bun and Kiviet 2001, 500 replications),
against 0.8129 uncorrected, and beta(PUB) = -0.00098 (SE 0.00394,
p = 0.80).

## Citation

Please cite the article and this deposit. The citation metadata of the
deposit are in `CITATION.cff`; each release is archived on Zenodo with
its own DOI.

## Licence

See `LICENSE.txt`.

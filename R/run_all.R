# Runs the whole pipeline in order:
#
#   cd R
#   Rscript run_all.R
#
# Each script runs in its own R process. Results go to ../output,
# figures to ../figures, and the console output of each script to
# ../log/<script>.log. At the end the R session information, with
# the version of every package the scripts use, is written to
# ../output/sessionInfo.txt.
#
# No network access is needed: scripts 11, 22 and 37 read World Bank data
# through the WDI package only if their extract is missing from
# ../data.
#
# R messages are produced in English (LANGUAGE=en) whatever the
# locale of the machine; this affects the logs only.

if (!file.exists("00_common.R"))
  stop("run this from the R/ directory: Rscript run_all.R")

dir.create(file.path("..", "log"), showWarnings = FALSE, recursive = TRUE)

STEPS <- c(
  "01_data.R"                , # panel with the external share, cross-sectional targets
  "02_shape_post_pca.R"      , # shape features of the 2005-2024 trajectory, PCA
  "03_pls_predictive.R"      , # predictive PLS, regional holdout, paired VIP
  "04_panel_mixed.R"         , # mixed model on the raw shares
  "05_plots.R"               , # Supplementary Figures S1-S4
  "08_compositional_ilr.R"   , # ILR coordinates: PLS and mixed model
  "09_robustness.R"          , # initial window, influential countries
  "10_publication_plots.R"   , # Supplementary Figures S5, S6
  "11_shock_variables.R"     , # HIV prevalence and homicide rate
  "15_gmm_robustness.R"      , # System GMM with and without year effects
  "16_longrun.R"             , # Table 4, long-run multipliers, Fieller sets
  "17_ilr_sysgmm.R"          , # System GMM on the ILR coordinates
  "18_zeros.R"               , # zeros in the simplex, detection limit
  "19_inference.R"           , # within/between SD, MDE, equivalence tests
  "20_ladder_windows.R"      , # estimator ladder, 2000-2024 and 2000-2019
  "21_exogeneity.R"          , # exogeneity classification, longer lags
  "22_level_spending.R"      , # spending as % of GDP and per head
  "23_loro_pls.R"            , # leave-one-region-out, VIP bootstrap
  "24_bayes.R"               , # Bayesian hierarchical model
  "25_missingness.R"         , # missingness, gaps, inclusion rule
  "26_lsdvc.R"               , # Kiviet-Bruno bias correction
  "27_table2.R"              , # Table 2
  "28_figures.R"             , # Figures 1, 2 and S7-S12
  "29_text_numbers.R"        , # numbers quoted in the text, Table 1 body
  "30_outlier.R"             , # Central African Republic 2009
  "31_ar3.R"                 , # the AR(3) rejection
  "32_simplex_evidence.R"    , # information criteria, raw shares against ILR
  "33_latex_tables.R"        , # LaTeX for Tables 1 and 4
  "35_hiv.R"                 , # HIV samples A, B, C; Supplementary Table S1
  "36_fwl.R"                 , # what the year effects take out of static FE
  "37_total_spending.R"      , # total health spending per head; Table 4
  "38_income.R"              , # public share by income group in 2000
  "39_planning.R"            , # a planned ten-point shift; Table 5
  "34_summary.R"               # headline numbers and samples, from the above
)

# Order constraints: 01 before everything; 02 before 03, 08, 09; 03
# before 05; 04 and 08 before 32; 08 and 09 before 10; 11 before 15,
# 24, 27, 29, 33, 34 and 35; 16 before 19, 28 and 29; 18, 19, 23, 24
# and 26 before 28; 28 before 29; 34 last. 37 and 38 need only 01;
# 39 reads the output of 16, 21, 26, 35, 37 and 38.

t0 <- Sys.time()
for (f in STEPS) {
  cat("\n---", f, "\n")
  ti  <- Sys.time()
  log <- file.path("..", "log", sub("\\.R$", ".log", f))
  st  <- system2("Rscript", f, stdout = log, stderr = log,
                 env = "LANGUAGE=en")
  cat(sprintf("%-28s %-6s %5.1f min\n", f,
              if (st == 0) "ok" else paste0("exit ", st),
              as.numeric(difftime(Sys.time(), ti, units = "mins"))))
}
cat(sprintf("\nTotal: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# session information, with every package the scripts load
PKGS <- c("dplyr", "tidyr", "plm", "lme4", "lmerTest", "rstanarm",
          "posterior", "bayesplot", "performance", "pls", "compositions",
          "zCompositions", "changepoint", "tsfeatures", "WDI",
          "countrycode", "MASS", "ggplot2", "patchwork", "forcats",
          "scales")
suppressPackageStartupMessages(
  for (p in PKGS) library(p, character.only = TRUE))
writeLines(c(sprintf("Pipeline run of %s", format(t0, "%Y-%m-%d %H:%M %Z")),
             "", capture.output(sessionInfo())),
           file.path("..", "output", "sessionInfo.txt"))

cat("Results in ../output, figures in ../figures, logs in ../log\n")

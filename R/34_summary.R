# 34: headline numbers and the sample behind each analysis
#
# Collects from the outputs of the other scripts the numbers the
# article quotes most often, so that text and estimates can be checked
# in one place, and lists the observations and countries behind each
# analysis. Runs last.
#
# Reads  output/ (the files named below)
# Writes output/34_headline.rds, 34_headline.csv, 34_samples.csv
source("00_common.R")
rd <- function(f) { p <- file.path(OUT, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }

S <- list()
sr <- rd("16_sysgmm_shortrun_ci.csv"); lr <- rd("16_sysgmm_longrun.csv")
S$alpha      <- sr$beta[sr$term == "lag(LE, 1)"]
S$alpha_se   <- sr$se[sr$term == "lag(LE, 1)"]
S$pub_sr     <- sr$beta[sr$term == "PUB"];      S$pub_sr_se <- sr$se[sr$term == "PUB"]
S$imm_sr     <- sr$beta[sr$term == "IMM_DPT"];  S$imm_sr_se <- sr$se[sr$term == "IMM_DPT"]
# long-run multipliers: Fieller sets only (script 16)
pl <- lr[lr$term == "PUB", ]; il <- lr[lr$term == "IMM_DPT", ]
S$pub_lr              <- pl$beta_LR
S$pub_lr_f90          <- c(pl$fieller90_lo, pl$fieller90_hi)
S$pub_lr_f95          <- pl$fieller_type
S$fieller_bounded_to  <- pl$level_bounded
S$pub_wholeline_from  <- pl$level_wholeline
S$imm_lr              <- il$beta_LR
S$imm_lr_f90          <- c(il$fieller90_lo, il$fieller90_hi)
S$imm_lr_f95          <- c(il$fieller_lo, il$fieller_hi)
S$imm_lr_f95_type     <- il$fieller_type

lad <- rd("28_ladder_impact.csv")
g <- function(p) lad$beta[grepl(p, lad$rung)]
S$fe_oneway <- g("^1"); S$fe_twoway <- g("^2")
S$share_year_effects <- (S$fe_twoway - S$fe_oneway) / (S$pub_sr - S$fe_oneway)

sd  <- rd("19_within_between_sd.csv")
S$sd_within_pub <- sd$sd_within[sd$variable == "Public share (% CHE)"]
mde <- rd("19_mde.csv")
# long-run MDE as the impact MDE over 1 - alpha (script 19)
mlr <- rd("19_mde_longrun.csv")
S$mde_sr <- mde$mde80[1]; S$mde_lr <- mde$mde80_exact[2]
S$mde_lr_alpha90 <- mlr$mde80[mlr$alpha_at %in% c("90% CI lower", "90% CI upper")]
S$mde_lr_min95   <- mlr$mde80[mlr$alpha_at == "95% CI lower"]
tt <- rd("19_tost.csv")
S$tost_sr_min <- suppressWarnings(min(tt$margin[tt$horizon == "short run (impact)" & tt$equivalent]))
S$tost_lr_any <- any(tt$equivalent[tt$method == "Fieller"])

R5 <- rd("20_ladder.csv")
pre <- R5[R5$window == "2000-2019" & grepl("^5", R5$rung), ]
S$pre_pub <- pre$beta; S$pre_se <- pre$se; S$pre_p <- pre$p; S$pre_alpha <- pre$alpha

ex <- rd("21_exogeneity.csv")
S$pub_endog <- ex$beta[ex$term == "PUB" & grepl("^C ", ex$spec)]
S$imm_endog <- ex$beta[ex$term == "IMM_DPT" & grepl("^C ", ex$spec)]
S$imm_endog_p <- ex$p[ex$term == "IMM_DPT" & grepl("^C ", ex$spec)]

zs <- rd("18_delta_sensitivity.csv")
S$ilr1_min <- min(zs$ILR1); S$ilr1_max <- max(zs$ILR1)
S$ilr1_ever_sig <- any(zs$p_ILR1 < 0.05)
S$imm_always_sig <- all(zs$p_IMM < 0.05)

lo <- rd("23_loro_decomposed.csv")
S$loro_neg_total <- sum(lo$R2_total < 0); S$loro_n <- nrow(lo)
S$loro_neg_given <- sum(lo$R2_given_region_mean < 0)
S$loro_pct_level <- median(lo$pct_mse_from_level)
S$loro_cor <- mean(lo$cor_pred_true, na.rm = TRUE)

bt <- rd("23_vip_bootstrap.csv")
S$vip_w_iid <- mean(bt$width[grepl("^A", bt$scheme)])
S$vip_w_clu <- mean(bt$width[grepl("^B", bt$scheme)])
S$vip_w_two <- mean(bt$width[grepl("^C", bt$scheme)])
S$vip_any_sig <- any(bt$excludes_zero)

cv <- rd("24_convergence.csv")
if (!is.null(cv)) { S$rhat_max <- max(cv$rhat, na.rm = TRUE)
                    S$ess_min  <- min(cv$ess_bulk, na.rm = TRUE) }
ad <- rd("24_adapt_delta.csv")
if (!is.null(ad)) { S$divergent <- ad$divergent
                    S$adapt_shift_sd <- ad$max_shift_in_sd[2] }

S$zero_alpha_range   <- range(zs$alpha)
S$zero_bounded_range <- range(zs$level_bounded)

# HIV: samples A (analytical panel), C (HIV reported, with HIV)
hv <- rd("35_hiv.csv"); hs <- rd("35_hiv_samples.csv")
gA <- hv[hv$sample == "A" & hv$estimator == "System GMM" & hv$term == "PUB", ]
gC <- hv[hv$sample == "C" & hv$estimator == "System GMM" & hv$term == "PUB", ]
S$hiv_pub_sr     <- c(gC$estimate, gC$se)
S$hiv_shift_se   <- (gC$estimate - gA$estimate) / gA$se
S$hiv_alpha      <- hs$alpha[hs$sample == "C"]
S$hiv_bounded_to <- hs$level_bounded[hs$sample == "C"]
S$hiv_within_fe  <- unlist(hv[hv$sample == "C" & hv$term == "HIV_pct" &
                              hv$estimator == "FE, country effects", c("estimate", "se")])
S$hiv_within_mm  <- unlist(hv[hv$sample == "C" & hv$term == "HIV_pct_wn",
                              c("estimate", "se")])

# what the year effects take out of the static estimate (script 36)
fw <- rd("36_fwl_periods.csv")
S$fwl_by_period <- setNames(fw$contribution[match(c("2000-2009", "2010-2019", "2020-2024"),
                                                  fw$period)],
                            c("2000-09", "2010-19", "2020-24"))
S$fwl_pre_move  <- fw$contribution[fw$period == "re-estimated on 2000-2019"]

# distributed lag, lags 0 to 3 (script 21)
dlc <- rd("21_dlag_cumulative.csv"); dlc <- dlc[dlc$lags == "0-3", ]
S$dl_cumulative <- setNames(dlc$cumulative, dlc$driver)
S$dl_p          <- setNames(dlc$p, dlc$driver)

cat("=== Headline numbers ===\n\n")
f <- function(k, fmt, ...) cat(sprintf(paste0("%-46s", fmt, "\n"), k, ...))
cat("-- persistence and the long run --\n")
f("persistence alpha", "%.4f (SE %.4f)", S$alpha, S$alpha_se)
f("1 - alpha", "%.4f, 95%% CI [%.4f, %.4f]", 1-S$alpha,
  1-S$alpha-1.96*S$alpha_se, 1-S$alpha+1.96*S$alpha_se)
f("long-run multiplier on the public share", "%+.4f", S$pub_lr)
f("  its 90% Fieller interval", "[%+.4f, %+.4f]", S$pub_lr_f90[1], S$pub_lr_f90[2])
f("  its 95% Fieller set", "%s", S$pub_lr_f95)
f("Fieller sets are intervals up to", "%.1f%%", 100*S$fieller_bounded_to)
f("  public share: whole line from", "%.2f%%", 100*S$pub_wholeline_from)
f("alpha across zero treatments (ILR model)", "%.3f to %.3f", S$zero_alpha_range[1], S$zero_alpha_range[2])
f("  sets are intervals up to", "%.1f%% to %.1f%%", 100*S$zero_bounded_range[1], 100*S$zero_bounded_range[2])
cat("\n-- static and dynamic estimates on the public share --\n")
f("static FE, country effects only", "%+.4f", S$fe_oneway)
f("static FE, country and year effects", "%+.4f", S$fe_twoway)
f("System GMM", "%+.4f", S$pub_sr)
f("share of the move due to year effects alone", "%.0f%%", 100*S$share_year_effects)
cat("\n-- what the design can and cannot detect --\n")
f("within-country SD of the public share", "%.2f pp", S$sd_within_pub)
f("1-SD move at the one-way FE estimate", "%.3f years", S$fe_oneway*S$sd_within_pub)
f("MDE (80% power), impact", "%.4f years/pp", S$mde_sr)
f("MDE (80% power), long run, at alpha-hat", "%.4f years/pp", S$mde_lr)
f("  over the 90% interval for alpha", "%.4f to %.4f", S$mde_lr_alpha90[1], S$mde_lr_alpha90[2])
f("  over the 95% interval for alpha", "at least %.4f, no upper bound", S$mde_lr_min95)
f("TOST equivalence, impact", "established from margin %.4f", S$tost_sr_min)
f("TOST equivalence, long run (Fieller)", "%s", ifelse(S$tost_lr_any, "established", "not established at any margin"))
cat("\n-- robustness --\n")
f("System GMM 2000-2019", "%+.4f (SE %.4f, p %.3f)", S$pre_pub, S$pre_se, S$pre_p)
f("public share, drivers treated as endogenous", "%+.4f", S$pub_endog)
f("ILR1 across all zero treatments", "[%+.4f, %+.4f]", S$ilr1_min, S$ilr1_max)
f("  significant at 5% under any treatment", "%s", ifelse(S$ilr1_ever_sig, "yes", "no"))
cat("\n-- immunisation --\n")
f("impact effect", "%+.4f (SE %.4f)", S$imm_sr, S$imm_sr_se)
f("long-run multiplier", "%+.4f", S$imm_lr)
f("  90% Fieller interval", "[%+.4f, %+.4f]", S$imm_lr_f90[1], S$imm_lr_f90[2])
f("  95% Fieller set", "%s", fmt_fieller(S$imm_lr_f95[1], S$imm_lr_f95[2], S$imm_lr_f95_type))
f("  over 10 pp of coverage", "%+.2f years; 90%% [%+.2f, %+.2f]", 10*S$imm_lr,
  10*S$imm_lr_f90[1], 10*S$imm_lr_f90[2])
f("  over 10 pp, 95% set", "%s", fmt_fieller(10*S$imm_lr_f95[1], 10*S$imm_lr_f95[2],
                                             S$imm_lr_f95_type, 2))
f("significant under every zero rule", "%s", ifelse(S$imm_always_sig, "yes", "no"))
f("with drivers treated as endogenous", "%+.4f (p %.3f)", S$imm_endog, S$imm_endog_p)
cat("\n-- regional transfer and level-vs-shape --\n")
f("regions with negative out-of-region R2", "%d of %d", S$loro_neg_total, S$loro_n)
f("  still negative after granting the mean", "%d of %d", S$loro_neg_given, S$loro_n)
f("  median share of error that is level", "%.0f%%", S$loro_pct_level)
f("  mean within-region correlation", "%+.3f", S$loro_cor)
f("VIP interval width, countries iid", "%.2f", S$vip_w_iid)
f("VIP interval width, region cluster", "%.2f", S$vip_w_clu)
f("VIP interval width, two-stage", "%.2f", S$vip_w_two)
f("  any difference excluding zero", "%s", ifelse(S$vip_any_sig, "yes", "no"))
if (!is.null(cv)) { cat("\n-- Bayesian --\n")
  f("worst R-hat", "%.4f", S$rhat_max); f("smallest bulk ESS", "%.0f", S$ess_min) }
if (!is.null(ad)) {
  f("divergences at adapt_delta 0.95 / 0.99", "%d / %d", S$divergent[1], S$divergent[2])
  f("  largest shift in a posterior mean", "%.2f posterior SD", S$adapt_shift_sd) }
cat("\n-- HIV in the model (sample C against A) --\n")
f("public share, impact, with HIV", "%+.4f (SE %.4f)", S$hiv_pub_sr[1], S$hiv_pub_sr[2])
f("  shift from the analytical panel", "%+.2f SE", S$hiv_shift_se)
f("  alpha; Fieller sets intervals up to", "%.3f; %.1f%%", S$hiv_alpha, 100*S$hiv_bounded_to)
f("HIV within, FE clustered", "%+.2f (SE %.2f)", S$hiv_within_fe[1], S$hiv_within_fe[2])
f("HIV within, mixed model", "%+.2f (SE %.3f)", S$hiv_within_mm[1], S$hiv_within_mm[2])
cat("\n-- what the year effects take out of the static estimate --\n")
f("by period, year effects from their mean", "%s", paste(sprintf("%s %+.4f",
  names(S$fwl_by_period), S$fwl_by_period), collapse = ", "))
f("move re-estimated on 2000-2019", "%+.4f", S$fwl_pre_move)
cat("\n-- distributed lag, lags 0 to 3 --\n")
for (v in names(S$dl_cumulative))
  f(paste(" ", v), "%+.4f (p %.3f)", S$dl_cumulative[[v]], S$dl_p[[v]])

saveRDS(S, file.path(OUT, "34_headline.rds"))
write.csv(data.frame(quantity = names(S),
                     value = sapply(S, function(x) paste(format(x, digits = 5), collapse = "; "))),
          file.path(OUT, "34_headline.csv"), row.names = FALSE)

# ---- the sample behind each analysis ------------------------
pan <- build_panel(year_max = 2024, zero_rule = "multiplicative")
pre <- build_panel(year_max = 2019, zero_rule = "multiplicative")
hiv <- pan %>%
  inner_join(readRDS(file.path(OUT, "shocks_panel.rds")) %>%
               filter(!is.na(HIV_pct)) %>% dplyr::select(iso3c, year),
             by = c("iso3c", "year")) %>%
  filter(!is.na(region), region != "")
dl <- rd("21_dlag_cumulative.csv"); dl <- dl[dl$lags == "0-3", ][1, ]
gp <- rd("25_gap_sensitivity.csv"); gp <- gp[grepl("contiguous", gp$sample), ]
sw <- rd("sensitivity_window.csv")
t2 <- rd("27_table2.csv")

samp <- data.frame(
  analysis = c(
    "Static FE, LSDV, Kiviet-Bruno, difference and System GMM; mixed model; System GMM on ILR coordinates",
    "Distributed lag, lags 0 to 3",
    "Pre-pandemic window",
    "Longest contiguous run per country",
    "Bayesian hierarchical model; HIV in the panel models",
    "Predictive PLS, leave-one-region-out, paired VIP bootstrap",
    "Initial-window sensitivity, shortest window",
    "Initial-window sensitivity, longest window",
    "Cross-sectional R2 with HIV and homicide (Supplementary Table S2)"),
  design = c(rep("panel", 5), rep("cross-section", 4)),
  years = c("2000-2024", "2000-2024", "2000-2019", "2000-2024", "2000-2024",
            "initial window 2000-2005", "initial window 2000-2003",
            "initial window 2000-2007", "initial window 2000-2005"),
  n_obs = c(nrow(pan), dl$n_obs, nrow(pre), gp$n_obs, nrow(hiv),
            NA, NA, NA, NA),
  n_countries = c(n_distinct(pan$iso3c), dl$n_countries,
                  n_distinct(pre$iso3c), gp$n_countries, n_distinct(hiv$iso3c),
                  lo$n_test[1] + lo$n_train[1], sw$n_countries[sw$window == "W3"],
                  sw$n_countries[sw$window == "W7"], t2$n[1]),
  scripts = c("16-20, 26, 28", "21", "20, 36", "25", "24, 35",
              "03, 23", "09", "09", "27"))
stopifnot(sw$n_countries[sw$window == "W5"] == lo$n_test[1] + lo$n_train[1])

cat("\n-- the sample behind each analysis --\n")
print(samp %>% dplyr::select(design, years, n_obs, n_countries, scripts, analysis),
      row.names = FALSE, right = FALSE)
write.csv(samp, file.path(OUT, "34_samples.csv"), row.names = FALSE)

cat("\n-- countries per World Bank region, analytical panel --\n")
reg <- pan %>% dplyr::distinct(iso3c, region) %>%
  dplyr::count(region, name = "countries") %>% dplyr::arrange(dplyr::desc(countries))
print(as.data.frame(reg), row.names = FALSE)
cat(sprintf("total: %d countries\n", sum(reg$countries)))
cat("income groups represented (World Bank classification, any year):",
    paste(sort(unique(na.omit(pan$GNI_label))), collapse = ", "), "\n")

cat("\nSaved: 34_headline.rds, 34_headline.csv, 34_samples.csv\n")

# 35: HIV prevalence, in the sample and in the model (Sections 5.2
#     and 6.5, Supplementary Table S1)
#
# HIV prevalence is reported for fewer countries than the other
# drivers, so adding it as a regressor changes two things at once:
# the sample and the model. Three samples separate them:
#   A  the analytical panel of build_panel(), without HIV
#   B  the observations of A where HIV is reported, without HIV
#   C  the same observations as B, with HIV as a regressor
# A to B is the change of sample, B to C the change of model.
#
# Each sample is fitted with the System GMM of Table 3, the static
# within estimator with country effects and with country and year
# effects, and the mixed model of Figure 1. For the public share the
# script also repeats the equivalence tests of script 19 and the
# Fieller sets of script 16.
#
# Reads  output/shocks_panel.rds; data/ through build_panel()
# Writes output/35_hiv.csv, 35_hiv_samples.csv, 35_hiv_equivalence.csv,
#        35_hiv_table.tex

source("00_common.R")
suppressPackageStartupMessages({ library(plm); library(lme4) })

shk  <- readRDS(file.path(OUT, "shocks_panel.rds")) %>%
  dplyr::select(iso3c, year, HIV_pct)
base <- build_panel(year_max = 2024, zero_rule = "multiplicative") %>%
  left_join(shk, by = c("iso3c", "year"))

samples <- list(
  A = list(dat = base, hiv = FALSE,
           label = "A  analytical panel, without HIV"),
  B = list(dat = base %>% filter(!is.na(HIV_pct)), hiv = FALSE,
           label = "B  HIV reported, without HIV"),
  C = list(dat = base %>% filter(!is.na(HIV_pct)), hiv = TRUE,
           label = "C  HIV reported, with HIV"))

# ---- 1. what requiring HIV removes --------------------------
cat("=== Samples ===\n")
for (s in names(samples)) {
  d <- samples[[s]]$dat
  cat(sprintf("%-36s %5d obs, %3d countries\n", samples[[s]]$label,
              nrow(d), n_distinct(d$iso3c)))
}
cov_hiv <- base %>% group_by(iso3c) %>%
  summarise(n = n(), n_hiv = sum(!is.na(HIV_pct)), .groups = "drop")
gone <- cov_hiv %>% filter(n_hiv == 0)
part <- cov_hiv %>% filter(n_hiv > 0, n_hiv < n)
cat(sprintf("\nCountries with no HIV value: %d (%d obs)\n",
            nrow(gone), sum(gone$n)))
cat(" ", paste(gone$iso3c, collapse = " "), "\n")
cat(sprintf("Countries with some years missing: %d (%d obs)\n",
            nrow(part), sum(part$n - part$n_hiv)))

consec <- function(d) d %>% group_by(iso3c) %>% arrange(year, .by_group = TRUE) %>%
  summarise(p = sum(diff(year) == 1), g = sum(diff(year) > 1),
            .groups = "drop")
for (s in c("A", "B")) {
  k <- consec(samples[[s]]$dat)
  cat(sprintf("Sample %s: %d consecutive-year pairs, %d internal gaps\n",
              s, sum(k$p), sum(k$g)))
}
hw <- samples$C$dat %>% group_by(iso3c) %>%
  mutate(w = HIV_pct - mean(HIV_pct)) %>% ungroup()
cat(sprintf("HIV prevalence in C: mean %.2f, SD %.2f, within-country SD %.2f pp\n",
            mean(hw$HIV_pct), sd(hw$HIV_pct), sd(hw$w)))

# ---- 2. the four estimators ---------------------------------
tidy <- function(est, se, est_name, s) {
  data.frame(sample = s, estimator = est_name, term = names(est),
             estimate = unname(est), se = unname(se),
             z = unname(est / se), p = unname(2 * pnorm(-abs(est / se))),
             row.names = NULL)
}

fit_gmm <- function(d, X, s) {
  pd <- pdata.frame(d, index = c("iso3c", "year"))
  f  <- as.formula(paste("LE ~ lag(LE, 1) +", paste(X, collapse = " + "),
                         "| lag(LE, 2:5)"))
  m  <- pgmm(f, data = pd, effect = "twoways", model = "twosteps",
             transformation = "ld", collapse = TRUE)
  b  <- coef(m); V <- pgmm_vcov(m)
  keep <- c("lag(LE, 1)", X)
  co <- tidy(b[keep], sqrt(diag(V)[keep]), "System GMM", s)

  lr   <- longrun(b, V, alpha_name = "lag(LE, 1)")
  lr90 <- longrun(b, V, alpha_name = "lag(LE, 1)", level = 0.90)
  lv   <- fieller_levels(b, V, "lag(LE, 1)")
  lr <- data.frame(term = lr$term, beta_LR = lr$beta_LR,
                   fieller_lo = lr$fieller_lo, fieller_hi = lr$fieller_hi,
                   fieller_type = lr$fieller_type,
                   fieller90_lo = lr90$fieller_lo,
                   fieller90_hi = lr90$fieller_hi,
                   fieller90_type = lr90$fieller_type,
                   level_wholeline = lv$level_wholeline[match(lr$term, lv$term)])
  co <- left_join(co, lr, by = "term")

  sg  <- sargan(m)
  dg <- data.frame(sample = s, alpha = unname(b["lag(LE, 1)"]),
                   se_alpha = sqrt(V["lag(LE, 1)", "lag(LE, 1)"]),
                   level_bounded = lv$level_bounded[1],
                   n_instr = as.numeric(sg$parameter) + length(b),
                   hansen_p = as.numeric(sg$p.value),
                   ar1_p = as.numeric(mtest(m, order = 1)$p.value),
                   ar2_p = as.numeric(mtest(m, order = 2)$p.value))
  list(co = co, diag = dg, b = b, V = V)
}

fit_fe <- function(d, X, s, effect) {
  pd <- pdata.frame(d, index = c("iso3c", "year"))
  m  <- plm(as.formula(paste("LE ~", paste(X, collapse = " + "))),
            data = pd, model = "within", effect = effect)
  V  <- plm::vcovHC(m, cluster = "group", method = "arellano")
  lab <- if (effect == "individual") "FE, country effects" else
    "FE, country and year effects"
  tidy(coef(m)[X], sqrt(diag(V)[X]), lab, s)
}

fit_mixed <- function(d, hiv, s) {
  bw <- c("log_GNI", "EDU_GROSS", "IMM_DPT", "SANITATION",
          "ILR1", "ILR2", "ILR3", if (hiv) "HIV_pct")
  dm <- d %>% group_by(iso3c) %>%
    mutate(across(all_of(bw), list(bw = ~ mean(.x), wn = ~ .x - mean(.x)))) %>%
    ungroup() %>% mutate(year_c = year - 2012)
  rhs <- paste(c("year_c", "is_covid_year",
                 as.vector(rbind(paste0(bw, "_bw"), paste0(bw, "_wn"))),
                 "(1 | iso3c)"), collapse = " + ")
  m  <- lmer(as.formula(paste("LE ~", rhs)), data = dm, REML = TRUE,
             control = lmerControl(calc.derivs = FALSE))
  ss <- summary(m)$coefficients
  ss <- ss[grepl("_(bw|wn)$", rownames(ss)), , drop = FALSE]
  tidy(setNames(ss[, 1], rownames(ss)), ss[, 2], "mixed model", s)
}

co <- list(); dg <- list(); gm <- list()
for (s in names(samples)) {
  d <- samples[[s]]$dat
  X <- c(KEY_X, if (samples[[s]]$hiv) "HIV_pct")
  g <- fit_gmm(d, X, s)
  gm[[s]] <- g
  dg[[s]] <- g$diag %>% mutate(n_obs = nrow(d),
                               n_countries = n_distinct(d$iso3c), .after = 1)
  co[[s]] <- bind_rows(g$co,
                       fit_fe(d, X, s, "individual"),
                       fit_fe(d, X, s, "twoways"),
                       fit_mixed(d, samples[[s]]$hiv, s))
}
co <- bind_rows(co)
dg <- bind_rows(dg)

cat("\n=== Coefficients by sample: estimate (SE) ===\n")
show <- function(est, tm) {
  r <- co %>% filter(estimator == est, term == tm)
  v <- setNames(sprintf("%+.4f (%.4f)", r$estimate, r$se), r$sample)
  cat(sprintf("%-30s %-13s %-18s %-18s %-18s\n", est, tm,
              ifelse(is.na(v["A"]), "", v["A"]), ifelse(is.na(v["B"]), "", v["B"]),
              ifelse(is.na(v["C"]), "", v["C"])))
}
cat(sprintf("%-30s %-13s %-18s %-18s %-18s\n", "", "", "A", "B", "C"))
for (tm in c("lag(LE, 1)", "PUB", "OOP_share", "IMM_DPT", "HIV_pct"))
  show("System GMM", tm)
for (est in c("FE, country effects", "FE, country and year effects"))
  for (tm in c("PUB", "HIV_pct")) show(est, tm)
for (tm in c("ILR1_wn", "ILR1_bw", "IMM_DPT_wn", "HIV_pct_wn", "HIV_pct_bw"))
  show("mixed model", tm)

cat("\n=== System GMM diagnostics ===\n")
print(dg %>% mutate(across(c(alpha, se_alpha, level_bounded, hansen_p,
                             ar1_p, ar2_p), ~ round(.x, 4))))

# ---- 3. the public share: long run and equivalence ----------
# margins as in script 19; 0.0275 is the static within estimate
# with country effects on the analytical panel
margins <- c(0.0125, 0.025, 0.0275, 0.05)
tost_p <- function(b, se, m) max(1 - pnorm((b + m) / se), pnorm((b - m) / se))
tost_fieller_p <- function(b, V, m) {
  den <- 1 - b["lag(LE, 1)"]; vb <- V["PUB", "PUB"]
  va  <- V["lag(LE, 1)", "lag(LE, 1)"]; cab <- V["PUB", "lag(LE, 1)"]
  t   <- function(th) (b["PUB"] - th * den) / sqrt(vb + 2 * th * cab + th^2 * va)
  unname(max(1 - pnorm(t(-m)), pnorm(t(m))))
}
eq <- bind_rows(lapply(names(samples), function(s) {
  r <- co %>% filter(sample == s, estimator == "System GMM", term == "PUB")
  data.frame(sample = s, margin = margins,
             p_tost_impact = sapply(margins, function(m) tost_p(r$estimate, r$se, m)),
             p_tost_longrun_fieller = sapply(margins, function(m)
               tost_fieller_p(gm[[s]]$b, gm[[s]]$V, m)))
})) %>%
  mutate(equiv_impact = p_tost_impact < 0.05,
         equiv_longrun_fieller = p_tost_longrun_fieller < 0.05)

cat("\n=== Public share, System GMM: long run and equivalence ===\n")
for (s in names(samples)) {
  r <- co %>% filter(sample == s, estimator == "System GMM", term == "PUB")
  e <- eq %>% filter(sample == s)
  cat(sprintf("%s: beta_LR %+.4f; Fieller 95%% %s; 90%% %s\n", s, r$beta_LR,
              fmt_fieller(r$fieller_lo, r$fieller_hi, r$fieller_type),
              fmt_fieller(r$fieller90_lo, r$fieller90_hi, r$fieller90_type)))
  cat(sprintf("   impact equivalence from margin %s; long-run (Fieller) %s\n",
              ifelse(any(e$equiv_impact),
                     sprintf("%.4f", min(e$margin[e$equiv_impact])), "none"),
              ifelse(any(e$equiv_longrun_fieller),
                     sprintf("from margin %.4f",
                             min(e$margin[e$equiv_longrun_fieller])),
                     "at no margin")))
}

# ---- 4. sample C against sample A -------------------------
# Three comparisons on the public share: the shift of the impact
# coefficient in standard errors of A (more than half an SE?), the
# equivalence tests, and the shape of the Fieller sets for the
# long-run multiplier at 95 and 90 per cent. B against A isolates the
# change of sample.
pubg <- co %>% filter(estimator == "System GMM", term == "PUB")
cmp <- function(s, ref = "A") {
  x <- pubg[pubg$sample == s, ]; r <- pubg[pubg$sample == ref, ]
  ex <- eq[eq$sample == s, ];   er <- eq[eq$sample == ref, ]
  data.frame(comparison = paste(s, "vs", ref),
             shift_in_ref_se = (x$estimate - r$estimate) / r$se,
             tost_changes = !identical(ex$equiv_impact, er$equiv_impact) ||
               !identical(ex$equiv_longrun_fieller, er$equiv_longrun_fieller),
             fieller_shape_changes = x$fieller_type != r$fieller_type,
             fieller90_shape_changes = x$fieller90_type != r$fieller90_type)
}
ch <- bind_rows(cmp("B"), cmp("C"), cmp("C", "B"))
cat("\n=== Change in the public-share results ===\n")
print(ch %>% mutate(shift_in_ref_se = round(shift_in_ref_se, 2)))
cc <- ch[ch$comparison == "C vs A", ]
cat(sprintf("\nC against A: shift %.2f SE (%s half an SE); TOST %s;\n",
            cc$shift_in_ref_se,
            ifelse(abs(cc$shift_in_ref_se) > 0.5, "more than", "within"),
            ifelse(cc$tost_changes, "changes", "unchanged")))
cat(sprintf("Fieller set at 95%% %s, at 90%% %s.\n",
            ifelse(cc$fieller_shape_changes, "changes shape", "same shape"),
            ifelse(cc$fieller90_shape_changes, "changes shape", "same shape")))
cat(sprintf("In C every Fieller set is an interval only up to the %.1f%% level.\n",
            100 * dg$level_bounded[dg$sample == "C"]))

write.csv(co, file.path(OUT, "35_hiv.csv"), row.names = FALSE)
write.csv(dg, file.path(OUT, "35_hiv_samples.csv"), row.names = FALSE)
write.csv(eq, file.path(OUT, "35_hiv_equivalence.csv"), row.names = FALSE)

# ---- 5. Supplementary Table S1, in LaTeX -------------------
cell <- function(est, tm, s) {
  r <- co %>% filter(estimator == est, term == tm, sample == s)
  if (!nrow(r)) return("---")
  sprintf("$%+.4f$ (%.4f)", r$estimate, r$se)
}
fset <- function(s) {
  r <- co %>% filter(estimator == "System GMM", term == "PUB", sample == s)
  if (r$fieller90_type == "bounded")
    sprintf("$[%+.3f,\\ %+.3f]$", r$fieller90_lo, r$fieller90_hi) else
    "whole line"
}
ROWS <- list(
  list("System GMM", "lag(LE, 1)", "Life expectancy, lagged one year"),
  list("System GMM", "PUB",        "Public share of CHE"),
  list("System GMM", "OOP_share",  "Out-of-pocket share of CHE"),
  list("System GMM", "IMM_DPT",    "Immunisation"),
  list("System GMM", "HIV_pct",    "HIV prevalence"),
  list("FE, country effects", "PUB", "Public share, FE with country effects"),
  list("FE, country and year effects", "PUB",
       "Public share, FE with country and year effects"),
  list("mixed model", "ILR1_wn", "Public balance (ILR$_1$), within, mixed model"),
  list("mixed model", "HIV_pct_wn", "HIV prevalence, within, mixed model"))
tex <- c("\\begin{tabular}{@{\\extracolsep\\fill}lccc@{}}", "\\toprule",
         "& A & B & C \\\\",
         "& Analytical panel & HIV reported & HIV reported \\\\",
         "& without HIV & without HIV & with HIV \\\\", "\\midrule")
for (r in ROWS)
  tex <- c(tex, sprintf("%s & %s & %s & %s \\\\", r[[3]],
                        cell(r[[1]], r[[2]], "A"), cell(r[[1]], r[[2]], "B"),
                        cell(r[[1]], r[[2]], "C")))
tex <- c(tex, "\\midrule",
  sprintf("Public share, long run: 90\\%% Fieller interval & %s & %s & %s \\\\",
          fset("A"), fset("B"), fset("C")),
  sprintf("Fieller sets are intervals up to & %s \\\\",
          paste(sprintf("%.1f\\%%", 100 * dg$level_bounded), collapse = " & ")),
  sprintf("Hansen $J$, $p$ & %s \\\\",
          paste(sprintf("%.2f", dg$hansen_p), collapse = " & ")),
  sprintf("AR(2), $p$ & %s \\\\",
          paste(sprintf("%.2f", dg$ar2_p), collapse = " & ")),
  sprintf("Observations & %s \\\\", paste(dg$n_obs, collapse = " & ")),
  sprintf("Countries & %s \\\\", paste(dg$n_countries, collapse = " & ")),
  "\\botrule", "\\end{tabular}")
writeLines(tex, file.path(OUT, "35_hiv_table.tex"))
cat("\n=== Supplementary table, LaTeX ===\n")
cat(tex, sep = "\n")

cat("\nSaved: 35_hiv.csv, 35_hiv_samples.csv, 35_hiv_equivalence.csv,\n")
cat("35_hiv_table.tex\n")

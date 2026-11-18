# 37: the public share with and without total health spending per head
#
# The public share of current health expenditure (CHE) says nothing
# about the amount spent: a share of 50 per cent is 500 per head when
# CHE is 1000 and 1000 per head when CHE is 2000. This script adds
# log CHE per head to the reference System GMM of Table 3, keeping
# the public and out-of-pocket shares as they are, and compares the
# coefficient on the public share on identical country-years
# (Section 6.3, Table 4):
#   A  the reference model on the analytical panel
#   B  the same model on the observations where CHE per head is
#      observed
#   C  the rows of B, with log CHE per head added
# Outcome, the other regressors, the two-way effects, the estimator
# and the observations are the same in B and C; only the regressor
# and its instruments change. In the reference model the drivers are
# strictly exogenous and instrument themselves. CHE per head is
# entered under three classifications, because spending may respond
# to mortality shocks within the year: endogenous, with collapsed
# instruments from t-2 to t-4 (the main comparison), predetermined,
# from t-1 to t-3, and strictly exogenous, like the other drivers.
# The same comparison is made for the two-way static fixed effects
# rung of Table 2, and CHE per head in current US dollars replaces
# the PPP series in a further variant.
#
# SH.XPD.CHEX.PP.CD  current health expenditure per head, PPP,
#                    current international dollars (main measure)
# SH.XPD.CHEX.PC.CD  current health expenditure per head, current US$
# Both come from the WHO Global Health Expenditure Database through
# the World Development Indicators. With year effects in every
# model, the common movement of the international-dollar price level
# is absorbed, so the within-country variation of the PPP series is
# close to a change in real terms; the US-dollar series also moves
# with exchange rates.
#
# Reads  data/wdi_total_spending.rds (downloaded through WDI if
#        missing), and data/ through build_panel()
# Writes output/37_total_spending.csv, 37_total_spending_diag.csv,
#        37_total_spending_samples.csv, 37_total_spending_rows.csv
#        (the country-years of each sample), 37_total_spending_table.tex

source("00_common.R")
suppressPackageStartupMessages(library(plm))

CACHE <- file.path(DATA, "wdi_total_spending.rds")

if (file.exists(CACHE)) {
  w <- readRDS(CACHE); cat("using cached WDI extract\n")
} else {
  cat("downloading from World Bank Open Data ...\n")
  suppressPackageStartupMessages(library(WDI))
  w <- WDI(country = "all",
           indicator = c(CHE_PPP = "SH.XPD.CHEX.PP.CD",
                         CHE_USD = "SH.XPD.CHEX.PC.CD"),
           start = 2000, end = 2024, extra = FALSE)
  saveRDS(w, CACHE)
}
w <- w %>% dplyr::select(iso3c, year, CHE_PPP, CHE_USD) %>%
  filter(!is.na(iso3c))

# coverage of the extract, before the merge with the panel
cov <- w %>% filter(!is.na(CHE_PPP))
cat(sprintf("extract coverage, CHE per head PPP: %d codes (countries and World Bank aggregates), years %d-%d, %d code-years\n",
            n_distinct(cov$iso3c), min(cov$year), max(cov$year), nrow(cov)))

# ---- samples -----------------------------------------------
datA <- build_panel(zero_rule = "multiplicative") %>%
  left_join(w, by = c("iso3c", "year"))
stopifnot(nrow(datA) == nrow(build_panel(zero_rule = "multiplicative")))

# the US-dollar series of the extract is the one already in the panel
chk <- datA %>% filter(!is.na(CHE_USD), !is.na(HEALTH_PC))
cat(sprintf("CHE per head in US$, extract against panel: max abs difference %.3g on %d obs\n",
            max(abs(chk$CHE_USD - chk$HEALTH_PC)), nrow(chk)))

datA <- datA %>%
  mutate(log_CHE     = ifelse(CHE_PPP > 0, log(CHE_PPP), NA_real_),
         log_CHE_USD = ifelse(CHE_USD > 0, log(CHE_USD), NA_real_))

cat(sprintf("\nanalytical panel: %d obs, %d countries\n",
            nrow(datA), n_distinct(datA$iso3c)))
cat(sprintf("CHE per head, PPP: missing %d, zero or negative %d\n",
            sum(is.na(datA$CHE_PPP)), sum(datA$CHE_PPP <= 0, na.rm = TRUE)))
cat(sprintf("CHE per head, US$: missing %d, zero or negative %d\n",
            sum(is.na(datA$CHE_USD)), sum(datA$CHE_USD <= 0, na.rm = TRUE)))

# B and C use the same rows: those with log CHE (PPP) observed. No
# value is replaced or imputed to take the logarithm.
datB <- datA %>% filter(!is.na(log_CHE)) %>% arrange(iso3c, year)
cat(sprintf("common sample of B and C: %d obs, %d countries (%d obs dropped)\n",
            nrow(datB), n_distinct(datB$iso3c), nrow(datA) - nrow(datB)))
# the US-dollar variant is fitted on the rows of B where it is observed
datU <- datB %>% filter(!is.na(log_CHE_USD))
cat(sprintf("US$ variant: %d obs, %d countries\n",
            nrow(datU), n_distinct(datU$iso3c)))

pdA <- pdata.frame(datA, index = c("iso3c", "year"))
pdB <- pdata.frame(datB, index = c("iso3c", "year"))
pdU <- pdata.frame(datU, index = c("iso3c", "year"))

# ---- the scale variable ------------------------------------
dec <- function(x, g) {
  mu <- tapply(x, g, mean); cm <- as.numeric(mu[as.character(g)])
  c(mean = mean(x), sd_between = sd(mu), sd_within = sd(x - cm))
}
cat("\n=== log CHE per head (PPP), within and between ===\n")
print(round(dec(datB$log_CHE, datB$iso3c), 3))
wc <- datB %>% group_by(iso3c) %>%
  mutate(a = log_CHE - mean(log_CHE), b = PUB - mean(PUB),
         g = log_GNI - mean(log_GNI)) %>% ungroup()
cat(sprintf("within-country correlation with the public share: %.3f\n",
            cor(wc$a, wc$b)))
cat(sprintf("within-country correlation with log income      : %.3f\n",
            cor(wc$a, wc$g)))
# residual within variation of log CHE after the two-way effects and
# the other regressors, the variation that identifies its coefficient
r_che <- resid(plm(log_CHE ~ log_GNI + EDU_GROSS + PUB + OOP_share +
                     IMM_DPT + SANITATION, data = pdB,
                   model = "within", effect = "twoways"))
cat(sprintf("residual within SD of log CHE (two-way, given the drivers): %.3f\n",
            sd(as.numeric(r_che))))

# ---- System GMM fits ---------------------------------------
# The formulas carry the lag numbers literally, because pgmm
# evaluates them in its own environment.
f_ref <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION | lag(LE, 2:5)
f_end <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION + log_CHE |
         lag(LE, 2:5) + lag(log_CHE, 2:4)
f_pre <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION + log_CHE |
         lag(LE, 2:5) + lag(log_CHE, 1:3)
f_exo <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION + log_CHE | lag(LE, 2:5)
f_usd <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION + log_CHE_USD |
         lag(LE, 2:5) + lag(log_CHE_USD, 2:4)

gmm <- function(f, pd, tr = "ld")
  pgmm(f, data = pd, effect = "twoways", model = "twosteps",
       transformation = tr, collapse = TRUE)

# TOST against +/- margin, as in script 19
tost_p <- function(b, se, m)
  max(1 - pnorm((b + m) / se), pnorm((b - m) / se))

# one row per regressor, with the diagnostics of the fit
gmm_rows <- function(f, pd, label, che_class) {
  m  <- gmm(f, pd)
  md <- gmm(f, pd, tr = "d")           # difference twin
  b  <- coef(m); V <- pgmm_vcov(m)
  sg <- sargan(m); sgd <- sargan(md)
  dih_C  <- as.numeric(sg$statistic - sgd$statistic)
  dih_df <- as.numeric(sg$parameter - sgd$parameter)
  ar <- sapply(1:3, function(k)
    tryCatch(as.numeric(mtest(m, order = k)$p.value),
             error = function(e) NA_real_))
  lr95 <- longrun(b, V, "lag(LE, 1)", level = 0.95)
  lr90 <- longrun(b, V, "lag(LE, 1)", level = 0.90)
  fl   <- fieller_levels(b, V, "lag(LE, 1)")
  a <- unname(b["lag(LE, 1)"]); se_a <- sqrt(V["lag(LE, 1)", "lag(LE, 1)"])
  terms <- setdiff(names(b), "lag(LE, 1)")
  out <- data.frame(
    spec = label, estimator = "System GMM", che_class = che_class,
    term = terms, beta = unname(b[terms]),
    se = sqrt(diag(V)[terms]),
    alpha = a, se_alpha = se_a,
    n_obs = nrow(pd), n_countries = n_distinct(pd$iso3c),
    n_inst = as.numeric(sg$parameter) + length(b),
    hansen = as.numeric(sg$statistic), hansen_df = as.numeric(sg$parameter),
    hansen_p = as.numeric(sg$p.value),
    diffhansen = dih_C, diffhansen_df = dih_df,
    diffhansen_p = pchisq(dih_C, dih_df, lower.tail = FALSE),
    ar1_p = ar[1], ar2_p = ar[2], ar3_p = ar[3],
    row.names = NULL)
  out <- out %>%
    left_join(lr95 %>% dplyr::select(term, beta_LR,
                                     fieller95_lo = fieller_lo,
                                     fieller95_hi = fieller_hi,
                                     fieller95_type = fieller_type),
              by = "term") %>%
    left_join(lr90 %>% dplyr::select(term, fieller90_lo = fieller_lo,
                                     fieller90_hi = fieller_hi,
                                     fieller90_type = fieller_type),
              by = "term") %>%
    left_join(fl, by = "term")
  out
}

res_gmm <- bind_rows(
  gmm_rows(f_ref, pdA, "A reference model, analytical panel", "not included"),
  gmm_rows(f_ref, pdB, "B reference model, common sample", "not included"),
  gmm_rows(f_end, pdB, "C log CHE per head (PPP), endogenous", "endogenous, t-2 to t-4"),
  gmm_rows(f_pre, pdB, "C log CHE per head (PPP), predetermined", "predetermined, t-1 to t-3"),
  gmm_rows(f_exo, pdB, "C log CHE per head (PPP), strictly exogenous", "strictly exogenous"),
  gmm_rows(f_usd, pdU, "C log CHE per head (US$), endogenous", "endogenous, t-2 to t-4")
)

# ---- static two-way fixed effects --------------------------
fe_rows <- function(f, pd, label) {
  m <- plm(f, data = pd, model = "within", effect = "twoways")
  V <- plm::vcovHC(m, cluster = "group")
  b <- coef(m)
  data.frame(spec = label, estimator = "Two-way fixed effects",
             che_class = NA_character_, term = names(b),
             beta = unname(b), se = sqrt(diag(V)[names(b)]),
             n_obs = nrow(pd), n_countries = n_distinct(pd$iso3c),
             row.names = NULL)
}
s_ref <- LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION
s_che <- LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION +
              log_CHE
res_fe <- bind_rows(
  fe_rows(s_ref, pdA, "A reference model, analytical panel"),
  fe_rows(s_ref, pdB, "B reference model, common sample"),
  fe_rows(s_che, pdB, "C log CHE per head (PPP)")
)

res <- bind_rows(res_gmm, res_fe) %>%
  mutate(z = beta / se, p = 2 * pnorm(-abs(z)),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         mde80 = (qnorm(0.975) + qnorm(0.80)) * se)
res$tost_p_0125 <- ifelse(res$term == "PUB",
                          mapply(tost_p, res$beta, res$se, 0.0125), NA_real_)
res$tost_p_025  <- ifelse(res$term == "PUB",
                          mapply(tost_p, res$beta, res$se, 0.025), NA_real_)

# B and C must use exactly the same country-years
stopifnot(all(res$n_obs[grepl("^[BC] .*PPP|^B ", res$spec)] == nrow(datB)))

# ---- console report ----------------------------------------
cat("\n=== Public share of CHE: A, B and C ===\n")
print(res %>% filter(term == "PUB") %>%
        mutate(across(c(beta, se, lo, hi, mde80), ~ round(.x, 4)),
               across(c(p, tost_p_0125, tost_p_025), ~ round(.x, 3)),
               alpha = round(alpha, 4)) %>%
        dplyr::select(estimator, spec, beta, se, lo, hi, p, mde80,
                      tost_p_0125, tost_p_025, alpha, n_obs, n_countries))

cat("\n=== log CHE per head ===\n")
print(res %>% filter(term %in% c("log_CHE", "log_CHE_USD")) %>%
        mutate(across(c(beta, se, lo, hi), ~ round(.x, 4)),
               p = round(p, 3)) %>%
        dplyr::select(estimator, spec, beta, se, lo, hi, p))
b_c <- res %>% filter(spec == "C log CHE per head (PPP), endogenous",
                      term == "log_CHE")
cat(sprintf("a 10%% rise in CHE per head, endogenous classification: %+.4f years within the year\n",
            log(1.1) * b_c$beta))

cat("\n=== System GMM diagnostics ===\n")
diag <- res_gmm %>% group_by(spec) %>% slice(1) %>% ungroup() %>%
  dplyr::select(spec, che_class, n_obs, n_countries, n_inst, alpha,
                se_alpha, ar1_p, ar2_p, ar3_p, hansen, hansen_df,
                hansen_p, diffhansen, diffhansen_df, diffhansen_p,
                level_bounded) %>%
  arrange(spec)
print(diag %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
      width = 200)

cat("\n=== Long-run multiplier of the public share ===\n")
print(res_gmm %>% filter(term == "PUB") %>%
        mutate(f90 = mapply(fmt_fieller, fieller90_lo, fieller90_hi,
                            fieller90_type, 3),
               f95 = mapply(fmt_fieller, fieller95_lo, fieller95_hi,
                            fieller95_type, 3),
               beta_LR = round(beta_LR, 4),
               level_bounded = round(level_bounded, 3)) %>%
        dplyr::select(spec, beta_LR, f90, f95, level_bounded))

cat("\n=== Change from B to C in the public-share coefficient ===\n")
chg <- function(est, lab_c) {
  rb <- res %>% filter(estimator == est, term == "PUB",
                       spec == "B reference model, common sample")
  rc <- res %>% filter(estimator == est, term == "PUB", spec == lab_c)
  cat(sprintf("%-22s %-48s %+.4f -> %+.4f, change %+.4f (%.2f SE of B); SE ratio C/B %.2f\n",
              est, lab_c, rb$beta, rc$beta, rc$beta - rb$beta,
              (rc$beta - rb$beta) / rb$se, rc$se / rb$se))
}
for (lab in grep("^C .*PPP", unique(res_gmm$spec), value = TRUE))
  chg("System GMM", lab)
chg("Two-way fixed effects", "C log CHE per head (PPP)")

# ---- files -------------------------------------------------
write.csv(res, file.path(OUT, "37_total_spending.csv"), row.names = FALSE)
write.csv(diag, file.path(OUT, "37_total_spending_diag.csv"), row.names = FALSE)
smp <- data.frame(sample = c("A analytical panel", "B and C common sample",
                             "US$ variant"),
                  n_obs = c(nrow(datA), nrow(datB), nrow(datU)),
                  n_countries = c(n_distinct(datA$iso3c),
                                  n_distinct(datB$iso3c),
                                  n_distinct(datU$iso3c)))
write.csv(smp, file.path(OUT, "37_total_spending_samples.csv"), row.names = FALSE)

# the country-years behind each sample
key <- function(d) paste(d$iso3c, d$year)
rows <- datA %>%
  dplyr::select(iso3c, country, year) %>%
  mutate(in_A = TRUE,
         in_B_and_C = key(datA) %in% key(datB),
         in_USD_variant = key(datA) %in% key(datU))
write.csv(rows, file.path(OUT, "37_total_spending_rows.csv"), row.names = FALSE)
cat(sprintf("country-years written: %d (A), %d (B and C), %d (US$ variant)\n",
            sum(rows$in_A), sum(rows$in_B_and_C), sum(rows$in_USD_variant)))

# LaTeX body of Table 4. Columns: (1) the reference model, which on
# the common sample is model B; (2) to (4) model C with log CHE per
# head endogenous, predetermined and strictly exogenous.
sA <- "A reference model, analytical panel"
sB <- "B reference model, common sample"
gcols <- c(sB, "C log CHE per head (PPP), endogenous",
           "C log CHE per head (PPP), predetermined",
           "C log CHE per head (PPP), strictly exogenous")
fcols <- c(sB, "C log CHE per head (PPP)", NA, NA)
get <- function(est, sp, tn) {
  if (is.na(sp)) return(NULL)
  r <- res %>% filter(estimator == est, spec == sp, term == tn)
  if (nrow(r) == 0) NULL else r
}
bse <- function(r, digits = 4) {
  if (is.null(r)) return(c("", ""))
  c(sprintf("$%+.*f$", digits, r$beta), sprintf("(%.*f)", digits, r$se))
}
row2 <- function(label, est, cols, tn, digits = 4) {
  cells <- lapply(cols, function(sp) bse(get(est, sp, tn), digits))
  c(paste(label, paste(sapply(cells, `[`, 1), collapse = " & "), sep = " & "),
    paste("", paste(sapply(cells, `[`, 2), collapse = " & "), sep = " & "))
}
rowf <- function(label, f) paste(label, paste(sapply(gcols, f), collapse = " & "), sep = " & ")
dgv <- function(sp, col) diag[[col]][diag$spec == sp]
tex <- c(
  "\\multicolumn{5}{@{}l}{\\textit{System GMM}} \\\\",
  paste0(row2("Public share of CHE", "System GMM", gcols, "PUB"), " \\\\"),
  paste0(rowf("\\quad 95\\% CI", function(sp) {
    r <- get("System GMM", sp, "PUB"); sprintf("$[%+.4f, %+.4f]$", r$lo, r$hi)}), " \\\\"),
  paste0(rowf("\\quad Equivalence $p$", function(sp) {
    r <- get("System GMM", sp, "PUB"); sprintf("%.3f", r$tost_p_0125)}), " \\\\"),
  paste0(row2("Log CHE per head", "System GMM", gcols, "log_CHE", 3), " \\\\"),
  paste0(rowf("Persistence $\\hat{\\alpha}$", function(sp)
    sprintf("%.3f (%.3f)", dgv(sp, "alpha"), dgv(sp, "se_alpha"))), " \\\\"),
  paste0(rowf("Instruments", function(sp) sprintf("%d", as.integer(dgv(sp, "n_inst")))), " \\\\"),
  paste0(rowf("AR(2) $p$", function(sp) sprintf("%.2f", dgv(sp, "ar2_p"))), " \\\\"),
  paste0(rowf("Hansen $J$ $p$", function(sp) sprintf("%.2f", dgv(sp, "hansen_p"))), " \\\\"),
  paste0(rowf("Difference-in-Hansen $p$", function(sp) sprintf("%.2f", dgv(sp, "diffhansen_p"))), " \\\\"),
  "\\midrule",
  "\\multicolumn{5}{@{}l}{\\textit{Static fixed effects, country and year effects}} \\\\",
  paste0(row2("Public share of CHE", "Two-way fixed effects", fcols, "PUB"), " \\\\"),
  paste0(row2("Log CHE per head", "Two-way fixed effects", fcols, "log_CHE", 3), " \\\\"),
  "\\midrule",
  sprintf("Observations & %s \\\\", paste(rep(nrow(datB), 4), collapse = " & ")),
  sprintf("Countries & %s \\\\", paste(rep(n_distinct(datB$iso3c), 4), collapse = " & "))
)
# model B is model A when no observation is lost
if (nrow(datB) == nrow(datA))
  stopifnot(isTRUE(all.equal(res$beta[res$spec == sA], res$beta[res$spec == sB])))
writeLines(tex, file.path(OUT, "37_total_spending_table.tex"))

cat("\nSaved: 37_total_spending.csv, 37_total_spending_diag.csv,\n")
cat("       37_total_spending_samples.csv, 37_total_spending_rows.csv,\n")
cat("       37_total_spending_table.tex\n")

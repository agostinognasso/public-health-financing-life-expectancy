# 21: exogeneity classification, spending per head, country trends
#     and longer distributed lags
#
# A regressor left out of the GMM instrument block of pgmm is treated
# as strictly exogenous and instruments itself. The System GMM is
# re-estimated with the drivers (A) strictly exogenous,
# (B) predetermined, instrumented from t-1, (C) endogenous,
# instrumented from t-2, and (D) with only the financing block
# endogenous (Table 2, Section 5.2). The script also replaces the
# share with log public spending per head (E, F), adds
# country-specific linear trends to the static model (G), and sums
# the distributed lag over 3, 5 and 8 lags (H).
#
# Reads  data/ through build_panel()
# Writes output/21_exogeneity.csv, 21_dlag_cumulative.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

dat <- build_panel(zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c", "year"))

sysfit <- function(f, tr = "ld")
  pgmm(f, data = pd, effect = "twoways", model = "twosteps",
       transformation = tr, collapse = TRUE)

grab <- function(m, label) {
  if (is.null(m)) return(NULL)
  b <- coef(m); V <- pgmm_vcov(m)
  keep <- intersect(names(b),
                    c("lag(LE, 1)", "log_GNI", "EDU_GROSS", "PUB",
                      "OOP_share", "IMM_DPT", "SANITATION", "log_PUB_PC"))
  a <- b["lag(LE, 1)"]
  sg <- tryCatch(sargan(m), error = function(e) NULL)
  ar2 <- tryCatch(mtest(m, order = 2), error = function(e) NULL)
  data.frame(spec = label, term = keep, beta = b[keep],
             se = sqrt(diag(V)[keep]),
             alpha = as.numeric(a),
             beta_LR = b[keep] / (1 - as.numeric(a)),
             n_inst = dim(m$W[[1]])[2],
             hansen_p = ifelse(is.null(sg), NA, sg$p.value),
             ar2_p = ifelse(is.null(ar2), NA, ar2$p.value),
             row.names = NULL) %>%
    mutate(p = 2 * pnorm(-abs(beta / se)),
           lo = beta - 1.96 * se, hi = beta + 1.96 * se)
}

# ---- A. drivers strictly exogenous -------------------------
fA <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION | lag(LE, 2:5)

# ---- B. drivers predetermined: instruments from lag 1 -------
fB <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION |
      lag(LE, 2:5) + lag(PUB, 1:3) + lag(OOP_share, 1:3) +
      lag(log_GNI, 1:3) + lag(EDU_GROSS, 1:3) + lag(IMM_DPT, 1:3) +
      lag(SANITATION, 1:3)

# ---- C. drivers endogenous: instruments from lag 2 ----------
fC <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION |
      lag(LE, 2:5) + lag(PUB, 2:4) + lag(OOP_share, 2:4) +
      lag(log_GNI, 2:4) + lag(EDU_GROSS, 2:4) + lag(IMM_DPT, 2:4) +
      lag(SANITATION, 2:4)

# ---- D. only the financing block endogenous ----------------
fD <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION |
      lag(LE, 2:5) + lag(PUB, 2:4) + lag(OOP_share, 2:4)

res <- bind_rows(
  grab(tryCatch(sysfit(fA), error = function(e) NULL),
       "A drivers strictly exogenous"),
  grab(tryCatch(sysfit(fB), error = function(e) NULL),
       "B drivers predetermined (instruments from t-1)"),
  grab(tryCatch(sysfit(fC), error = function(e) NULL),
       "C drivers endogenous (instruments from t-2)"),
  grab(tryCatch(sysfit(fD), error = function(e) NULL),
       "D financing block endogenous only")
)

cat("=== Effect of the exogeneity classification on beta(PUB) ===\n")
print(res %>% filter(term == "PUB") %>%
        mutate(across(c(beta, se, lo, hi, beta_LR, alpha), ~ round(.x, 4)),
               across(c(p, hansen_p, ar2_p), ~ round(.x, 3))) %>%
        dplyr::select(spec, beta, se, p, lo, hi, alpha, beta_LR,
                      n_inst, hansen_p, ar2_p))

cat("\n=== Immunisation under the same classifications ===\n")
print(res %>% filter(term == "IMM_DPT") %>%
        mutate(across(c(beta, se, lo, hi, beta_LR), ~ round(.x, 4)),
               p = round(p, 4)) %>%
        dplyr::select(spec, beta, se, p, lo, hi, beta_LR))

cat("\n=== Full coefficient table, all classifications ===\n")
print(res %>% mutate(across(c(beta, se, lo, hi, beta_LR), ~ round(.x, 4)),
                     p = round(p, 3)) %>%
        dplyr::select(spec, term, beta, se, p, lo, hi))

# ------------------------------------------------------------
# E. public health spending per head, instead of its share
#    of current health expenditure
# ------------------------------------------------------------
cat("\n\n=== Public spending as a level, not a share ===\n")
cat("The share rises when public spending grows and also when private\n")
cat("spending falls; log public spending per head moves only with the first.\n\n")
cat(sprintf("log_PUB_PC available for %d of %d observations (%.1f%%)\n",
            sum(!is.na(dat$log_PUB_PC)), nrow(dat),
            100 * mean(!is.na(dat$log_PUB_PC))))
cat(sprintf("correlation with the public share, within country: %.3f\n",
            {w <- dat %>% filter(!is.na(log_PUB_PC)) %>% group_by(iso3c) %>%
                   mutate(a = log_PUB_PC - mean(log_PUB_PC),
                          b = z_PUB - mean(z_PUB)) %>% ungroup()
             cor(w$a, w$b)}))

fE <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + log_PUB_PC + OOP_share +
           IMM_DPT + SANITATION | lag(LE, 2:5)
fE2 <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + log_PUB_PC + OOP_share +
            IMM_DPT + SANITATION |
       lag(LE, 2:5) + lag(log_PUB_PC, 2:4) + lag(OOP_share, 2:4)
resE <- bind_rows(
  grab(tryCatch(sysfit(fE),  error = function(e) NULL), "E log public $/head, exogenous"),
  grab(tryCatch(sysfit(fE2), error = function(e) NULL), "F log public $/head, endogenous")
)
print(resE %>% filter(term == "log_PUB_PC") %>%
        mutate(across(c(beta, se, lo, hi, beta_LR, alpha), ~ round(.x, 4)),
               p = round(p, 4)) %>%
        dplyr::select(spec, beta, se, p, lo, hi, alpha, beta_LR))
cat("\nThe coefficient is in years of life expectancy per log point of\n")
cat("public health spending per head; a 10 per cent rise is log(1.1)\n")
cat("times the coefficient.\n")
be <- resE$beta[resE$term == "log_PUB_PC" & grepl("^E", resE$spec)]
cat(sprintf("A 10%% rise in public spending per head: %+.4f years (short run)\n",
            log(1.1) * be))

# static FE for the level variable, one-way and two-way
m1 <- plm(LE ~ log_GNI + EDU_GROSS + log_PUB_PC + OOP_share + IMM_DPT +
               SANITATION, data = pd, model = "within")
m2 <- plm(LE ~ log_GNI + EDU_GROSS + log_PUB_PC + OOP_share + IMM_DPT +
               SANITATION, data = pd, model = "within", effect = "twoways")
cat(sprintf("\nstatic FE one-way  : %+.4f (SE %.4f)\n", coef(m1)["log_PUB_PC"],
            sqrt(plm::vcovHC(m1, cluster = "group")["log_PUB_PC", "log_PUB_PC"])))
cat(sprintf("static FE two-way  : %+.4f (SE %.4f)\n", coef(m2)["log_PUB_PC"],
            sqrt(plm::vcovHC(m2, cluster = "group")["log_PUB_PC", "log_PUB_PC"])))
cat("HEALTH_PC is in current US dollars, so without year effects the\n")
cat("estimate also picks up common price changes. Script 22 repeats the\n")
cat("ladder on public spending as a share of GDP.\n")

# ------------------------------------------------------------
# G. country-specific linear trends
# ------------------------------------------------------------
cat("\n\n=== Country-specific linear trends ===\n")
cat("The level-equation instruments of System GMM assume mean\n")
cat("stationarity; country-specific linear trends are a check on it.\n")
dat_tr <- dat %>% group_by(iso3c) %>% mutate(t = year - min(year)) %>% ungroup()
pdt <- pdata.frame(dat_tr, index = c("iso3c", "year"))
mt1 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION +
                factor(iso3c):t, data = pdt, model = "within")
cat(sprintf("static FE with country trends: beta(PUB) = %+.4f\n",
            coef(mt1)["PUB"]))
fe1 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION,
           data = pd, model = "within")
fe2 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION,
           data = pd, model = "within", effect = "twoways")
cat(sprintf("(without trends: one-way FE %+.4f, two-way FE %+.4f)\n",
            coef(fe1)["PUB"], coef(fe2)["PUB"]))

# ------------------------------------------------------------
# H. deeper distributed lag
# ------------------------------------------------------------
cat("\n\n=== Distributed lag at three, five and eight lags ===\n")
DL_X <- c("PUB", "log_GNI", "IMM_DPT", "SANITATION", "EDU_GROSS", "OOP_share")
dl <- list()
for (K in c(3, 5, 8)) {
  m <- tryCatch(plm(as.formula(sprintf(
        "LE ~ lag(PUB, 0:%d) + lag(log_GNI, 0:%d) + lag(IMM_DPT, 0:%d) +
              lag(SANITATION, 0:%d) + lag(EDU_GROSS, 0:%d) + lag(OOP_share, 0:%d)",
        K, K, K, K, K, K)), data = pd, model = "within", effect = "twoways"),
      error = function(e) NULL)
  if (is.null(m)) { cat(sprintf("lags 0-%d: failed\n", K)); next }
  V <- plm::vcovHC(m, cluster = "group", method = "arellano")
  idx <- grep("^lag\\(PUB", names(coef(m)))
  R <- rep(0, length(coef(m))); R[idx] <- 1
  cum <- sum(coef(m)[idx]); se <- sqrt(as.numeric(t(R) %*% V %*% R))
  cat(sprintf("lags 0-%d (%d coefficients): cumulative %+.4f (SE %.4f, p = %.3f), 95%% CI [%+.4f, %+.4f]\n",
              K, length(idx), cum, se, 2 * pnorm(-abs(cum / se)),
              cum - 1.96 * se, cum + 1.96 * se))

  # the same cumulative sum for every driver, from the same fit
  for (v in DL_X) {
    iv <- grep(paste0("^lag\\(", v, ","), names(coef(m)))
    Rv <- rep(0, length(coef(m))); Rv[iv] <- 1
    cv <- sum(coef(m)[iv]); sv <- sqrt(as.numeric(t(Rv) %*% V %*% Rv))
    dl[[paste(K, v)]] <- data.frame(lags = sprintf("0-%d", K), driver = v,
                                    n_coef = length(iv), cumulative = cv,
                                    se = sv, p = 2 * pnorm(-abs(cv / sv)),
                                    lo = cv - 1.96 * sv, hi = cv + 1.96 * sv,
                                    n_obs = nobs(m),
                                    n_countries = pdim(m)$nT$n)
  }
}
dl <- bind_rows(dl)
cat("\nCumulative effect of every driver (two-way FE, SE clustered by country):\n")
print(dl %>% mutate(across(c(cumulative, se, lo, hi), ~ round(.x, 4)),
                    p = round(p, 3)) %>%
        dplyr::select(lags, driver, cumulative, se, p, lo, hi, n_obs),
      row.names = FALSE)
cat("\nThe individual lag coefficients are nearly collinear, so only their\n")
cat("sums are reported.\n")

write.csv(bind_rows(res, resE), file.path(OUT, "21_exogeneity.csv"),
          row.names = FALSE)
write.csv(dl, file.path(OUT, "21_dlag_cumulative.csv"), row.names = FALSE)
cat("\nSaved: 21_exogeneity.csv, 21_dlag_cumulative.csv\n")

# 22: public spending as a level rather than a share
#
# SH.XPD.GHED.CH.ZS is a share of current health expenditure: it
# rises when public spending grows and also when private spending
# falls. The ladder is repeated on two measures whose denominator
# lies outside the health accounts (Section 6.1):
#   SH.XPD.GHED.GD.ZS  public health spending, % of GDP
#   SH.XPD.GHED.PC.CD  public health spending per head, current US$
#
# Reads  data/wdi_level_spending.rds (downloaded through WDI if
#        missing), and data/ through build_panel()
# Writes output/22_level_ladder.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(WDI))

CACHE <- file.path(DATA, "wdi_level_spending.rds")

if (file.exists(CACHE)) {
  w <- readRDS(CACHE); cat("using cached WDI extract\n")
} else {
  cat("downloading from World Bank Open Data ...\n")
  w <- WDI(country = "all",
           indicator = c(PUB_GDP = "SH.XPD.GHED.GD.ZS",
                         PUB_PCAP = "SH.XPD.GHED.PC.CD",
                         CHE_GDP = "SH.XPD.CHEX.GD.ZS"),
           start = 2000, end = 2024, extra = FALSE)
  saveRDS(w, CACHE)
}
w <- w %>% dplyr::select(iso3c, year, PUB_GDP, PUB_PCAP, CHE_GDP) %>%
  filter(!is.na(iso3c))

dat <- build_panel(zero_rule = "multiplicative") %>%
  left_join(w, by = c("iso3c", "year")) %>%
  mutate(log_PUB_PCAP = ifelse(PUB_PCAP > 0, log(PUB_PCAP), NA_real_))

cat(sprintf("\nmerged: PUB_GDP present for %d of %d obs (%.1f%%)\n",
            sum(!is.na(dat$PUB_GDP)), nrow(dat), 100*mean(!is.na(dat$PUB_GDP))))

cat("\n=== The share and the level are different regressors ===\n")
wc <- dat %>% filter(!is.na(PUB_GDP)) %>% group_by(iso3c) %>%
  mutate(sh = z_PUB - mean(z_PUB), gd = PUB_GDP - mean(PUB_GDP),
         pc = log_PUB_PCAP - mean(log_PUB_PCAP, na.rm = TRUE)) %>% ungroup()
cat(sprintf("within-country correlation, public share of CHE vs %% of GDP : %.3f\n",
            cor(wc$sh, wc$gd, use = "complete.obs")))
cat(sprintf("within-country correlation, public share of CHE vs log per head: %.3f\n",
            cor(wc$sh, wc$pc, use = "complete.obs")))

dec <- function(v) {
  x <- dat[[v]]; g <- dat$iso3c; ok <- !is.na(x); x <- x[ok]; g <- g[ok]
  mu <- tapply(x, g, mean); cm <- as.numeric(mu[as.character(g)])
  c(mean = mean(x), sd_between = sd(mu), sd_within = sd(x - cm))
}
cat("\nwithin/between SD:\n")
print(round(rbind(`public share of CHE` = dec("z_PUB"),
                  `public % of GDP`     = dec("PUB_GDP"),
                  `log public per head` = dec("log_PUB_PCAP")), 3))

# ---- the ladder on each measure ----------------------------
run_ladder <- function(focal, label) {
  X <- c("log_GNI", "EDU_GROSS", focal, "OOP_share", "IMM_DPT", "SANITATION")
  d <- dat %>% filter(!is.na(.data[[focal]]))
  pd <- pdata.frame(d, index = c("iso3c", "year"))
  rhs <- paste(X, collapse = " + ")

  m1 <- plm(as.formula(paste("LE ~", rhs)), data = pd, model = "within")
  m2 <- plm(as.formula(paste("LE ~", rhs)), data = pd, model = "within",
            effect = "twoways")
  m3 <- plm(as.formula(paste("LE ~ lag(LE,1) +", rhs)), data = pd,
            model = "within", effect = "twoways")
  fg <- as.formula(paste("LE ~ lag(LE,1) +", rhs, "| lag(LE, 2:5)"))
  fg2 <- as.formula(paste("LE ~ lag(LE,1) +", rhs,
                          "| lag(LE, 2:5) + lag(", focal, ", 2:4)"))
  m4 <- tryCatch(pgmm(fg, data = pd, effect = "twoways", model = "twosteps",
                      transformation = "ld", collapse = TRUE),
                 error = function(e) NULL)
  m5 <- tryCatch(pgmm(fg2, data = pd, effect = "twoways", model = "twosteps",
                      transformation = "ld", collapse = TRUE),
                 error = function(e) NULL)

  row <- function(m, rung, gmm = FALSE) {
    if (is.null(m)) return(NULL)
    b <- coef(m)
    V <- if (gmm) pgmm_vcov(m) else plm::vcovHC(m, cluster = "group",
                                                method = "arellano")
    a <- if ("lag(LE, 1)" %in% names(b)) b[["lag(LE, 1)"]] else NA_real_
    data.frame(measure = label, rung = rung, beta = b[[focal]],
               se = sqrt(V[focal, focal]), alpha = a,
               beta_LR = ifelse(is.na(a), b[[focal]], b[[focal]] / (1 - a)),
               n = nrow(d), row.names = NULL)
  }
  bind_rows(row(m1, "1 FE one-way"), row(m2, "2 FE two-way"),
            row(m3, "3 LSDV two-way"),
            row(m4, "4 System GMM, exogenous", TRUE),
            row(m5, "5 System GMM, endogenous", TRUE)) %>%
    mutate(p = 2 * pnorm(-abs(beta / se)),
           lo = beta - 1.96 * se, hi = beta + 1.96 * se)
}

res <- bind_rows(
  run_ladder("z_PUB",        "public share of CHE"),
  run_ladder("PUB_GDP",      "public health spending, % of GDP"),
  run_ladder("log_PUB_PCAP", "log public health spending per head")
)

cat("\n=== The ladder on three measures of public spending ===\n")
print(res %>% mutate(across(c(beta, se, alpha, beta_LR, lo, hi), ~ round(.x, 4)),
                     p = round(p, 3)) %>%
        dplyr::select(measure, rung, beta, se, p, lo, hi, alpha, beta_LR))

cat("\n=== Two-way FE on the level measures, in policy units ===\n")
g <- function(meas, rung, col = "beta")
  res[[col]][res$measure == meas & res$rung == rung]
b_gdp_fe <- g("public health spending, % of GDP", "2 FE two-way")
s_gdp_fe <- g("public health spending, % of GDP", "2 FE two-way", "se")
cat(sprintf("public %% of GDP, two-way FE: %+.4f years per pp of GDP (SE %.4f)\n",
            b_gdp_fe, s_gdp_fe))
cat(sprintf("  1 pp of GDP: %+.3f years [%+.3f, %+.3f]\n",
            b_gdp_fe, b_gdp_fe - 1.96 * s_gdp_fe, b_gdp_fe + 1.96 * s_gdp_fe))
b_pc_fe <- g("log public health spending per head", "2 FE two-way")
s_pc_fe <- g("log public health spending per head", "2 FE two-way", "se")
cat(sprintf("log per head, two-way FE: %+.4f years per log point (SE %.4f)\n",
            b_pc_fe, s_pc_fe))
# 0.0953 = log(1.1)
cat(sprintf("  a 10%% rise: %+.4f years [%+.4f, %+.4f]\n",
            0.0953 * b_pc_fe, 0.0953 * (b_pc_fe - 1.96 * s_pc_fe),
            0.0953 * (b_pc_fe + 1.96 * s_pc_fe)))

cat("\n=== Sign of the coefficient across measures ===\n")
for (meas in unique(res$measure)) {
  sub <- res %>% filter(measure == meas)
  cat(sprintf("%-38s one-way %+.4f | two-way %+.4f | SysGMM %+.4f\n", meas,
              sub$beta[sub$rung == "1 FE one-way"],
              sub$beta[sub$rung == "2 FE two-way"],
              ifelse(any(sub$rung == "4 System GMM, exogenous"),
                     sub$beta[sub$rung == "4 System GMM, exogenous"], NA)))
}
cat("\nrungs with p < 0.05:\n")
print(res %>% filter(p < 0.05) %>% dplyr::select(measure, rung, beta, se, p) %>%
        mutate(across(c(beta, se), ~ round(.x, 4)), p = round(p, 3)),
      row.names = FALSE)

write.csv(res, file.path(OUT, "22_level_ladder.csv"), row.names = FALSE)
cat("\nSaved: 22_level_ladder.csv\n")

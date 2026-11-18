# 29: numbers quoted in the text
#
# Quantities the article gives in prose rather than in a table: the
# within-country coefficients and t statistics of the mixed model
# (Section 5.1), its intraclass correlation and R2, the full body of
# Table 1, and the cross-sectional R2 with and without HIV.
#
# Reads  output/16_sysgmm_shortrun_ci.csv, 28_ladder_impact.csv,
#        driver_initial.rds, le_targets.rds, shape_pcs.rds,
#        shocks_country.rds; data/ through build_panel()
# Writes output/29_table1_full.csv, 29_within_tstats.csv

source("00_common.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

dat <- build_panel(zero_rule = "multiplicative")

# ---- 1. the t statistics in Section 5.1 --------------------
cat("=== Section 5.1: within-country coefficients of the mixed model ===\n")
dm <- dat %>% group_by(iso3c) %>%
  mutate(across(c(log_GNI, EDU_GROSS, IMM_DPT, SANITATION, ILR1, ILR2, ILR3),
                list(bw = ~ mean(.x), wn = ~ .x - mean(.x)))) %>%
  ungroup() %>% mutate(year_c = year - 2012)
mm <- lmer(LE ~ year_c + is_covid_year +
             log_GNI_bw + log_GNI_wn + EDU_GROSS_bw + EDU_GROSS_wn +
             IMM_DPT_bw + IMM_DPT_wn + SANITATION_bw + SANITATION_wn +
             ILR1_bw + ILR1_wn + ILR2_bw + ILR2_wn + ILR3_bw + ILR3_wn +
             (1 | iso3c), data = dm, REML = TRUE,
           control = lmerControl(calc.derivs = FALSE))
ss <- as.data.frame(summary(mm)$coefficients)
ss$term <- rownames(ss); rownames(ss) <- NULL
wn <- ss %>% filter(grepl("_wn$", term)) %>%
  transmute(term, beta = Estimate, se = `Std. Error`, t = `t value`) %>%
  arrange(desc(abs(t)))
cat("within-country coefficients, ordered by |t|:\n")
print(wn %>% mutate(across(c(beta, se, t), ~ round(.x, 3))))

# ---- 2. ICC and R2 -----------------------------------------
vc <- as.data.frame(VarCorr(mm))
icc <- vc$vcov[vc$grp == "iso3c"] / sum(vc$vcov)
cat(sprintf("\nintraclass correlation: %.3f\n", icc))
r2 <- performance::r2_nakagawa(mm)
cat(sprintf("marginal R2 %.3f, conditional R2 %.3f\n",
            r2$R2_marginal, r2$R2_conditional))

cs <- readRDS(file.path(OUT, "driver_initial.rds")) %>%
  inner_join(readRDS(file.path(OUT, "le_targets.rds")), by = "iso3c") %>%
  inner_join(readRDS(file.path(OUT, "shape_pcs.rds")),  by = "iso3c") %>%
  na.omit()
ctrl <- "log_GNI_init + EDU_init + IMM_init + SAN_init"

# ---- 3. the body of Table 1 --------------------------------
cat("\n\n=== Table 1: n, quantiles, and the SD decomposition ===\n")
V <- c(LE = "Life expectancy (years)", log_GNI = "Log income per head",
       z_PUB = "Public share (% CHE)", z_OOP = "Out-of-pocket (% CHE)",
       z_PNO = "Private non-OOP (% CHE)", z_EXT = "External (% CHE)",
       EDU_GROSS = "Secondary enrolment (%)", IMM_DPT = "Immunisation DPT (%)",
       SANITATION = "Basic sanitation (%)")
t1 <- bind_rows(lapply(names(V), function(v) {
  x <- dat[[v]]; g <- dat$iso3c; ok <- !is.na(x); x <- x[ok]; g <- g[ok]
  mu <- tapply(x, g, mean); cm <- as.numeric(mu[as.character(g)])
  data.frame(variable = V[[v]], n = length(x), countries = length(unique(g)),
             min = min(x), q1 = quantile(x, .25), median = median(x),
             q3 = quantile(x, .75), max = max(x), mean = mean(x),
             sd_total = sd(x), sd_between = sd(mu), sd_within = sd(x - cm),
             row.names = NULL)
}))
print(t1 %>% mutate(across(where(is.numeric), ~ round(.x, 2))))

# ---- 4. explained variance with and without the HIV control
cat("\n\n=== Cross-sectional R2, HIV control in and out ===\n")
hiv <- readRDS(file.path(OUT, "shocks_country.rds"))
csh <- cs %>% inner_join(hiv, by = "iso3c")
cat(sprintf("cross-section for the HIV comparison: %d countries\n", nrow(csh)))
hivvar <- intersect(c("HIV_peak", "HIV_pct_peak", "HIV_max"), names(csh))[1]
if (!is.na(hivvar)) {
  f0 <- as.formula(sprintf("LE_shift ~ %s + PUB_init + OOP_init + PRIV_NONOOP_init + EXT_init", ctrl))
  f1 <- as.formula(sprintf("LE_shift ~ %s + PUB_init + OOP_init + PRIV_NONOOP_init + EXT_init + %s", ctrl, hivvar))
  s0 <- summary(lm(f0, data = csh)); s1 <- summary(lm(f1, data = csh))
  cat(sprintf("HIV variable used: %s\n", hivvar))
  cat(sprintf("R2      : %.3f -> %.3f  (gain %.1f points)\n",
              s0$r.squared, s1$r.squared, 100*(s1$r.squared - s0$r.squared)))
  cat(sprintf("adj. R2 : %.3f -> %.3f  (gain %.1f points)\n",
              s0$adj.r.squared, s1$adj.r.squared,
              100*(s1$adj.r.squared - s0$adj.r.squared)))
} else cat("HIV peak variable not found in shocks_country.rds; names:",
           paste(names(csh), collapse = ", "), "\n")

# ---- 5. headline coefficients at full precision -------------
cat("\n\n=== Coefficients at four decimals ===\n")
sr  <- read.csv(file.path(OUT, "16_sysgmm_shortrun_ci.csv"))
lad <- read.csv(file.path(OUT, "28_ladder_impact.csv"))
cat(sprintf("static one-way FE on PUB : %+.4f\n",
            lad$beta[grepl("country effects only", lad$rung)]))
cat(sprintf("System GMM on PUB        : %+.4f\n",
            sr$beta[sr$term == "PUB"]))
cat(sprintf("System GMM on IMM        : %+.4f\n",
            sr$beta[sr$term == "IMM_DPT"]))
cat(sprintf("persistence alpha        : %+.4f\n",
            sr$beta[sr$term == "lag(LE, 1)"]))

write.csv(t1,   file.path(OUT, "29_table1_full.csv"), row.names = FALSE)
write.csv(wn,   file.path(OUT, "29_within_tstats.csv"), row.names = FALSE)
cat("\nSaved: 29_table1_full.csv, 29_within_tstats.csv\n")

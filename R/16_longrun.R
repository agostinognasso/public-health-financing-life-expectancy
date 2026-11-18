# 16: System GMM of Table 3 and its long-run multipliers
#
# The static within coefficient is comparable with the long-run
# multiplier beta/(1-alpha) of the dynamic model, not with the
# one-year impact coefficient. The script fits the System GMM of
# Table 3 and reports, for every regressor, beta_LR with a
# delta-method SE, Fieller sets at 90 and 95 per cent, and the
# confidence levels at which each set stops being an interval and
# becomes the whole line. Cov(beta, alpha) is saved for script 19.
#
# Reads  data/ through build_panel()
# Writes output/16_sysgmm_shortrun_ci.csv, 16_sysgmm_longrun.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
cat("Panel:", nrow(dat), "obs,", n_distinct(dat$iso3c), "countries,",
    min(dat$year), "-", max(dat$year), "\n")

pd <- pdata.frame(dat, index = c("iso3c", "year"))

f_raw <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
              IMM_DPT + SANITATION | lag(LE, 2:5)

# ---- 1. the System GMM of Table 3 --------------------------
m_sys <- pgmm(f_raw, data = pd, effect = "twoways",
              model = "twosteps", transformation = "ld", collapse = TRUE)
s_sys <- summary(m_sys, robust = TRUE)

cat("=== System GMM, raw shares ===\n")
print(round(s_sys$coefficients, 4))

b <- coef(m_sys)
V <- pgmm_vcov(m_sys)
stopifnot(all(names(b) %in% rownames(V)))

# ---- 2. short-run coefficients with 95% CI ----------------
sr <- data.frame(
  term = names(b),
  beta = as.numeric(b),
  se   = sqrt(diag(V)[names(b)]),
  row.names = NULL
) %>%
  mutate(z = beta / se,
         p = 2 * pnorm(-abs(z)),
         lo = beta - 1.96 * se,
         hi = beta + 1.96 * se)

cat("\n=== Short-run coefficients with 95% CI ===\n")
print(sr %>% mutate(across(c(beta, se, lo, hi), ~ round(.x, 4)),
                    z = round(z, 2), p = fmt_p(p)))

# ---- 3. long-run multipliers -------------------------------
lr <- longrun(b, V, alpha_name = "lag(LE, 1)")
alpha <- attr(lr, "alpha"); se_a <- attr(lr, "se_alpha")

# the same sets at 90 per cent, the levels at which they change
# shape, and Cov(beta, alpha) for the equivalence tests of script 19
lr90 <- longrun(b, V, alpha_name = "lag(LE, 1)", level = 0.90)
stopifnot(identical(lr$term, lr90$term))
lr$fieller90_lo   <- lr90$fieller_lo
lr$fieller90_hi   <- lr90$fieller_hi
lr$fieller90_type <- lr90$fieller_type
lr <- merge(lr, fieller_levels(b, V, "lag(LE, 1)"), by = "term",
            sort = FALSE)
lr <- lr[match(lr90$term, lr$term), ]
rownames(lr) <- NULL
lr$cov_alpha <- unname(V[lr$term, "lag(LE, 1)"])

cat("\n=== Persistence ===\n")
cat(sprintf("alpha = %.4f (SE %.4f)\n", alpha, se_a))
cat(sprintf("1 - alpha = %.4f (SE %.4f), 95%% CI [%.4f, %.4f]\n",
            1 - alpha, se_a, 1 - alpha - 1.96 * se_a, 1 - alpha + 1.96 * se_a))
cat(sprintf("implied half-life = %.1f years\n", log(0.5) / log(alpha)))
if (1 - alpha - 1.96 * se_a <= 0)
  cat("note: the 95% CI for (1-alpha) contains zero, so the denominator\n",
      "      of the long-run multiplier may be zero.\n", sep = "")

cat("\n=== Long-run multipliers beta/(1-alpha) ===\n")
print(lr %>% mutate(across(c(beta_SR, se_SR, beta_LR, se_LR,
                             lo_delta, hi_delta), ~ round(.x, 4)),
                    fieller_lo = round(fieller_lo, 4),
                    fieller_hi = round(fieller_hi, 4)) %>%
        dplyr::select(term:fieller_type))

cat("\n=== Fieller sets at 90 and 95 per cent ===\n")
cat(sprintf("Every set is an interval up to the %.1f%% level, where 1 - alpha\n",
            100 * lr$level_bounded[1]))
cat("stops being significantly different from zero; above it a set is\n")
cat("two half-lines, and above level_wholeline the whole line.\n\n")
for (tn in KEY_X) {
  r <- lr[lr$term == tn, ]
  cat(sprintf("%-11s beta_LR %+.4f\n", tn, r$beta_LR))
  cat(sprintf("            90%%: %s\n",
              fmt_fieller(r$fieller90_lo, r$fieller90_hi, r$fieller90_type)))
  cat(sprintf("            95%%: %s\n",
              fmt_fieller(r$fieller_lo, r$fieller_hi, r$fieller_type)))
  cat(sprintf("            whole line from the %.2f%% level\n",
              100 * r$level_wholeline))
}

# ---- 4. static within against the long-run interval -------
# static within estimator, country effects only
m_fe <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
              SANITATION, data = pd, model = "within",
            effect = "individual")
b_static <- unname(coef(m_fe)["PUB"])
pub <- lr[lr$term == "PUB", ]

cat("\n=== Static within estimate against the long-run intervals ===\n")
cat(sprintf("static FE within beta(PUB)      = %+.4f\n", b_static))
cat(sprintf("long-run beta_LR(PUB)           = %+.4f (delta SE %.4f)\n",
            pub$beta_LR, pub$se_LR))
cat(sprintf("delta-method 95%% CI             = [%+.4f, %+.4f]  contains it? %s\n",
            pub$lo_delta, pub$hi_delta,
            ifelse(b_static >= pub$lo_delta & b_static <= pub$hi_delta, "yes", "no")))
cat(sprintf("Fieller 95%% set                 = %s\n",
            fmt_fieller(pub$fieller_lo, pub$fieller_hi, pub$fieller_type)))
cat(sprintf("Fieller 90%% set                 = %s  contains it? %s\n",
            fmt_fieller(pub$fieller90_lo, pub$fieller90_hi, pub$fieller90_type),
            ifelse(b_static >= pub$fieller90_lo & b_static <= pub$fieller90_hi,
                   "yes", "no")))

cat("\n--- implied effect of a 20 percentage point shift in the public share ---\n")
cat(sprintf("point estimate : %+.2f years\n", 20 * pub$beta_LR))
cat(sprintf("delta 95%% CI   : [%+.2f, %+.2f] years\n",
            20 * pub$lo_delta, 20 * pub$hi_delta))

imm <- lr[lr$term == "IMM_DPT", ]
cat("\n--- immunisation: 20 pp of coverage ---\n")
cat(sprintf("long-run beta_LR = %+.4f years/pp, 20 pp -> %+.2f years\n",
            imm$beta_LR, 20 * imm$beta_LR))
cat(sprintf("delta 95%% CI on 20 pp: [%+.2f, %+.2f] years\n",
            20 * imm$lo_delta, 20 * imm$hi_delta))
cat(sprintf("10 pp -> %+.2f years [%+.2f, %+.2f]\n",
            10 * imm$beta_LR, 10 * imm$lo_delta, 10 * imm$hi_delta))
cat(sprintf("Fieller 95%% set on 10 pp: %s years\n",
            fmt_fieller(10 * imm$fieller_lo, 10 * imm$fieller_hi,
                        imm$fieller_type, digits = 2)))
cat(sprintf("Fieller 90%% set on 10 pp: %s years\n",
            fmt_fieller(10 * imm$fieller90_lo, 10 * imm$fieller90_hi,
                        imm$fieller90_type, digits = 2)))

write.csv(sr, file.path(OUT, "16_sysgmm_shortrun_ci.csv"), row.names = FALSE)
write.csv(lr, file.path(OUT, "16_sysgmm_longrun.csv"), row.names = FALSE)

cat("\nSaved: 16_sysgmm_shortrun_ci.csv, 16_sysgmm_longrun.csv\n")

# 19: within/between variation, minimum detectable effects and
#     equivalence tests for the public share
#
#   1. within and between standard deviations of every variable
#      (Table 1)
#   2. the variation in the public share left after country and year
#      effects and the controls
#   3. the minimum detectable effect at 80 per cent power and two
#      one-sided equivalence tests (TOST), at impact and in the long
#      run. In the long run both are given by the delta method and
#      in a form that does not rely on it: the impact MDE divided by
#      1 - alpha, and TOST on the Fieller pivot (Sections 4.3, 5.3).
#
# Reads  output/16_sysgmm_shortrun_ci.csv, 16_sysgmm_longrun.csv;
#        data/ through build_panel()
# Writes output/19_within_between_sd.csv, 19_tost.csv, 19_mde.csv,
#        19_mde_longrun.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c", "year"))

# static within estimate with country effects, the reference value
# for the one-SD move below and for one of the equivalence margins
m_fe1 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
               SANITATION, data = pd, model = "within",
             effect = "individual")
b_fe1   <- unname(coef(m_fe1)["PUB"])
b_naive <- abs(b_fe1)

VARS <- c(LE = "Life expectancy (years)",
          log_GNI = "Log income per head",
          z_PUB = "Public share (% CHE)",
          z_OOP = "Out-of-pocket (% CHE)",
          z_PNO = "Private non-OOP (% CHE)",
          z_EXT = "External (% CHE)",
          EDU_GROSS = "Secondary enrolment (%)",
          IMM_DPT = "Immunisation, DPT (%)",
          SANITATION = "Basic sanitation (%)")

# ---- 1. within / between decomposition (Table 1) -----------
cat("=== Within and between standard deviations ===\n")
dec <- lapply(names(VARS), function(v) {
  x <- dat[[v]]; g <- dat$iso3c
  ok <- !is.na(x)
  x <- x[ok]; g <- g[ok]
  mu <- tapply(x, g, mean)
  cm <- as.numeric(mu[as.character(g)])
  data.frame(variable = VARS[[v]], n = length(x),
             n_country = length(unique(g)),
             mean = mean(x), sd_total = sd(x),
             sd_between = sd(mu), sd_within = sd(x - cm),
             share_within = sd(x - cm)^2 / sd(x)^2,
             row.names = NULL)
})
dec <- bind_rows(dec)
print(dec %>% mutate(across(c(mean, sd_total, sd_between, sd_within),
                            ~ round(.x, 3)),
                     share_within = round(share_within, 3)))

sd_w_pub <- dec$sd_within[dec$variable == "Public share (% CHE)"]
cat(sprintf("\nWithin-country SD of the public share: %.2f percentage points\n",
            sd_w_pub))
cat(sprintf("A one-within-SD move at the static estimate (%+.4f): %+.3f years\n",
            b_fe1, b_fe1 * sd_w_pub))

# ---- 2. identifying variation ------------------------------
cat("\n=== Variation in the public share after the fixed effects ===\n")
aux <- lm(z_PUB ~ factor(iso3c) + factor(year), data = dat)
r_pub <- residuals(aux)
cat(sprintf("SD of the public share: raw %.2f, within-country %.2f,\n",
            sd(dat$z_PUB), sd_w_pub))
cat(sprintf("after country and year effects %.2f pp (R2 of the auxiliary = %.3f)\n",
            sd(r_pub), summary(aux)$r.squared))
aux2 <- lm(z_PUB ~ factor(iso3c) + factor(year) + log_GNI + EDU_GROSS +
                   IMM_DPT + SANITATION + z_OOP, data = dat)
cat(sprintf("after the full control set %.2f pp (R2 = %.3f)\n",
            sd(residuals(aux2)), summary(aux2)$r.squared))

# ---- 3. minimum detectable effect and equivalence tests -----
sr <- read.csv(file.path(OUT, "16_sysgmm_shortrun_ci.csv"))
lr <- read.csv(file.path(OUT, "16_sysgmm_longrun.csv"))

se_pub_sr <- sr$se[sr$term == "PUB"]
b_pub_sr  <- sr$beta[sr$term == "PUB"]
se_pub_lr <- lr$se_LR[lr$term == "PUB"]
b_pub_lr  <- lr$beta_LR[lr$term == "PUB"]
a_hat     <- sr$beta[sr$term == "lag(LE, 1)"]
se_a      <- sr$se[sr$term == "lag(LE, 1)"]
cov_pub_a <- lr$cov_alpha[lr$term == "PUB"]

mde <- function(se, power = 0.80, alpha = 0.05)
  (qnorm(1 - alpha / 2) + qnorm(power)) * se

cat("\n=== Minimum detectable effect on the public share ===\n")
cat(sprintf("short run: SE = %.4f -> MDE(80%%) = %.4f years/pp",
            se_pub_sr, mde(se_pub_sr)))
cat(sprintf("  = %.2f years over 20 pp\n", 20 * mde(se_pub_sr)))
cat(sprintf("long run : SE = %.4f -> MDE(80%%) = %.4f years/pp",
            se_pub_lr, mde(se_pub_lr)))
cat(sprintf("  = %.2f years over 20 pp\n", 20 * mde(se_pub_lr)))
cat(sprintf("\nlong-run MDE / |static within estimate| = %.4f / %.4f = %.1f\n",
            mde(se_pub_lr), b_naive, mde(se_pub_lr) / b_naive))

# The delta-method SE treats beta/(1-alpha) as normal, which it is
# not when 1-alpha is imprecise. With 1-alpha non-zero,
# beta/(1-alpha) = 0 exactly when beta = 0, and the Fieller set
# excludes zero exactly when the impact coefficient is significant,
# so the test of a zero long-run effect is the test of a zero impact
# effect and, at a given alpha, the long-run MDE is the impact MDE
# divided by 1 - alpha. The MDE is reported along the interval for
# alpha.
cat("\n=== Long-run MDE as the impact MDE over (1 - alpha) ===\n")
z90 <- qnorm(0.95); z95 <- qnorm(0.975)
mde_lr <- data.frame(
  alpha_at = c("point estimate", "90% CI lower", "90% CI upper",
               "95% CI lower", "95% CI upper"),
  alpha    = a_hat + c(0, -z90, z90, -z95, z95) * se_a)
mde_lr$one_minus_alpha <- 1 - mde_lr$alpha
mde_lr$mde80 <- ifelse(mde_lr$one_minus_alpha > 0,
                       mde(se_pub_sr) / mde_lr$one_minus_alpha, Inf)
mde_lr$mde80_over_20pp <- 20 * mde_lr$mde80
print(mde_lr %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
cat(sprintf("\nAt alpha-hat the long-run MDE is %.4f years/pp (%.2f years over 20 pp);\n",
            mde_lr$mde80[1], mde_lr$mde80_over_20pp[1]))
cat(sprintf("the delta-method approximation above gives %.4f.\n", mde(se_pub_lr)))
cat(sprintf("Over the 90%% interval for alpha it runs from %.4f to %.4f.\n",
            mde_lr$mde80[2], mde_lr$mde80[3]))
cat(sprintf("Over the 95%% interval it is at least %.4f and has no upper bound,\n",
            mde_lr$mde80[4]))
cat("because that interval reaches alpha = 1.\n")

# TOST: two one-sided tests against +/- margin
tost <- function(b, se, margin) {
  t_lo <- (b - (-margin)) / se      # H0: b <= -margin
  t_hi <- (b -   margin)  / se      # H0: b >=  margin
  p_lo <- 1 - pnorm(t_lo)
  p_hi <- pnorm(t_hi)
  data.frame(margin = margin, beta = b, se = se,
             p_lower = p_lo, p_upper = p_hi,
             p_tost = max(p_lo, p_hi),
             equivalent = max(p_lo, p_hi) < 0.05)
}

margins <- sort(c(0.0125, 0.025, 0.05, round(b_naive, 4)))
cat("\n=== Equivalence tests (TOST), 5% level ===\n")
cat("Margins in years of life expectancy per percentage point. Three are\n")
cat("0.25, 0.50 and 1.00 years over a 20 pp change in the public share;\n")
cat(sprintf("the fourth, %.4f, is the static within estimate.\n\n",
            b_naive))

# The same tests on the long-run multiplier with the Fieller pivot
# in place of the delta-method SE. Under H0: beta/(1-alpha) = th0,
# (beta - th0 (1-alpha)) / sd(beta - th0 (1-alpha)) is standard
# normal however imprecise 1-alpha is. The one-sided tests are
# signed on the assumption that 1 - alpha > 0. Equivalence at 5 per
# cent then holds when the 90 per cent Fieller set of script 16 lies
# inside the margins; the stopifnot() below checks that the two agree.
den_a <- 1 - a_hat
f_t <- function(th0)
  (b_pub_sr - th0 * den_a) /
    sqrt(se_pub_sr^2 + 2 * th0 * cov_pub_a + th0^2 * se_a^2)
tost_fieller <- function(margin) {
  p_lo <- 1 - pnorm(f_t(-margin))   # H0: long-run effect <= -margin
  p_hi <- pnorm(f_t(margin))        # H0: long-run effect >=  margin
  data.frame(margin = margin, beta = b_pub_lr, se = NA_real_,
             p_lower = p_lo, p_upper = p_hi,
             p_tost = max(p_lo, p_hi),
             equivalent = max(p_lo, p_hi) < 0.05)
}

tt <- bind_rows(
  bind_rows(lapply(margins, function(m) tost(b_pub_sr, se_pub_sr, m))) %>%
    mutate(horizon = "short run (impact)", method = "Wald", .before = 1),
  bind_rows(lapply(margins, function(m) tost(b_pub_lr, se_pub_lr, m))) %>%
    mutate(horizon = "long run (beta/(1-alpha))", method = "delta",
           .before = 1),
  bind_rows(lapply(margins, tost_fieller)) %>%
    mutate(horizon = "long run (beta/(1-alpha))", method = "Fieller",
           .before = 1)
)
print(tt %>% mutate(across(c(beta, se), ~ round(.x, 4)),
                    across(c(p_lower, p_upper, p_tost), ~ round(.x, 4))))

f90 <- lr[lr$term == "PUB", c("fieller90_lo", "fieller90_hi", "fieller90_type")]
inside <- f90$fieller90_type == "bounded" &
  -margins < f90$fieller90_lo & f90$fieller90_hi < margins
stopifnot(identical(inside, tt$equivalent[tt$method == "Fieller"]))
cat(sprintf("\n90%% Fieller set for the long-run multiplier: %s;\n",
            fmt_fieller(f90$fieller90_lo, f90$fieller90_hi, f90$fieller90_type)))
cat("the Fieller TOST agrees with it at every margin.\n")

cat("\nEquivalence at 5 per cent:\n")
reading <- function(h, m, lab) {
  eq <- tt %>% filter(horizon == h, method == m, equivalent)
  if (nrow(eq))
    cat(sprintf("- %-22s: established at margins >= %.4f\n",
                lab, min(eq$margin))) else
    cat(sprintf("- %-22s: not established at any margin\n", lab))
}
reading("short run (impact)", "Wald", "short run")
reading("long run (beta/(1-alpha))", "delta", "long run, delta method")
reading("long run (beta/(1-alpha))", "Fieller", "long run, Fieller")

write.csv(dec, file.path(OUT, "19_within_between_sd.csv"), row.names = FALSE)
write.csv(tt,  file.path(OUT, "19_tost.csv"), row.names = FALSE)
write.csv(data.frame(
  horizon = c("short run", "long run"),
  beta = c(b_pub_sr, b_pub_lr), se = c(se_pub_sr, se_pub_lr),
  mde80 = c(mde(se_pub_sr), mde(se_pub_lr)),
  mde80_over_20pp = 20 * c(mde(se_pub_sr), mde(se_pub_lr)),
  mde80_exact = c(mde(se_pub_sr), mde_lr$mde80[1])),
  file.path(OUT, "19_mde.csv"), row.names = FALSE)
write.csv(mde_lr, file.path(OUT, "19_mde_longrun.csv"), row.names = FALSE)
cat("\nSaved: 19_within_between_sd.csv, 19_tost.csv, 19_mde.csv,\n")
cat("19_mde_longrun.csv\n")

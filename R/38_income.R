# 38: the public-share coefficient by income group
#
# One pre-specified split: countries classified by the World Bank as
# low or lower-middle income in 2000, the first year of the window,
# against those classified upper-middle or high income. The
# classification is fixed at 2000 so that it does not depend on later
# income growth. The public share enters together with its
# interaction with the low and lower-middle indicator, which is
# constant within country, in the reference System GMM of Table 3
# (drivers strictly exogenous, the interaction treated like the
# public share) and in the two-way static fixed effects rung. The
# script reports the coefficient for each group, their difference,
# the smallest equivalence margin that each group and the difference
# establish at the 5 per cent level, and the long-run multiplier of
# each group with Fieller sets (Section 5.3, Table 5).
#
# Reads  data/df_panel_drivers.rds, and data/ through build_panel()
# Writes output/38_income.csv, 38_income_groups.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

LOW_LABELS <- c("Low income", "Lower-middle income")

drivers <- readRDS(file.path(DATA, "df_panel_drivers.rds"))
g2000 <- drivers %>%
  filter(year == 2000) %>%
  transmute(iso3c = iso3c.x, label2000 = GNI_label) %>%
  distinct()
stopifnot(!any(duplicated(g2000$iso3c)))

dat <- build_panel(zero_rule = "multiplicative") %>%
  left_join(g2000, by = "iso3c")
cat(sprintf("panel: %d obs, %d countries; income group in 2000 missing for %d obs\n",
            nrow(dat), n_distinct(dat$iso3c), sum(is.na(dat$label2000))))
stopifnot(!any(is.na(dat$label2000)))

dat <- dat %>%
  mutate(LOWINC     = as.integer(label2000 %in% LOW_LABELS),
         PUB_LOWINC = PUB * LOWINC)

cat("\n=== Countries by income group in 2000 ===\n")
print(table(dat %>% distinct(iso3c, label2000) %>% pull(label2000)))

# identifying variation of the public share in each group
grp <- dat %>%
  group_by(iso3c) %>% mutate(pub_w = PUB - mean(PUB)) %>% ungroup() %>%
  group_by(LOWINC) %>%
  summarise(countries = n_distinct(iso3c), obs = n(),
            mean_PUB = mean(PUB), sd_within_PUB = sd(pub_w),
            .groups = "drop") %>%
  mutate(group = ifelse(LOWINC == 1, "low and lower-middle income in 2000",
                        "upper-middle and high income in 2000"))
cat("\n=== Groups ===\n")
print(grp %>% mutate(across(c(mean_PUB, sd_within_PUB), ~ round(.x, 2))))

pd <- pdata.frame(dat, index = c("iso3c", "year"))

# ---- reference System GMM with the interaction --------------
f_int <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + PUB_LOWINC +
              OOP_share + IMM_DPT + SANITATION | lag(LE, 2:5)
m  <- pgmm(f_int, data = pd, effect = "twoways", model = "twosteps",
           transformation = "ld", collapse = TRUE)
md <- pgmm(f_int, data = pd, effect = "twoways", model = "twosteps",
           transformation = "d", collapse = TRUE)
b <- coef(m); V <- pgmm_vcov(m)
sg <- sargan(m); sgd <- sargan(md)
ar2 <- mtest(m, order = 2)
dih <- as.numeric(sg$statistic - sgd$statistic)
dih_df <- as.numeric(sg$parameter - sgd$parameter)
cat("\n=== System GMM with the interaction ===\n")
print(round(summary(m, robust = TRUE)$coefficients, 4))
cat(sprintf("instruments %d, AR(2) p = %.3f, Hansen p = %.3f, Difference-in-Hansen p = %.3f\n",
            as.integer(sg$parameter) + length(b), ar2$p.value, sg$p.value,
            pchisq(dih, dih_df, lower.tail = FALSE)))

# group coefficients as linear combinations, with their covariance
# with the persistence parameter, for the Fieller sets
A <- rbind(alpha = c("lag(LE, 1)" = 1, PUB = 0, PUB_LOWINC = 0),
           high  = c(0, 1, 0),
           low   = c(0, 1, 1))
nm <- c("lag(LE, 1)", "PUB", "PUB_LOWINC")
bg <- as.numeric(A %*% b[nm]); names(bg) <- c("lag(LE, 1)", "PUB_high", "PUB_low")
Vg <- A %*% V[nm, nm] %*% t(A); dimnames(Vg) <- list(names(bg), names(bg))
lr95 <- longrun(bg, Vg, "lag(LE, 1)", level = 0.95)
lr90 <- longrun(bg, Vg, "lag(LE, 1)", level = 0.90)

z90 <- qnorm(0.95)
gmm_rows <- data.frame(
  estimator = "System GMM",
  group = c("upper-middle and high income in 2000",
            "low and lower-middle income in 2000",
            "difference, low minus high"),
  beta = c(bg["PUB_high"], bg["PUB_low"], b["PUB_LOWINC"]),
  se   = c(sqrt(Vg["PUB_high", "PUB_high"]), sqrt(Vg["PUB_low", "PUB_low"]),
           sqrt(V["PUB_LOWINC", "PUB_LOWINC"])),
  alpha = unname(b["lag(LE, 1)"]),
  se_alpha = sqrt(V["lag(LE, 1)", "lag(LE, 1)"]),
  beta_LR = c(lr95$beta_LR, NA),
  fieller90_lo = c(lr90$fieller_lo, NA), fieller90_hi = c(lr90$fieller_hi, NA),
  fieller90_type = c(lr90$fieller_type, NA),
  fieller95_lo = c(lr95$fieller_lo, NA), fieller95_hi = c(lr95$fieller_hi, NA),
  fieller95_type = c(lr95$fieller_type, NA),
  n_inst = as.integer(sg$parameter) + length(b),
  ar2_p = ar2$p.value, hansen_p = sg$p.value,
  diffhansen_p = pchisq(dih, dih_df, lower.tail = FALSE),
  row.names = NULL)

# ---- two-way static fixed effects ---------------------------
fe <- plm(LE ~ log_GNI + EDU_GROSS + PUB + PUB_LOWINC + OOP_share +
            IMM_DPT + SANITATION, data = pd, model = "within",
          effect = "twoways")
Vf <- plm::vcovHC(fe, cluster = "group")
bf <- coef(fe)
L <- rbind(c(1, 0), c(1, 1), c(0, 1))
nmf <- c("PUB", "PUB_LOWINC")
bfg <- as.numeric(L %*% bf[nmf])
Vfg <- L %*% Vf[nmf, nmf] %*% t(L)
fe_rows <- data.frame(
  estimator = "Two-way fixed effects",
  group = gmm_rows$group, beta = bfg, se = sqrt(diag(Vfg)),
  row.names = NULL)

res <- bind_rows(gmm_rows, fe_rows) %>%
  mutate(z = beta / se, p = 2 * pnorm(-abs(z)),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         margin_5pct = abs(beta) + z90 * se)

cat("\n=== Public share by income group in 2000 ===\n")
print(res %>% mutate(across(c(beta, se, lo, hi, margin_5pct), ~ round(.x, 4)),
                     p = round(p, 3)) %>%
        dplyr::select(estimator, group, beta, se, lo, hi, p, margin_5pct))

cat("\n=== Long-run multiplier by group, System GMM ===\n")
for (i in 1:2) {
  r <- gmm_rows[i, ]
  cat(sprintf("%-40s beta_LR %+.4f\n", r$group, r$beta_LR))
  cat(sprintf("%-40s 90%%: %s\n", "",
              fmt_fieller(r$fieller90_lo, r$fieller90_hi, r$fieller90_type, 3)))
  cat(sprintf("%-40s 95%%: %s\n", "",
              fmt_fieller(r$fieller95_lo, r$fieller95_hi, r$fieller95_type, 3)))
}

write.csv(res, file.path(OUT, "38_income.csv"), row.names = FALSE)
write.csv(grp, file.path(OUT, "38_income_groups.csv"), row.names = FALSE)
cat("\nSaved: 38_income.csv, 38_income_groups.csv\n")

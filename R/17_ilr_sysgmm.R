# 17: System GMM on the ILR coordinates
#
# The System GMM of Table 3 on the raw shares (A) and the same model
# with the three ILR coordinates in place of the shares (B), each with
# AR(1)-AR(3), Hansen J and Difference-in-Hansen; the long-run
# multipliers of B; and the map that converts a move along ILR1 into
# percentage points of the public share (Section 3.2).
#
# Reads  data/ through build_panel()
# Writes output/17_sysgmm_raw_vs_ilr.csv, 17_ilr_longrun.csv,
#        17_ilr_interpretive_map.csv, 17_diagnostics.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c", "year"))
cat("Panel:", nrow(dat), "obs,", n_distinct(dat$iso3c), "countries\n\n")

CORE <- c("lag(LE, 1)", "log_GNI", "EDU_GROSS", "PUB", "OOP_share",
          "IMM_DPT", "SANITATION", "ILR1", "ILR2", "ILR3")

fit_sys <- function(f) pgmm(f, data = pd, effect = "twoways",
                            model = "twosteps", transformation = "ld",
                            collapse = TRUE)

# ------------------------------------------------------------
# diagnostics, including AR(3) and Difference-in-Hansen with df
# ------------------------------------------------------------
diagnose <- function(m, m_diff, label) {
  s  <- summary(m, robust = TRUE)
  sg <- sargan(m)
  ar <- lapply(1:3, function(k)
    tryCatch(mtest(m, order = k), error = function(e) NULL))
  n_sys  <- length(m$W[[1]][1, ])
  n_inst_sys  <- dim(m$W[[1]])[2]
  sg_d <- sargan(m_diff)
  n_inst_diff <- dim(m_diff$W[[1]])[2]
  C  <- unname(sg$statistic - sg_d$statistic)
  df <- unname(sg$parameter - sg_d$parameter)
  pC <- if (!is.na(df) && df > 0) pchisq(C, df, lower.tail = FALSE) else NA

  cat("\n--- diagnostics:", label, "---\n")
  cat(sprintf("Hansen J        = %.3f  df = %d  p = %.4f\n",
              sg$statistic, sg$parameter, sg$p.value))
  for (k in 1:3) if (!is.null(ar[[k]]))
    cat(sprintf("Arellano-Bond AR(%d) = %+.3f  p = %.4f\n",
                k, ar[[k]]$statistic, ar[[k]]$p.value))
  cat(sprintf("Difference-in-Hansen C = %.3f  df = %d  p = %.4f\n", C, df, pC))
  cat(sprintf("instruments: system %d, difference %d\n",
              n_inst_sys, n_inst_diff))
  data.frame(spec = label,
             hansen_J = sg$statistic, hansen_df = sg$parameter,
             hansen_p = sg$p.value,
             ar1 = ar[[1]]$statistic, ar1_p = ar[[1]]$p.value,
             ar2 = ifelse(is.null(ar[[2]]), NA, ar[[2]]$statistic),
             ar2_p = ifelse(is.null(ar[[2]]), NA, ar[[2]]$p.value),
             ar3 = ifelse(is.null(ar[[3]]), NA, ar[[3]]$statistic),
             ar3_p = ifelse(is.null(ar[[3]]), NA, ar[[3]]$p.value),
             dhansen_C = C, dhansen_df = df, dhansen_p = pC,
             n_inst_sys = n_inst_sys, n_inst_diff = n_inst_diff,
             row.names = NULL)
}

tidy_core <- function(m, label) {
  b <- coef(m); V <- pgmm_vcov(m)
  keep <- intersect(names(b), CORE)
  data.frame(spec = label, term = keep,
             beta = b[keep], se = sqrt(diag(V)[keep]), row.names = NULL) %>%
    mutate(z = beta / se, p = 2 * pnorm(-abs(z)),
           lo = beta - 1.96 * se, hi = beta + 1.96 * se)
}

# ---- specification A: raw shares ---------------------------
fA  <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
            IMM_DPT + SANITATION | lag(LE, 2:5)
fAd <- fA
mA  <- fit_sys(fA)
mAd <- pgmm(fAd, data = pd, effect = "twoways", model = "twosteps",
            transformation = "d", collapse = TRUE)

# ---- specification B: ILR coordinates ----------------------
fB  <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + ILR1 + ILR2 + ILR3 +
            IMM_DPT + SANITATION | lag(LE, 2:5)
mB  <- fit_sys(fB)
mBd <- pgmm(fB, data = pd, effect = "twoways", model = "twosteps",
            transformation = "d", collapse = TRUE)

cat("=== A. System GMM, raw shares ===\n")
print(round(summary(mA, robust = TRUE)$coefficients[1:7, ], 4))
dA <- diagnose(mA, mAd, "raw shares")

cat("\n\n=== B. System GMM, ILR coordinates ===\n")
print(round(summary(mB, robust = TRUE)$coefficients[1:8, ], 4))
dB <- diagnose(mB, mBd, "ILR coordinates")

tab <- bind_rows(tidy_core(mA, "raw shares"), tidy_core(mB, "ILR"))
cat("\n=== Core coefficients, both specifications, 95% CI ===\n")
print(tab %>% mutate(across(c(beta, se, lo, hi), ~ round(.x, 4)),
                     z = round(z, 2), p = fmt_p(p)))

# ---- long-run multipliers for the ILR specification --------
lrB <- longrun(coef(mB), pgmm_vcov(mB), "lag(LE, 1)") %>%
  filter(term %in% CORE)
cat("\n=== Long-run multipliers, ILR specification ===\n")
print(lrB %>% mutate(across(c(beta_SR, se_SR, beta_LR, se_LR,
                              lo_delta, hi_delta, fieller_lo, fieller_hi),
                            ~ round(.x, 4))))

# ------------------------------------------------------------
# From ILR1 to percentage points of the public share
#
# Holding ILR2 and ILR3 fixed keeps the ratios among OOP, PNO and
# EXT fixed, so the geometric mean of the non-public parts is
# proportional to (1 - u) where u is the public share as a
# fraction. Then
#     ILR1 = sqrt(3/4) * (logit(u) - log c)
#     d ILR1 / du = sqrt(3/4) / (u (1-u))
# so a one percentage point rise in the public share (du = 0.01),
# with the non-public composition held fixed, is a move of
#     sqrt(3/4) * 0.01 / (u (1-u))
# along ILR1, and the effect in years per percentage point is
# beta_ILR1 times that factor.
# ------------------------------------------------------------
u_q <- quantile(dat$z_PUB / 100, c(0.25, 0.5, 0.75))
k_of_u <- function(u) sqrt(3 / 4) * 0.01 / (u * (1 - u))

bI  <- tab$beta[tab$spec == "ILR" & tab$term == "ILR1"]
seI <- tab$se[tab$spec == "ILR"  & tab$term == "ILR1"]
lrI <- lrB$beta_LR[lrB$term == "ILR1"]; lrIse <- lrB$se_LR[lrB$term == "ILR1"]

map <- data.frame(
  quantile   = names(u_q),
  public_share = round(100 * as.numeric(u_q), 1),
  dILR1_per_pp = round(k_of_u(as.numeric(u_q)), 4),
  effect_SR_per_pp = round(bI * k_of_u(as.numeric(u_q)), 5),
  SR_lo = round((bI - 1.96 * seI) * k_of_u(as.numeric(u_q)), 5),
  SR_hi = round((bI + 1.96 * seI) * k_of_u(as.numeric(u_q)), 5),
  effect_LR_per_pp = round(lrI * k_of_u(as.numeric(u_q)), 5),
  LR_lo = round((lrI - 1.96 * lrIse) * k_of_u(as.numeric(u_q)), 5),
  LR_hi = round((lrI + 1.96 * lrIse) * k_of_u(as.numeric(u_q)), 5)
)
cat("\n=== One percentage point of public share, in ILR1 terms ===\n")
cat("(years of life expectancy per pp, non-public composition held fixed)\n")
print(map)

cat("\nThe raw-share System GMM gives",
    sprintf("%+.4f", tab$beta[tab$spec == "raw shares" & tab$term == "PUB"]),
    "years/pp (short run).\n")

write.csv(tab,  file.path(OUT, "17_sysgmm_raw_vs_ilr.csv"), row.names = FALSE)
write.csv(lrB,  file.path(OUT, "17_ilr_longrun.csv"),        row.names = FALSE)
write.csv(map,  file.path(OUT, "17_ilr_interpretive_map.csv"), row.names = FALSE)
write.csv(bind_rows(dA, dB), file.path(OUT, "17_diagnostics.csv"), row.names = FALSE)
cat("\nSaved: 17_sysgmm_raw_vs_ilr.csv, 17_ilr_longrun.csv,\n")
cat("17_ilr_interpretive_map.csv, 17_diagnostics.csv\n")

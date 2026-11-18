# 31: the AR(3) rejection (Sections 5.2 and 6.7)
#
# The Arellano-Bond moment conditions need no second-order serial
# correlation in the differenced residuals; AR(3) is not among them.
# With eps_it serially uncorrelated, Delta eps_it is MA(1), so AR(2)
# and higher orders should be insignificant in differences. An AR(3)
# rejection points to E(eps_it eps_i,t-2) != 0, under which lags 2
# and 3 of the outcome are not valid instruments and the set has to
# start at t-4.
#
# The script re-estimates the System GMM and the difference GMM with
# instruments starting at t-2 to t-5, adds a second lag of the
# outcome, repeats the ILR specification at two depths, and checks
# the AR tests for a pure AR(1) in LE, by end year of the window, and
# the autocorrelation of the differenced within residuals.
#
# Reads  data/ through build_panel()
# Writes output/31_ar3_windows.csv, 31_ar3_depth_sysgmm.csv,
#        31_ar3_depth_diffgmm.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c", "year"))
FOC <- "PUB"

fit_depth <- function(lo, hi, transformation = "ld") {
  fg <- as.formula(sprintf(
    "LE ~ lag(LE,1) + log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
     SANITATION | lag(LE, %d:%d)", lo, hi))
  m <- tryCatch(pgmm(fg, data = pd, effect = "twoways", model = "twosteps",
                     transformation = transformation, collapse = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- coef(m); V <- pgmm_vcov(m)
  s <- summary(m, robust = TRUE)
  ar <- sapply(1:4, function(o)
    tryCatch(plm::mtest(m, order = o)$p.value, error = function(e) NA_real_))
  sg <- tryCatch(s$sargan$p.value, error = function(e) NA_real_)
  sg_stat <- tryCatch(unname(s$sargan$statistic), error = function(e) NA_real_)
  sg_df   <- tryCatch(unname(s$sargan$parameter), error = function(e) NA_real_)
  a  <- unname(b["lag(LE, 1)"]); sa <- sqrt(V["lag(LE, 1)", "lag(LE, 1)"])
  bp <- unname(b[FOC]);          sb <- sqrt(V[FOC, FOC])
  lr <- longrun(b[c("lag(LE, 1)", FOC)],
                V[c("lag(LE, 1)", FOC), c("lag(LE, 1)", FOC)], "lag(LE, 1)")
  data.frame(
    depth = sprintf("t-%d to t-%d", lo, hi),
    lo = lo, hi = hi,
    n_inst = length(m$W[[1]][1, ]),
    alpha = a, se_alpha = sa, t_1ma = (1 - a) / sa,
    beta = bp, se_beta = sb, p_beta = 2 * pnorm(-abs(bp / sb)),
    imm = unname(b["IMM_DPT"]), se_imm = sqrt(V["IMM_DPT", "IMM_DPT"]),
    p_imm = 2 * pnorm(-abs(unname(b["IMM_DPT"]) / sqrt(V["IMM_DPT", "IMM_DPT"]))),
    ar1 = ar[1], ar2 = ar[2], ar3 = ar[3], ar4 = ar[4],
    hansen_chi2 = sg_stat, hansen_df = sg_df, hansen_p = sg,
    beta_LR = lr$beta_LR, fieller = lr$fieller_type,
    fieller_lo = lr$fieller_lo, fieller_hi = lr$fieller_hi,
    row.names = NULL)
}

cat("=== System GMM: instrument depth and the AR(3) rejection ===\n\n")
tab <- bind_rows(lapply(list(c(2,5), c(3,6), c(4,7), c(5,8), c(4,9)),
                        function(z) fit_depth(z[1], z[2])))
print(tab %>% dplyr::select(depth, n_inst, alpha, se_alpha, beta, se_beta,
                            p_beta, imm, p_imm, ar2, ar3, ar4, hansen_p) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

cat("\n=== Inference on the long run, by instrument depth ===\n")
print(tab %>% dplyr::select(depth, alpha, se_alpha, t_1ma, beta_LR, fieller,
                            fieller_lo, fieller_hi) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# ---- difference GMM at the same depths ---------------------
cat("\n=== The same sequence in difference GMM ===\n")
tabd <- bind_rows(lapply(list(c(2,5), c(3,6), c(4,7), c(5,8)),
                         function(z) fit_depth(z[1], z[2], "d")))
print(tabd %>% dplyr::select(depth, n_inst, alpha, beta, se_beta, p_beta,
                             ar2, ar3, ar4, hansen_p) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# ---- a second lag of the outcome ---------------------------
cat("\n=== Alternative: a second lag of the outcome in the model ===\n")
fg2 <- LE ~ lag(LE,1) + lag(LE,2) + log_GNI + EDU_GROSS + PUB + OOP_share +
            IMM_DPT + SANITATION | lag(LE, 3:6)
m2 <- tryCatch(pgmm(fg2, data=pd, effect="twoways", model="twosteps",
                    transformation="ld", collapse=TRUE), error=function(e) NULL)
if (!is.null(m2)) {
  b <- coef(m2); V <- pgmm_vcov(m2)
  a1 <- unname(b["lag(LE, 1)"]); a2 <- unname(b["lag(LE, 2)"])
  # long-run denominator for an AR(2): 1 - a1 - a2
  g <- c(-1, -1); idx <- c("lag(LE, 1)", "lag(LE, 2)")
  se_den <- sqrt(as.numeric(t(g) %*% V[idx, idx] %*% g))
  bp <- unname(b[FOC]); sb <- sqrt(V[FOC,FOC])
  ar <- sapply(1:4, function(o)
    tryCatch(plm::mtest(m2, order=o)$p.value, error=function(e) NA_real_))
  cat(sprintf("alpha1 = %.4f, alpha2 = %+.4f, sum = %.4f\n", a1, a2, a1+a2))
  cat(sprintf("1 - alpha1 - alpha2 = %.4f (SE %.4f, t = %.2f)\n",
              1-a1-a2, se_den, (1-a1-a2)/se_den))
  cat(sprintf("beta(PUB) = %+.5f (SE %.5f, p = %.3f)\n", bp, sb,
              2*pnorm(-abs(bp/sb))))
  cat(sprintf("long-run multiplier = %+.4f\n", bp/(1-a1-a2)))
  cat(sprintf("AR(1) %.3f  AR(2) %.3f  AR(3) %.3f  AR(4) %.3f\n",
              ar[1], ar[2], ar[3], ar[4]))
}

# ---- the ILR specification ---------------------------------
cat("\n=== The ILR specification at the two depths ===\n")
for (z in list(c(2,5), c(4,7))) {
  fgi <- as.formula(sprintf(
    "LE ~ lag(LE,1) + log_GNI + EDU_GROSS + ILR1 + ILR2 + ILR3 + IMM_DPT +
     SANITATION | lag(LE, %d:%d)", z[1], z[2]))
  mi <- tryCatch(pgmm(fgi, data=pd, effect="twoways", model="twosteps",
                      transformation="ld", collapse=TRUE), error=function(e) NULL)
  if (is.null(mi)) next
  b <- coef(mi); V <- pgmm_vcov(mi)
  ar <- sapply(1:4, function(o)
    tryCatch(plm::mtest(mi, order=o)$p.value, error=function(e) NA_real_))
  cat(sprintf("t-%d to t-%d: alpha = %.4f, ILR1 = %+.5f (SE %.5f), ",
              z[1], z[2], unname(b["lag(LE, 1)"]), unname(b["ILR1"]),
              sqrt(V["ILR1","ILR1"])))
  cat(sprintf("AR(2) %.3f AR(3) %.3f AR(4) %.3f\n", ar[2], ar[3], ar[4]))
}

# ---- Sources of the AR(3) rejection ----------------------
ar_p <- function(m) sapply(1:4, function(o)
  tryCatch(plm::mtest(m, order = o)$p.value, error = function(e) NA_real_))

cat("\n=== Pure AR(1) in life expectancy, no covariates ===\n")
m0 <- pgmm(LE ~ lag(LE,1) | lag(LE,2:5), data = pd, effect = "twoways",
           model = "twosteps", transformation = "ld", collapse = TRUE)
a0 <- ar_p(m0)
cat(sprintf("alpha = %.4f, ", unname(coef(m0)["lag(LE, 1)"])))
cat(sprintf("AR(2) %.3f AR(3) %.3f AR(4) %.3f\n", a0[2], a0[3], a0[4]))

cat("\n=== Diagnostics by end year of the window ===\n")
fgw <- LE ~ lag(LE,1) + log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
            SANITATION | lag(LE, 2:5)
wtab <- bind_rows(lapply(2018:2024, function(ymax) {
  d <- build_panel(year_max = ymax, zero_rule = "multiplicative")
  pw <- pdata.frame(d, index = c("iso3c","year"))
  m <- tryCatch(pgmm(fgw, data = pw, effect = "twoways", model = "twosteps",
                     transformation = "ld", collapse = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- coef(m); V <- pgmm_vcov(m); a <- ar_p(m)
  data.frame(window = sprintf("2000-%d", ymax),
             alpha = unname(b["lag(LE, 1)"]),
             ar1 = a[1], ar2 = a[2], ar3 = a[3], ar4 = a[4],
             beta = unname(b[FOC]),
             p_beta = 2*pnorm(-abs(b[FOC]/sqrt(V[FOC,FOC]))), row.names = NULL)
}))
print(wtab %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

cat("\n=== Autocorrelation of the within-country residuals ===\n")
mw <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION,
          data = pd, model = "within", effect = "twoways")
rr <- residuals(mw); ix <- attr(rr, "index")
rd <- data.frame(i = ix[[1]], t = as.integer(as.character(ix[[2]])),
                 e = as.numeric(rr)) %>% arrange(i, t) %>%
      group_by(i) %>% mutate(d = e - dplyr::lag(e)) %>% ungroup()
ac <- numeric(4)
for (kk in 1:4) {
  lk <- rd %>% group_by(i) %>% mutate(dl = dplyr::lag(d, kk)) %>% ungroup()
  ac[kk] <- cor(lk$d, lk$dl, use = "complete.obs")
  cat(sprintf("  differenced residuals, lag %d: %+.3f\n", kk, ac[kk]))
}
cat(sprintf("largest in absolute value beyond lag 1: lag %d\n",
            1 + which.max(abs(ac[-1]))))

write.csv(wtab, file.path(OUT, "31_ar3_windows.csv"), row.names = FALSE)
write.csv(tab,  file.path(OUT, "31_ar3_depth_sysgmm.csv"), row.names = FALSE)
write.csv(tabd, file.path(OUT, "31_ar3_depth_diffgmm.csv"), row.names = FALSE)
cat("\nSaved: 31_ar3_windows.csv, 31_ar3_depth_sysgmm.csv,\n")
cat("31_ar3_depth_diffgmm.csv\n")

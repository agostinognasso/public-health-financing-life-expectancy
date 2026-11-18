# 20: the estimator ladder on 2000-2024 and on 2000-2019
#
# The coefficient on the public share under (i) static FE with
# country and year effects, (ii) LSDV with the lagged outcome,
# (iii) LSDV corrected by the half-panel jackknife of Dhaene and
# Jochmans (2015), (iv) difference GMM and (v) System GMM. With
# T = 25 the Nickell bias of LSDV is O(1/T), so the intermediate
# rungs separate the effect of adding the lagged outcome from the
# effect of instrumenting. Every rung is fitted on the full window
# and on 2000-2019 (Section 6.1), and the distributed lag of
# Section 5.2 on both windows.
#
# Reads  data/ through build_panel()
# Writes output/20_ladder.csv, 20_dlag_windows.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

FOC <- "PUB"   # focal regressor

fit_ladder <- function(dat, label) {
  pd <- pdata.frame(dat, index = c("iso3c", "year"))
  X  <- c("log_GNI", "EDU_GROSS", "PUB", "OOP_share", "IMM_DPT", "SANITATION")
  rhs <- paste(X, collapse = " + ")

  out <- list()

  # (i) static FE, two-way, country-clustered
  m1 <- plm(as.formula(paste("LE ~", rhs)), data = pd,
            model = "within", effect = "twoways")
  V1 <- plm::vcovHC(m1, cluster = "group", method = "arellano")
  out$fe <- data.frame(rung = "1 static FE (two-way)",
                       alpha = NA_real_,
                       beta = coef(m1)[FOC], se = sqrt(V1[FOC, FOC]))

  # (ii) LSDV: FE with the lagged dependent variable
  m2 <- plm(as.formula(paste("LE ~ lag(LE, 1) +", rhs)), data = pd,
            model = "within", effect = "twoways")
  V2 <- plm::vcovHC(m2, cluster = "group", method = "arellano")
  a2 <- coef(m2)["lag(LE, 1)"]
  out$lsdv <- data.frame(rung = "2 LSDV (lagged outcome)",
                         alpha = a2,
                         beta = coef(m2)[FOC], se = sqrt(V2[FOC, FOC]))

  # (iii) half-panel jackknife bias correction (Dhaene-Jochmans 2015)
  yrs <- sort(unique(dat$year))
  cut <- yrs[floor(length(yrs) / 2)]
  h1 <- dat %>% filter(year <= cut)
  h2 <- dat %>% filter(year >  cut)
  fit_half <- function(d) {
    p <- pdata.frame(d, index = c("iso3c", "year"))
    tryCatch(coef(plm(as.formula(paste("LE ~ lag(LE, 1) +", rhs)),
                      data = p, model = "within", effect = "twoways")),
             error = function(e) NULL)
  }
  c1 <- fit_half(h1); c2 <- fit_half(h2)
  if (!is.null(c1) && !is.null(c2)) {
    cf  <- coef(m2)
    nm  <- intersect(names(cf), intersect(names(c1), names(c2)))
    hpj <- 2 * cf[nm] - (c1[nm] + c2[nm]) / 2
    out$lsdvbc <- data.frame(rung = "3 LSDV bias-corrected (HPJ)",
                             alpha = hpj["lag(LE, 1)"],
                             beta = hpj[FOC], se = NA_real_)
  }

  # (iv) difference GMM, (v) System GMM
  fg <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
             IMM_DPT + SANITATION | lag(LE, 2:5)
  m4 <- tryCatch(pgmm(fg, data = pd, effect = "twoways", model = "twosteps",
                      transformation = "d", collapse = TRUE),
                 error = function(e) NULL)
  m5 <- tryCatch(pgmm(fg, data = pd, effect = "twoways", model = "twosteps",
                      transformation = "ld", collapse = TRUE),
                 error = function(e) NULL)
  for (nmm in c("dgmm", "sgmm")) {
    m <- if (nmm == "dgmm") m4 else m5
    if (is.null(m)) next
    b <- coef(m); V <- pgmm_vcov(m)
    out[[nmm]] <- data.frame(
      rung = ifelse(nmm == "dgmm", "4 difference GMM", "5 System GMM"),
      alpha = b["lag(LE, 1)"], beta = b[FOC], se = sqrt(V[FOC, FOC]))
  }

  tab <- bind_rows(out) %>%
    mutate(window = label,
           z = beta / se, p = 2 * pnorm(-abs(z)),
           lo = beta - 1.96 * se, hi = beta + 1.96 * se,
           beta_LR = ifelse(is.na(alpha), beta, beta / (1 - alpha)),
           .before = 1) %>%
    relocate(window, rung)
  rownames(tab) <- NULL
  attr(tab, "models") <- list(fe = m1, lsdv = m2, dgmm = m4, sgmm = m5)
  tab
}

# ---- run both windows --------------------------------------
d_full <- build_panel(year_max = 2024, zero_rule = "multiplicative")
d_pre  <- build_panel(year_max = 2019, zero_rule = "multiplicative")

cat("=== Samples ===\n")
cat(sprintf("2000-2024: %d obs, %d countries\n", nrow(d_full),
            n_distinct(d_full$iso3c)))
cat(sprintf("2000-2019: %d obs, %d countries\n", nrow(d_pre),
            n_distinct(d_pre$iso3c)))

L_full <- fit_ladder(d_full, "2000-2024")
L_pre  <- fit_ladder(d_pre,  "2000-2019")

cat("\n=== The ladder: coefficient on the public share ===\n")
lad <- bind_rows(L_full, L_pre)
print(lad %>% mutate(across(c(alpha, beta, se, lo, hi, beta_LR),
                            ~ round(.x, 4)),
                     z = round(z, 2), p = round(p, 4)) %>%
        dplyr::select(window, rung, alpha, beta, se, p, lo, hi, beta_LR))

cat("\n=== Step by step, full window ===\n")
g <- function(tab, r) tab$beta[grepl(r, tab$rung)]
b_fe   <- g(L_full, "^1"); b_lsdv <- g(L_full, "^2")
b_bc   <- g(L_full, "^3"); b_dg   <- g(L_full, "^4"); b_sg <- g(L_full, "^5")
cat(sprintf("static FE               %+.4f\n", b_fe))
cat(sprintf("+ lagged outcome (LSDV) %+.4f   change from the lagged outcome: %+.4f\n",
            b_lsdv, b_lsdv - b_fe))
cat(sprintf("+ HPJ bias correction   %+.4f   change from the correction:     %+.4f\n",
            b_bc, b_bc - b_lsdv))
cat(sprintf("difference GMM          %+.4f\n", b_dg))
cat(sprintf("System GMM              %+.4f   change from HPJ to System GMM:  %+.4f\n",
            b_sg, b_sg - b_bc))
cat(sprintf("\nTotal move from static FE to System GMM: %+.4f to %+.4f (%+.4f);\n",
            b_fe, b_sg, b_sg - b_fe))
cat(sprintf("the lagged-outcome step is %.0f%% of it, the HPJ-to-System-GMM step %.0f%%.\n",
            100 * (b_lsdv - b_fe) / (b_sg - b_fe),
            100 * (b_sg - b_bc) / (b_sg - b_fe)))

cat("\n=== 2000-2024 against 2000-2019 ===\n")
cmp <- lad %>%
  dplyr::select(rung, window, beta, se, p) %>%
  tidyr::pivot_wider(names_from = window, values_from = c(beta, se, p))
print(cmp %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

sg_pre <- L_pre %>% filter(grepl("^5", rung))
cat(sprintf("\nSystem GMM on 2000-2019: beta(PUB) = %+.4f (SE %.4f, p = %.3f)\n",
            sg_pre$beta, sg_pre$se, sg_pre$p))
cat(sprintf("alpha = %.4f, long-run multiplier = %+.4f\n",
            sg_pre$alpha, sg_pre$beta_LR))
cat(sprintf("p > 0.05 on 2000-2019: %s\n", ifelse(sg_pre$p > 0.05, "yes", "no")))

# ---- distributed lag cumulative effect, both windows --------
cat("\n=== Distributed lag cumulative effect on the public share ===\n")
dl <- lapply(list(`2000-2024` = d_full, `2000-2019` = d_pre), function(d) {
  p <- pdata.frame(d, index = c("iso3c", "year"))
  m <- plm(LE ~ lag(PUB, 0:3) + lag(log_GNI, 0:3) + lag(IMM_DPT, 0:3) +
                lag(SANITATION, 0:3) + lag(EDU_GROSS, 0:3) +
                lag(OOP_share, 0:3),
           data = p, model = "within", effect = "twoways")
  V <- plm::vcovHC(m, cluster = "group", method = "arellano")
  idx <- grep("^lag\\(PUB", names(coef(m)))
  R <- rep(0, length(coef(m))); R[idx] <- 1
  cum <- sum(coef(m)[idx]); se <- sqrt(as.numeric(t(R) %*% V %*% R))
  data.frame(cum = cum, se = se, z = cum / se,
             p = 2 * pnorm(-abs(cum / se)),
             lo = cum - 1.96 * se, hi = cum + 1.96 * se, n_lags = length(idx))
})
dlt <- bind_rows(dl, .id = "window")
print(dlt %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

write.csv(lad, file.path(OUT, "20_ladder.csv"), row.names = FALSE)
write.csv(dlt, file.path(OUT, "20_dlag_windows.csv"), row.names = FALSE)
cat("\nSaved: 20_ladder.csv, 20_dlag_windows.csv\n")

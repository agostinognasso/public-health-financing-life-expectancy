# 18: zeros in the financing simplex
#
# The log-ratio coordinates of equations (3)-(5) are undefined at
# zero, so zeros are replaced before the transform. The script counts
# how often the clipping of equations (1) and (2) binds, counts the
# zeros and the parts below the detection limit, relates their
# presence to income and region, and re-estimates the ILR System GMM
# under five detection limits, two replacement rules and the
# Bayesian-multiplicative rule of zCompositions (Section 3.2, Supplementary
# Section S5, Section 6.2 and Supplementary Figure S11).
#
# Reads  data/df_panel_drivers.rds, and data/ through build_panel()
# Writes output/18_zero_audit.csv, 18_delta_sensitivity.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(zCompositions))

# ---- 1. the clipping in equations (1) and (2) --------------
raw <- readRDS(file.path(DATA, "df_panel_drivers.rds"))
fin <- raw %>% filter(!is.na(PUB), !is.na(PRIV), !is.na(OOP_share))

n_eq2 <- sum(fin$PRIV - fin$OOP_share < 0)
n_eq4 <- sum(100 - fin$PUB - fin$PRIV < 0)
n_any <- sum((fin$PRIV - fin$OOP_share < 0) | (100 - fin$PUB - fin$PRIV < 0))

cat("=== Clipping in equations (1) and (2) ===\n")
cat(sprintf("records with the financing block observed : %d\n", nrow(fin)))
cat(sprintf("eq. (1) PNO = max(PRIV-OOP,0) binds       : %d (%.2f%%), max shortfall %.3f pp\n",
            n_eq2, 100 * n_eq2 / nrow(fin), max(fin$OOP_share - fin$PRIV)))
cat(sprintf("eq. (2) EXT = max(100-PUB-PRIV,0) binds   : %d (%.2f%%), max overshoot %.3f pp\n",
            n_eq4, 100 * n_eq4 / nrow(fin), max(fin$PUB + fin$PRIV - 100)))
cat(sprintf("either clipping binds                    : %d (%.2f%%)\n",
            n_any, 100 * n_any / nrow(fin)))
cat("\nThe base for these percentages is the financing records of the raw\n")
cat("panel, not the observations of the analytical panel built below, and\n")
cat("the two equations are counted separately.\n")

# ---- 2. zeros in the analytical panel ----------------------
d0 <- build_panel(year_max = 2024, zero_rule = "none")
parts <- c("c_PUB", "c_OOP", "c_PNO", "c_EXT")
lab   <- c("Public (PUB)", "Out-of-pocket (OOP)",
           "Private non-OOP (PNO)", "External (EXT)")

cat("\n\n=== Zeros in the analytical panel ===\n")
zt <- data.frame(component = lab,
                 n_zero = sapply(parts, function(v) sum(d0[[v]] <= 0)),
                 pct    = round(100 * sapply(parts,
                            function(v) mean(d0[[v]] <= 0)), 2),
                 min_positive = signif(sapply(parts,
                            function(v) min(d0[[v]][d0[[v]] > 0])), 3),
                 n_below_delta = sapply(parts, function(v) sum(d0[[v]] < 0.05)),
                 n_residue = sapply(parts, function(v) sum(d0[[v]] > 0 & d0[[v]] < 1e-9)),
                 row.names = NULL)
print(zt)
cat("n_below_delta: parts below the detection limit delta = 0.05 pp,\n")
cat("zeros included; n_residue: positive parts below 1e-9, which are\n")
cat("floating-point residues of the subtraction in equation (2).\n")
cat(sprintf("\nobservations with at least one zero: %d (%.2f%%)\n",
            sum(d0$n_zero > 0), 100 * mean(d0$n_zero > 0)))
cat(sprintf("countries with at least one zero   : %d of %d\n",
            n_distinct(d0$iso3c[d0$n_zero > 0]), n_distinct(d0$iso3c)))
cat(sprintf("observations with a part below delta = 0.05 (zeros included): %d (%.2f%%)\n",
            sum(d0$n_below > 0), 100 * mean(d0$n_below > 0)))
cat(sprintf("countries with a part below delta   : %d of %d\n",
            n_distinct(d0$iso3c[d0$n_below > 0]), n_distinct(d0$iso3c)))
cat(sprintf("countries with EXT zero in every year: %d\n",
            d0 %>% group_by(iso3c) %>%
              summarise(all0 = all(c_EXT <= 0), .groups = "drop") %>%
              summarise(n = sum(all0)) %>% pull(n)))

# ---- 3. where the zeros are --------------------------------
cat("\n=== Zeros by income and region ===\n")
sel <- d0 %>% mutate(has_zero = as.integer(n_zero > 0))
cat("mean log income per head, by presence of a zero component:\n")
print(sel %>% group_by(has_zero) %>%
        summarise(n = n(), log_GNI = round(mean(log_GNI), 3),
                  LE = round(mean(LE), 2), .groups = "drop"))
cat("\nshare of observations with a zero, by region:\n")
print(sel %>% group_by(region) %>%
        summarise(n = n(), pct_zero = round(100 * mean(has_zero), 1),
                  .groups = "drop") %>% arrange(desc(pct_zero)), n = 10)
cat("\nlogit of has_zero on log income:\n")
gl <- glm(has_zero ~ log_GNI, data = sel, family = binomial)
print(round(summary(gl)$coefficients, 4))

cat("\nThe same for the parts that are replaced, those below delta:\n")
selb <- d0 %>% mutate(has_below = as.integer(n_below > 0))
print(selb %>% group_by(has_below) %>%
        summarise(n = n(), log_GNI = round(mean(log_GNI), 3),
                  LE = round(mean(LE), 2), .groups = "drop"))
print(selb %>% group_by(region) %>%
        summarise(n = n(), pct_below = round(100 * mean(has_below), 1),
                  .groups = "drop") %>% arrange(desc(pct_below)), n = 10)

# ---- 4. sensitivity of the ILR1 coefficient to delta -------
cat("\n\n=== Sensitivity of the ILR System GMM to the zero treatment ===\n")
fB <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + ILR1 + ILR2 + ILR3 +
           IMM_DPT + SANITATION | lag(LE, 2:5)

# alpha moves with the zero treatment, and with it the confidence
# level up to which a Fieller set for a long-run multiplier is an
# interval; reported for the ILR1 multiplier at each treatment
fieller_cols <- function(b, V) {
  k   <- c("lag(LE, 1)", "ILR1")
  l95 <- longrun(b[k], V[k, k], alpha_name = "lag(LE, 1)")
  l90 <- longrun(b[k], V[k, k], alpha_name = "lag(LE, 1)", level = 0.90)
  lv  <- fieller_levels(b[k], V[k, k], "lag(LE, 1)")
  data.frame(se_alpha = sqrt(V["lag(LE, 1)", "lag(LE, 1)"]),
             level_bounded = lv$level_bounded,
             ILR1_fieller_type = l95$fieller_type,
             ILR1_fieller90_lo = l90$fieller_lo,
             ILR1_fieller90_hi = l90$fieller_hi,
             ILR1_fieller90_type = l90$fieller_type,
             ILR1_level_wholeline = lv$level_wholeline,
             row.names = NULL)
}

grid <- expand.grid(delta = c(0.01, 0.025, 0.05, 0.10, 0.25),
                    rule  = c("multiplicative", "additive"),
                    stringsAsFactors = FALSE)

res <- lapply(seq_len(nrow(grid)), function(i) {
  dd <- build_panel(year_max = 2024, delta = grid$delta[i],
                    zero_rule = grid$rule[i])
  pp <- pdata.frame(dd, index = c("iso3c", "year"))
  m  <- tryCatch(pgmm(fB, data = pp, effect = "twoways", model = "twosteps",
                      transformation = "ld", collapse = TRUE),
                 error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- coef(m); V <- pgmm_vcov(m)
  a <- b["lag(LE, 1)"]
  data.frame(rule = grid$rule[i], delta = grid$delta[i],
             alpha = a,
             ILR1 = b["ILR1"], se_ILR1 = sqrt(V["ILR1", "ILR1"]),
             ILR1_LR = b["ILR1"] / (1 - a),
             IMM = b["IMM_DPT"], se_IMM = sqrt(V["IMM_DPT", "IMM_DPT"]),
             fieller_cols(b, V),
             row.names = NULL)
})
sens <- bind_rows(res) %>%
  mutate(p_ILR1 = 2 * pnorm(-abs(ILR1 / se_ILR1)),
         p_IMM  = 2 * pnorm(-abs(IMM / se_IMM))) %>%
  relocate(p_ILR1, p_IMM, .after = se_IMM)
print(sens %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

# Bayesian-multiplicative replacement (zCompositions)
cat("\n--- Bayesian-multiplicative replacement, zCompositions::cmultRepl ---\n")
dz <- build_panel(year_max = 2024, zero_rule = "none")
cm <- as.matrix(dz[, parts]); cm[is.na(cm)] <- 0
cm[cm < 0.05] <- 0   # below the detection limit, as in build_panel()
bm <- tryCatch({
  rep <- cmultRepl(cm, label = 0, method = "GBM",
                   output = "prop", z.warning = 0.9, suppress.print = TRUE)
  rep <- rep / rowSums(rep) * 100
  colnames(rep) <- c("z_PUB", "z_OOP", "z_PNO", "z_EXT")
  dd <- dz; dd[, colnames(rep)] <- as.data.frame(rep)
  sbp <- matrix(c(+1,-1,-1,-1, 0,+1,+1,-1, 0,+1,-1,0), nrow = 3, byrow = TRUE,
                dimnames = list(c("ILR1","ILR2","ILR3"),
                                c("z_PUB","z_OOP","z_PNO","z_EXT")))
  ilr <- ilr_from_sbp(dd[, colnames(sbp)], sbp)
  dd$ILR1 <- ilr[, 1]; dd$ILR2 <- ilr[, 2]; dd$ILR3 <- ilr[, 3]
  pp <- pdata.frame(dd, index = c("iso3c", "year"))
  m <- pgmm(fB, data = pp, effect = "twoways", model = "twosteps",
            transformation = "ld", collapse = TRUE)
  b <- coef(m); V <- pgmm_vcov(m)
  data.frame(rule = "GBM (cmultRepl)", delta = NA,
             alpha = b["lag(LE, 1)"],
             ILR1 = b["ILR1"], se_ILR1 = sqrt(V["ILR1","ILR1"]),
             ILR1_LR = b["ILR1"] / (1 - b["lag(LE, 1)"]),
             IMM = b["IMM_DPT"], se_IMM = sqrt(V["IMM_DPT","IMM_DPT"]),
             fieller_cols(b, V),
             row.names = NULL)
}, error = function(e) { cat("cmultRepl failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(bm)) {
  bm <- bm %>% mutate(p_ILR1 = 2*pnorm(-abs(ILR1/se_ILR1)),
                      p_IMM = 2*pnorm(-abs(IMM/se_IMM))) %>%
    relocate(p_ILR1, p_IMM, .after = se_IMM)
  print(bm %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
  sens <- bind_rows(sens, bm)
}

cat("\nRange of ILR1 across all zero treatments: [",
    sprintf("%+.4f, %+.4f", min(sens$ILR1), max(sens$ILR1)), "]\n")
cat("ILR1 significant at 5% under any treatment:",
    ifelse(any(sens$p_ILR1 < 0.05), "yes", "no"), "\n")
cat("immunisation significant at 5% under every treatment:",
    ifelse(all(sens$p_IMM < 0.05), "yes", "no"),
    sprintf("(max p = %.3f)\n", max(sens$p_IMM)))

cat("\n=== Persistence and the Fieller sets across zero treatments ===\n")
print(sens %>%
        transmute(rule, delta, alpha = round(alpha, 4),
                  se_alpha = round(se_alpha, 4),
                  bounded_up_to = sprintf("%.1f%%", 100 * level_bounded),
                  ILR1_LR = round(ILR1_LR, 4),
                  ILR1_fieller90 = mapply(fmt_fieller, ILR1_fieller90_lo,
                                          ILR1_fieller90_hi, ILR1_fieller90_type),
                  ILR1_fieller95 = ILR1_fieller_type), row.names = FALSE)
cat(sprintf("\nalpha runs from %.3f to %.3f; Fieller sets are intervals up to a\n",
            min(sens$alpha), max(sens$alpha)))
cat(sprintf("level between %.1f%% and %.1f%%, depending on the zero treatment.\n",
            100 * min(sens$level_bounded), 100 * max(sens$level_bounded)))

write.csv(zt,   file.path(OUT, "18_zero_audit.csv"),  row.names = FALSE)
write.csv(sens, file.path(OUT, "18_delta_sensitivity.csv"), row.names = FALSE)
cat("\nSaved: 18_zero_audit.csv, 18_delta_sensitivity.csv\n")

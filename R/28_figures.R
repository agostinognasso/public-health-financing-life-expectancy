# 28: Figures 1 and 2, Supplementary Figures S7 to S12
#
#   Figure 1    mixed-model coefficients, between and within, with
#               95 per cent intervals
#   Figure 2    the estimator ladder, impact coefficient and long-run
#               multiplier
#   Figure S7   paired VIP bootstrap under three resampling schemes
#   Figure S8   Bayesian fixed effects (main fit of script 24)
#   Figure S9   region-specific within-country slopes (script 24)
#   Figure S10  leave-one-region-out, decomposed
#   Figure S11  sensitivity to the zero-replacement rule
#   Figure S12  minimum detectable effects and the intervals for the
#               public share
#
# Reads  output/16_sysgmm_longrun.csv, 18_delta_sensitivity.csv,
#        19_mde.csv, 19_mde_longrun.csv, 23_loro_decomposed.csv,
#        23_vip_bootstrap.csv, 24_fixed_effects.csv,
#        24_region_slopes.csv, 26_lsdvc.csv, 26_lsdvc_bootstrap.csv;
#        data/ through build_panel()
# Writes figures/Figure1, Figure2, FigureS7 to FigureS12;
#        output/28_ladder_impact.csv, 28_ladder_longrun.csv,
#        28_mixed_coefficients.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages({
  library(ggplot2); library(forcats); library(patchwork)
})

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3, colour = "grey92"),
          axis.line.x = element_line(linewidth = 0.4, colour = "black"),
          strip.text = element_text(face = "bold", size = 12, hjust = 0),
          legend.position = "none")
}
# No legends, titles or annotations inside the figures: the captions
# describe the symbols, panels carry only their letter, and variable
# codes are replaced by their names.
DRIVER_NAMES <- c(log_GNI = "Log income per head", EDU = "Secondary enrolment",
                  IMM = "Immunisation", SAN = "Basic sanitation",
                  PUB = "Public share", OOP = "Out-of-pocket share",
                  PRIV_NONOOP = "Private non-out-of-pocket share",
                  EXT = "External share", HIV = "HIV prevalence",
                  ILR1 = "Balance: public vs non-public",
                  ILR2 = "Balance: domestic private vs external",
                  ILR3 = "Balance: out-of-pocket vs prepaid")

dat <- build_panel(zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c", "year"))
rhs <- "log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION"

# Segments that draw a Fieller set inside [-L, L]. Each piece is
# split at its midpoint and drawn outwards, so that an arrowhead can
# mark an end that runs on to infinity or past the edge of the panel.
ARROW <- arrow(length = unit(0.07, "inches"), type = "closed")
fieller_segments <- function(y, lo, hi, type, L) {
  half <- function(y, a, b, open_a, open_b) {
    m <- (a + b) / 2
    rbind(data.frame(y = y, x = m, xend = a, arrow = open_a),
          data.frame(y = y, x = m, xend = b, arrow = open_b))
  }
  edge <- function(y, side)                # a set lying wholly off-panel
    data.frame(y = y, x = side * 0.9 * L, xend = side * L, arrow = TRUE)
  bind_rows(lapply(seq_along(y), function(i) switch(type[i],
    "bounded"    = half(y[i], max(lo[i], -L), min(hi[i], L),
                        lo[i] < -L, hi[i] > L),
    "whole line" = half(y[i], -L, L, TRUE, TRUE),
    "exclusive/unbounded" = rbind(
      if (lo[i] > -L) half(y[i], -L, min(lo[i], L), TRUE, lo[i] > L)
      else edge(y[i], -1),
      if (hi[i] <  L) half(y[i], max(hi[i], -L), L, hi[i] < -L, TRUE)
      else edge(y[i], 1)))))
}

# ------------------------------------------------------------
# the ladder, impact and long run
# ------------------------------------------------------------
rowify <- function(b, se, alpha, lab, kind) {
  lr <- if (is.na(alpha)) b else b / (1 - alpha)
  data.frame(rung = lab, kind = kind, beta = b, se = se, alpha = alpha,
             beta_LR = lr, row.names = NULL)
}

m1  <- plm(as.formula(paste("LE ~", rhs)), data = pd, model = "within")
m2  <- plm(as.formula(paste("LE ~", rhs)), data = pd, model = "within",
           effect = "twoways")
m3  <- plm(as.formula(paste("LE ~ lag(LE,1) +", rhs)), data = pd,
           model = "within", effect = "twoways")
V1 <- plm::vcovHC(m1, cluster = "group", method = "arellano")
V2 <- plm::vcovHC(m2, cluster = "group", method = "arellano")
V3 <- plm::vcovHC(m3, cluster = "group", method = "arellano")

fg  <- LE ~ lag(LE,1) + log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
            SANITATION | lag(LE, 2:5)
# every driver endogenous, instrumented from t-2: specification C of
# script 21, the bottom row of Table 2
fg2 <- LE ~ lag(LE,1) + log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT +
            SANITATION | lag(LE, 2:5) + lag(PUB, 2:4) + lag(OOP_share, 2:4) +
            lag(log_GNI, 2:4) + lag(EDU_GROSS, 2:4) + lag(IMM_DPT, 2:4) +
            lag(SANITATION, 2:4)
m4 <- pgmm(fg,  data = pd, effect = "twoways", model = "twosteps",
           transformation = "d",  collapse = TRUE)
m5 <- pgmm(fg,  data = pd, effect = "twoways", model = "twosteps",
           transformation = "ld", collapse = TRUE)
m6 <- pgmm(fg2, data = pd, effect = "twoways", model = "twosteps",
           transformation = "ld", collapse = TRUE)
V4 <- pgmm_vcov(m4); V5 <- pgmm_vcov(m5); V6 <- pgmm_vcov(m6)

lad <- bind_rows(
  rowify(coef(m1)["PUB"], sqrt(V1["PUB","PUB"]), NA,
         "1  FE, country effects only", "static"),
  rowify(coef(m2)["PUB"], sqrt(V2["PUB","PUB"]), NA,
         "2  FE, country and year effects", "static"),
  rowify(coef(m3)["PUB"], sqrt(V3["PUB","PUB"]), coef(m3)["lag(LE, 1)"],
         "3  LSDV, lagged outcome", "dynamic"),
  rowify(coef(m4)["PUB"], sqrt(V4["PUB","PUB"]), coef(m4)["lag(LE, 1)"],
         "5  Difference GMM", "dynamic"),
  rowify(coef(m5)["PUB"], sqrt(V5["PUB","PUB"]), coef(m5)["lag(LE, 1)"],
         "6  System GMM, drivers exogenous", "dynamic"),
  rowify(coef(m6)["PUB"], sqrt(V6["PUB","PUB"]), coef(m6)["lag(LE, 1)"],
         "7  System GMM, drivers endogenous", "dynamic")
)

# the Kiviet-Bruno rung, from script 26
lsd  <- read.csv(file.path(OUT, "26_lsdvc.csv"))
lsdb <- read.csv(file.path(OUT, "26_lsdvc_bootstrap.csv"))
names(lsdb) <- c("alpha_b", "beta_b")
lsd_row <- lsd[lsd$method == "LSDVC B1 (init BB)", ]
a_kb  <- lsd_row$alpha
b_kb  <- lsd_row$beta_PUB
sa_kb <- sd(lsdb$alpha_b)
sb_kb <- sd(lsdb$beta_b)
lad <- bind_rows(
  lad,
  rowify(b_kb, sb_kb, a_kb, "4  LSDV, Kiviet-Bruno correction", "dynamic")
) %>%
  mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se)

# display order, top to bottom; the first factor level sits at the bottom
RUNGS <- c("1  FE, country effects only",
           "2  FE, country and year effects",
           "3  LSDV, lagged outcome",
           "4  LSDV, Kiviet-Bruno correction",
           "5  Difference GMM",
           "6  System GMM, drivers exogenous",
           "7  System GMM, drivers endogenous")
# numeric y so the row order is unambiguous: the top row has the largest y
lad$ypos <- length(RUNGS) + 1L - match(lad$rung, RUNGS)
stopifnot(!any(is.na(lad$ypos)))

# long-run multipliers with delta-method intervals and Fieller type
lrtab <- bind_rows(Map(function(m, V, lab) {
    L <- longrun(coef(m), V, "lag(LE, 1)")
    cbind(rung = lab, L[L$term == "PUB", ])
  },
  list(m3, m4, m5, m6), list(V3, V4, V5, V6),
  RUNGS[c(3, 5, 6, 7)]))

# Kiviet-Bruno: long run and Fieller set from the bootstrap covariance
Vkb <- cov(cbind(lsdb$alpha_b, lsdb$beta_b))
kb <- longrun(c("lag(LE, 1)" = a_kb, PUB = b_kb),
              matrix(c(Vkb[1,1], Vkb[1,2], Vkb[2,1], Vkb[2,2]), 2, 2,
                     dimnames = list(c("lag(LE, 1)","PUB"),
                                     c("lag(LE, 1)","PUB"))),
              "lag(LE, 1)")
lrtab <- bind_rows(lrtab,
                   cbind(rung = "4  LSDV, Kiviet-Bruno correction",
                         kb[kb$term == "PUB", ]))
lrtab$ypos <- length(RUNGS) + 1L - match(lrtab$rung, RUNGS)
stopifnot(!any(is.na(lrtab$ypos)))
lrtab$alpha <- lad$alpha[match(lrtab$rung, lad$rung)]

# the axis spells out FE; the CSV output keeps the short rung names
RUNG_LABELS <- sub("FE,", "Fixed effects,", RUNGS, fixed = TRUE)
YSCALE <- function() scale_y_continuous(
  breaks = length(RUNGS):1, labels = RUNG_LABELS, limits = c(0.4, length(RUNGS) + 0.75))

cat("=== ladder, impact estimand ===\n")
print(lad %>% mutate(across(c(beta, se, lo, hi, alpha, beta_LR), ~ round(.x, 4))))
cat("\n=== ladder, long-run estimand ===\n")
print(lrtab %>% mutate(across(c(beta_SR, beta_LR, se_LR, lo_delta, hi_delta),
                              ~ round(.x, 4))) %>%
        dplyr::select(rung, ypos, beta_LR, se_LR, lo_delta, hi_delta, fieller_type))

PUZZLE <- coef(m1)[["PUB"]]

# ---- panel A: impact coefficients ---------------------------
pA <- ggplot(lad, aes(y = ypos, x = beta)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = PUZZLE, linetype = "dotted", colour = "grey35") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(aes(shape = kind, fill = kind), size = 2.6, colour = "black") +
  scale_shape_manual(values = c(static = 22, dynamic = 21), name = NULL) +
  scale_fill_manual(values = c(static = "grey85", dynamic = "grey30"), name = NULL) +
  YSCALE() +
  labs(x = "Years of life expectancy per percentage point", y = NULL) +
  coord_cartesian(xlim = c(-0.075, 0.055), clip = "off") +
  theme_pub(10) + theme(legend.position = "none",
                        plot.margin = margin(12, 10, 6, 6))

# ---- panel B: long-run multipliers, 95% Fieller sets --------
# Each set is drawn with its shape: an interval, two half-lines,
# or the whole real line. An arrow marks an end that runs on past
# the edge of the panel, to minus or plus infinity. A rung whose
# persistence estimate is one or more has no long-run multiplier to
# point to: its set is drawn and its point is not.
CAP <- 0.15
fsB <- fieller_segments(lrtab$ypos, lrtab$fieller_lo, lrtab$fieller_hi,
                        lrtab$fieller_type, CAP)
pB <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = PUZZLE, linetype = "dotted", colour = "grey35") +
  geom_segment(data = fsB %>% filter(!arrow),
               aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.5) +
  geom_segment(data = fsB %>% filter(arrow),
               aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.5,
               arrow = ARROW) +
  geom_point(data = lrtab %>% filter(alpha < 1), aes(x = beta_LR, y = ypos),
             size = 2.6, shape = 21, fill = "grey30", colour = "black") +
  YSCALE() +
  labs(x = "Years of life expectancy per percentage point", y = NULL) +
  coord_cartesian(xlim = c(-CAP, CAP)) +
  theme_pub(10) + theme(axis.text.y = element_blank(),
                        plot.margin = margin(12, 10, 6, 0))

F2 <- (pA | pB) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))
save_fig("Figure2", F2, 10.2, 4.4)

# ------------------------------------------------------------
# Figure 1: mixed-model coefficients with intervals
# ------------------------------------------------------------
suppressPackageStartupMessages(library(lme4))
dm <- dat %>% group_by(iso3c) %>%
  mutate(across(c(log_GNI, EDU_GROSS, IMM_DPT, SANITATION, ILR1, ILR2, ILR3),
                list(bw = ~ mean(.x), wn = ~ .x - mean(.x)))) %>%
  ungroup() %>% mutate(year_c = year - 2012)
fm <- LE ~ year_c + is_covid_year +
  log_GNI_bw + log_GNI_wn + EDU_GROSS_bw + EDU_GROSS_wn +
  IMM_DPT_bw + IMM_DPT_wn + SANITATION_bw + SANITATION_wn +
  ILR1_bw + ILR1_wn + ILR2_bw + ILR2_wn + ILR3_bw + ILR3_wn + (1 | iso3c)
mm <- lmer(fm, data = dm, REML = TRUE,
           control = lmerControl(calc.derivs = FALSE))
ss <- summary(mm)$coefficients
mix <- data.frame(term = rownames(ss), est = ss[, 1], se = ss[, 2],
                  row.names = NULL) %>%
  filter(grepl("_(bw|wn)$", term)) %>%
  mutate(variable = sub("_(bw|wn)$", "", term),
         dimension = ifelse(grepl("_bw$", term), "between countries",
                            "within countries"),
         lo = est - 1.96 * se, hi = est + 1.96 * se,
         variable = recode(variable,
           log_GNI = "Log income per head", EDU_GROSS = "Secondary enrolment",
           IMM_DPT = "Immunisation, DPT", SANITATION = "Basic sanitation",
           ILR1 = "Balance: public vs non-public",
           ILR2 = "Balance: domestic private vs external",
           ILR3 = "Balance: out-of-pocket vs prepaid"))
cat("\n=== mixed model coefficients on an interpretable scale ===\n")
print(mix %>% mutate(across(c(est, se, lo, hi), ~ round(.x, 4))) %>%
        dplyr::select(variable, dimension, est, se, lo, hi))

F1 <- ggplot(mix, aes(y = fct_rev(factor(variable)), x = est,
                      shape = dimension, fill = dimension)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, colour = "black",
             position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(`between countries` = 22,
                                `within countries` = 21), name = NULL) +
  scale_fill_manual(values = c(`between countries` = "grey85",
                               `within countries` = "grey30"), name = NULL) +
  labs(x = "Years of life expectancy per unit of the driver", y = NULL) +
  theme_pub(10)
save_fig("Figure1", F1, 8.6, 4.2)

# ------------------------------------------------------------
# Figure S11: zero-replacement sensitivity
# ------------------------------------------------------------
sens <- read.csv(file.path(OUT, "18_delta_sensitivity.csv"))
s1 <- sens %>% filter(!is.na(delta)) %>%
  mutate(lo = ILR1 - 1.96 * se_ILR1, hi = ILR1 + 1.96 * se_ILR1,
         rule = recode(rule, multiplicative = "multiplicative replacement",
                       additive = "additive replacement"))
gbm <- sens %>% filter(is.na(delta))
FS11 <- ggplot(s1, aes(x = factor(delta), y = ILR1, shape = rule, fill = rule)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45,
                position = position_dodge(width = 0.45)) +
  geom_point(size = 2.5, colour = "black",
             position = position_dodge(width = 0.45)) +
  { if (nrow(gbm)) geom_hline(yintercept = gbm$ILR1[1], linetype = "dashed",
                              colour = "grey40") } +
  scale_shape_manual(values = c(21, 22), name = NULL) +
  scale_fill_manual(values = c("grey30", "grey85"), name = NULL) +
  labs(x = expression(paste("Detection limit  ", delta,
                            "  (percentage points of current health expenditure)")),
       y = "Coefficient on the first log-ratio coordinate") +
  theme_pub(10)
save_fig("FigureS11", FS11, 8.4, 4.0)

# ------------------------------------------------------------
# Figure S10: leave-one-region-out, decomposed
# ------------------------------------------------------------
loro <- read.csv(file.path(OUT, "23_loro_decomposed.csv"))
l2 <- loro %>%
  dplyr::select(region, R2_total, R2_given_region_mean) %>%
  tidyr::pivot_longer(-region, names_to = "measure", values_to = "R2") %>%
  mutate(measure = recode(measure,
          R2_total = "as reported: benchmarked on the held-out region mean",
          R2_given_region_mean = "after granting the held-out region its mean"))
FS10 <- ggplot(l2, aes(y = fct_rev(factor(region)), x = R2,
                      shape = measure, fill = measure)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(size = 2.6, colour = "black") +
  scale_shape_manual(values = c(21, 22), name = NULL) +
  scale_fill_manual(values = c("grey30", "grey85"), name = NULL) +
  labs(x = expression(paste("Out-of-region ", R^2)), y = NULL) +
  theme_pub(10)
save_fig("FigureS10", FS10, 8.6, 4.0)

# ------------------------------------------------------------
# Figure S7: VIP bootstrap under three schemes
# ------------------------------------------------------------
bt <- read.csv(file.path(OUT, "23_vip_bootstrap.csv"))
bt$variable <- sub("_init$", "", bt$variable)
stopifnot(all(bt$variable %in% names(DRIVER_NAMES)))
bt$variable <- DRIVER_NAMES[bt$variable]
bt$panel <- substr(bt$scheme, 1, 1)          # schemes are labelled A, B, C
FS7 <- ggplot(bt, aes(y = fct_rev(factor(variable)), x = delta_vip)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45) +
  geom_point(size = 2.1, shape = 21, fill = "grey30", colour = "black") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Difference in variable importance, level minus shape", y = NULL) +
  theme_pub(9.5)
save_fig("FigureS7", FS7, 11.0, 3.6)

# ------------------------------------------------------------
# Figure S12: what the design can detect
# ------------------------------------------------------------
# Impact: Wald intervals. Long run: Fieller sets, and the MDE as the
# impact MDE over 1 - alpha (script 19), which has no upper bound
# over the 95 per cent interval for alpha. The dark band holds the
# effects too small to detect at any alpha in that interval, the
# light band those too small at the estimated alpha.
mdet <- read.csv(file.path(OUT, "19_mde.csv"))
mlr  <- read.csv(file.path(OUT, "19_mde_longrun.csv"))
plr  <- read.csv(file.path(OUT, "16_sysgmm_longrun.csv")) %>%
  filter(term == "PUB")
sri  <- mdet[mdet$horizon == "short run", ]
XL   <- 0.12
ivs <- data.frame(
  y     = c(2, 2, 1, 1),
  level = c("90", "95", "90", "95"),
  lo    = c(sri$beta - qnorm(0.95) * sri$se, sri$beta - qnorm(0.975) * sri$se,
            plr$fieller90_lo, plr$fieller_lo),
  hi    = c(sri$beta + qnorm(0.95) * sri$se, sri$beta + qnorm(0.975) * sri$se,
            plr$fieller90_hi, plr$fieller_hi),
  type  = c("bounded", "bounded", plr$fieller90_type, plr$fieller_type))
segS12 <- bind_rows(lapply(seq_len(nrow(ivs)), function(i)
  fieller_segments(ivs$y[i], ivs$lo[i], ivs$hi[i], ivs$type[i], XL) %>%
    mutate(level = ivs$level[i])))
band <- data.frame(
  y     = c(1, 2, 1),
  half  = c(mlr$mde80[mlr$alpha_at == "point estimate"], sri$mde80,
            mlr$mde80[mlr$alpha_at == "95% CI lower"]),
  shade = c("light", "dark", "dark"))
ptsS12 <- data.frame(y = c(2, 1), beta = c(sri$beta, plr$beta_LR))
LW <- c(`90` = 1.4, `95` = 0.45)

FS12 <- ggplot() +
  geom_rect(data = band, aes(xmin = -half, xmax = half,
                             ymin = y - 0.3, ymax = y + 0.3, fill = shade),
            colour = NA) +
  scale_fill_manual(values = c(light = "grey92", dark = "grey80"),
                    guide = "none") +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1) * abs(PUZZLE), linetype = "dotted",
             colour = "grey35") +
  geom_segment(data = segS12 %>% filter(!arrow),
               aes(x = x, xend = xend, y = y, yend = y, linewidth = level)) +
  geom_segment(data = segS12 %>% filter(arrow),
               aes(x = x, xend = xend, y = y, yend = y, linewidth = level),
               arrow = ARROW) +
  scale_linewidth_manual(values = LW, guide = "none") +
  geom_point(data = ptsS12, aes(x = beta, y = y), size = 2.8, shape = 21,
             fill = "grey25", colour = "black") +
  scale_y_continuous(breaks = c(2, 1),
                     labels = c("Impact (one year)", "Long run"),
                     limits = c(0.55, 2.45)) +
  labs(x = "Years of life expectancy per percentage point", y = NULL) +
  coord_cartesian(xlim = c(-XL, XL), clip = "off") +
  theme_pub(10) + theme(plot.margin = margin(8, 10, 6, 6))
save_fig("FigureS12", FS12, 8.8, 3.8)

# ------------------------------------------------------------
# Figure S8: Bayesian fixed effects, main fit of script 24
# ------------------------------------------------------------
bf <- read.csv(file.path(OUT, "24_fixed_effects.csv")) %>%
  filter(!term %in% c("year_c", "is_covid_year")) %>%
  mutate(var  = sub("_(bw|wn)_s$", "", term),
         type = ifelse(grepl("_bw_s$", term), "between", "within"))
ordS8 <- bf %>% group_by(var) %>% summarise(s = sum(abs(mean)), .groups = "drop") %>%
  arrange(desc(s)) %>% pull(var)
labS8 <- DRIVER_NAMES[c("log_GNI", "EDU", "IMM", "SAN", "HIV",
                        "ILR1", "ILR2", "ILR3")]
bf$var_lab <- factor(labS8[bf$var], levels = rev(labS8[ordS8]))
FS8 <- ggplot(bf, aes(x = mean, y = var_lab, shape = type, linetype = type)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = q025, xmax = q975), width = 0.2,
                position = position_dodge(0.6), colour = "black") +
  geom_point(size = 2.8, position = position_dodge(0.6),
             colour = "black", fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(between = 16, within = 21),
                     labels = c("Between (cross-country)", "Within (year-to-year)"),
                     name = NULL) +
  scale_linetype_manual(values = c(between = "solid", within = "22"),
                        labels = c("Between (cross-country)", "Within (year-to-year)"),
                        name = NULL) +
  labs(x = "Standardised posterior coefficient (95% credible interval)", y = NULL) +
  theme_pub(11)
save_fig("FigureS8", FS8, 9, 5.5)

# ------------------------------------------------------------
# Figure S9: region-specific within-country slopes, script 24
# ------------------------------------------------------------
rsS9 <- read.csv(file.path(OUT, "24_region_slopes.csv")) %>%
  mutate(driver = recode(term, IMM_wn_s = "Immunisation",
                         SAN_wn_s = "Sanitation", EDU_wn_s = "Education"),
         driver = factor(driver, levels = c("Immunisation", "Sanitation", "Education")))
rsS9$region <- fct_reorder(factor(rsS9$region), rsS9$mean, .fun = mean)
FS9 <- ggplot(rsS9, aes(x = mean, y = region, shape = driver, fill = driver)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, linewidth = 0.4,
                 position = position_dodge(0.7)) +
  geom_point(size = 2.6, colour = "black", position = position_dodge(0.7)) +
  scale_shape_manual(values = c(21, 22, 24), name = NULL) +
  scale_fill_manual(values = c("grey15", "grey55", "white"), name = NULL) +
  labs(x = "Within-country slope, global plus regional (standardised, 95% credible interval)",
       y = NULL) +
  theme_pub(11)
save_fig("FigureS9", FS9, 9.5, 5)

write.csv(lad,   file.path(OUT, "28_ladder_impact.csv"), row.names = FALSE)
write.csv(lrtab, file.path(OUT, "28_ladder_longrun.csv"), row.names = FALSE)
write.csv(mix,   file.path(OUT, "28_mixed_coefficients.csv"), row.names = FALSE)
cat("\nAll figures written to", FIG, "\n")

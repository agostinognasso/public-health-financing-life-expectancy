# 39: what the estimates imply for a planned shift in the public share
#
# Collects, for each set of assumptions reported in the article, the
# impact coefficient on the public share and translates it into the
# terms of a planned shift of ten percentage points (Section 7.2,
# Table 5):
#   - the smallest equivalence margin the estimate establishes at the
#     5 per cent level, |beta| + z(0.95) * SE, the margin at which the
#     two one-sided tests of script 19 would just reject;
#   - the 90 per cent interval for the change in life expectancy
#     within the year, 10 * (beta -/+ z(0.95) * SE);
#   - the 90 per cent Fieller set for the long-run change, 10 times
#     the set for the multiplier, which is undefined when the
#     estimated persistence is at least one.
# No model is fitted here; every number is read from the outputs of
# the scripts listed below. Effects are linear in the shift, so a
# shift of five or twenty points halves or doubles every entry.
#
# Reads  output/16_sysgmm_longrun.csv, 21_exogeneity.csv,
#        26_lsdvc.csv, 26_lsdvc_bootstrap.csv, 35_hiv.csv,
#        37_total_spending.csv, 38_income.csv
# Writes output/39_planning.csv, 39_planning_table.tex

source("00_common.R")

SHIFT <- 10
z90 <- qnorm(0.95)

rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE)
r16 <- rd("16_sysgmm_longrun.csv")
r21 <- rd("21_exogeneity.csv")
r26 <- rd("26_lsdvc.csv")
r26b <- rd("26_lsdvc_bootstrap.csv")
r35 <- rd("35_hiv.csv")
r37 <- rd("37_total_spending.csv")
r38 <- rd("38_income.csv")

row <- function(label, beta, se, alpha, f_lo, f_hi, f_type) {
  data.frame(spec = label, beta = beta, se = se, alpha = alpha,
             fieller90_lo = f_lo, fieller90_hi = f_hi,
             fieller90_type = f_type, stringsAsFactors = FALSE)
}

pub16 <- r16[r16$term == "PUB", ]
a21 <- function(s) r21[r21$spec == s & r21$term == "PUB", ]
p37 <- function(s) r37[r37$spec == s & r37$term == "PUB" &
                       r37$estimator == "System GMM", ]
p38 <- function(g) r38[r38$group == g & r38$estimator == "System GMM", ]
hiv <- r35[r35$sample == "C" & r35$estimator == "System GMM" & r35$term == "PUB", ]
stopifnot(nrow(pub16) == 1, nrow(hiv) == 1)

# Kiviet-Bruno: 90 per cent Fieller set from the bootstrap covariance,
# as script 28 does at 95 per cent
names(r26b) <- c("alpha_b", "beta_b")
kb <- r26[r26$method == "LSDVC B1 (init BB)", ]
stopifnot(nrow(kb) == 1)
Vkb <- cov(cbind(r26b$alpha_b, r26b$beta_b))
dimnames(Vkb) <- list(c("lag(LE, 1)", "PUB"), c("lag(LE, 1)", "PUB"))
kb90 <- longrun(c("lag(LE, 1)" = kb$alpha, PUB = kb$beta_PUB), Vkb,
                "lag(LE, 1)", level = 0.90)
kb90 <- kb90[kb90$term == "PUB", ]

pa <- a21("A drivers strictly exogenous")
stopifnot(isTRUE(all.equal(pa$beta, pub16$beta_SR)))
pb <- a21("B drivers predetermined (instruments from t-1)")
pc <- a21("C drivers endogenous (instruments from t-2)")
ce <- p37("C log CHE per head (PPP), endogenous")
cp <- p37("C log CHE per head (PPP), predetermined")
gl <- p38("low and lower-middle income in 2000")
gh <- p38("upper-middle and high income in 2000")

tab <- bind_rows(
  row("Reference specification",
      pub16$beta_SR, pub16$se_SR, pa$alpha,
      pub16$fieller90_lo, pub16$fieller90_hi, pub16$fieller90_type),
  row("+ expenditure per head, endogenous",
      ce$beta, ce$se, ce$alpha, ce$fieller90_lo, ce$fieller90_hi, ce$fieller90_type),
  row("+ expenditure per head, predetermined",
      cp$beta, cp$se, cp$alpha, cp$fieller90_lo, cp$fieller90_hi, cp$fieller90_type),
  row("Drivers predetermined", pb$beta, pb$se, pb$alpha, NA, NA, "undefined"),
  row("Drivers endogenous", pc$beta, pc$se, pc$alpha, NA, NA, "undefined"),
  row("With HIV prevalence, 1843 obs.",
      hiv$estimate, hiv$se, NA, hiv$fieller90_lo, hiv$fieller90_hi, hiv$fieller90_type),
  row("Bias-corrected LSDV",
      kb$beta_PUB, sd(r26b$beta_b), kb$alpha,
      kb90$fieller_lo, kb90$fieller_hi, kb90$fieller_type),
  row("Low, lower-middle income in 2000",
      gl$beta, gl$se, gl$alpha, gl$fieller90_lo, gl$fieller90_hi, gl$fieller90_type),
  row("Upper-middle, high income in 2000",
      gh$beta, gh$se, gh$alpha, gh$fieller90_lo, gh$fieller90_hi, gh$fieller90_type)
) %>%
  mutate(margin_5pct = abs(beta) + z90 * se,
         shift_lo = SHIFT * (beta - z90 * se),
         shift_hi = SHIFT * (beta + z90 * se),
         lr_lo = SHIFT * fieller90_lo,
         lr_hi = SHIFT * fieller90_hi)
# the persistence of the predetermined and endogenous classifications
# is above one, where the multiplier is not defined
stopifnot(all(tab$alpha[tab$fieller90_type == "undefined"] >= 1))

cat(sprintf("=== A planned shift of %d percentage points in the public share ===\n", SHIFT))
print(tab %>% mutate(across(c(beta, se, margin_5pct), ~ round(.x, 4)),
                     across(c(shift_lo, shift_hi, lr_lo, lr_hi), ~ round(.x, 2))) %>%
        dplyr::select(spec, beta, se, margin_5pct, shift_lo, shift_hi,
                      fieller90_type, lr_lo, lr_hi), width = 200)

write.csv(tab, file.path(OUT, "39_planning.csv"), row.names = FALSE)

# LaTeX body of Table 5
f2 <- function(x) sprintf("%+.2f", x)
lr_cell <- function(lo, hi, type) switch(type,
  "bounded"             = sprintf("$[%s, %s]$", f2(lo), f2(hi)),
  "exclusive/unbounded" = "unbounded",
  "whole line"          = "unbounded",
  "undefined"           = "not defined")
tex <- with(tab, sprintf("%s & $%+.4f$ (%.4f) & %.4f & $[%s, %s]$ & %s \\\\",
                         spec, beta, se, margin_5pct, f2(shift_lo), f2(shift_hi),
                         mapply(lr_cell, lr_lo, lr_hi, fieller90_type)))
writeLines(tex, file.path(OUT, "39_planning_table.tex"))
cat("\nSaved: 39_planning.csv, 39_planning_table.tex\n")

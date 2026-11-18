# 15: System GMM diagnostics with and without year effects
#
# The two-step System GMM of Table 3 with individual effects only and
# with two-way effects, instruments for the lagged outcome from t-2 to
# t-5, and with two-way effects and instruments from t-3 to t-5. For
# each: AR(1) and AR(2), the Hansen J test, Difference-in-Hansen, and
# the number of instruments (Section 6.2).
#
# In two-step estimation plm's sargan() uses the model's weight
# matrix, so it returns the robust Hansen J. plm has no
# Difference-in-Hansen test; it is computed here as
# C = J(system) - J(difference), with df = df(system) - df(difference).
#
# Reads  output/panel_with_ext.rds, shocks_panel.rds; data/df_le_long.rds
# Writes output/gmm_robustness_diagnostics.csv, gmm_robustness_coefs.csv

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(plm)
})

panel        <- readRDS(file.path(OUT, "panel_with_ext.rds"))
le_long      <- readRDS(file.path(DATA, "df_le_long.rds"))
shocks_panel <- readRDS(file.path(OUT, "shocks_panel.rds"))

dat <- panel %>%
  inner_join(le_long %>% select(iso3c, year, LE), by = c("iso3c","year")) %>%
  left_join(shocks_panel %>% select(iso3c, year, HIV_pct),
            by = c("iso3c","year")) %>%
  filter(!is.na(LE), !is.na(log_GNI), !is.na(EDU_GROSS),
         !is.na(PUB), !is.na(IMM_DPT), !is.na(SANITATION),
         !is.na(OOP_share)) %>%
  arrange(iso3c, year)

pd <- pdata.frame(dat, index = c("iso3c","year"))

# ---- diagnostics -------------------------------------------
diag_row <- function(model, m_diff_twin, label) {
  sg  <- sargan(model)            # Hansen J in two-step
  ab1 <- mtest(model, order = 1)
  ab2 <- mtest(model, order = 2)
  # Difference-in-Hansen: the system model against its difference twin
  dih_C <- dih_df <- dih_p <- NA
  if (!is.null(m_diff_twin)) {
    sgd <- sargan(m_diff_twin)
    dih_C  <- as.numeric(sg$statistic - sgd$statistic)
    dih_df <- as.numeric(sg$parameter - sgd$parameter)
    dih_p  <- pchisq(dih_C, dih_df, lower.tail = FALSE)
  }
  n_instr <- as.numeric(sg$parameter) + length(coef(model))
  data.frame(
    model       = label,
    n_instr     = n_instr,
    Hansen_chi2 = round(as.numeric(sg$statistic), 2),
    Hansen_df   = as.numeric(sg$parameter),
    Hansen_p    = round(as.numeric(sg$p.value), 4),
    AR1_p       = round(as.numeric(ab1$p.value), 4),
    AR2_p       = round(as.numeric(ab2$p.value), 4),
    DiffHansen_C  = round(dih_C, 2),
    DiffHansen_df = dih_df,
    DiffHansen_p  = round(dih_p, 4),
    stringsAsFactors = FALSE
  )
}

coef_row <- function(model, label) {
  co <- summary(model, robust = TRUE)$coefficients
  data.frame(model = label, term = rownames(co),
             beta = round(co[,1], 4), SE = round(co[,2], 4),
             z = round(co[,3], 3), p = round(co[,4], 4),
             row.names = NULL) %>%
    mutate(sig = ifelse(p<0.001,"***",ifelse(p<0.01,"**",
                 ifelse(p<0.05,"*",ifelse(p<0.1,".","")))))
}

# The formula is built with the lag numbers already in it, because
# pgmm evaluates it in its own environment and does not see the
# function arguments.
fit_gmm <- function(inst, eff, transf) {
  fml <- as.formula(sprintf(
    "LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION | lag(LE, %d:%d)",
    inst[1], inst[2]))
  pgmm(fml, data = pd, effect = eff, model = "twosteps",
       transformation = transf, collapse = TRUE)
}
fit_sys <- function(inst, eff) fit_gmm(inst, eff, "ld")
fit_dif <- function(inst, eff) fit_gmm(inst, eff, "d")

# ---- 1. individual against two-way effects -----------------
m_sys_ind  <- fit_sys(c(2,5), "individual")
m_dif_ind  <- fit_dif(c(2,5), "individual")
m_sys_tw   <- fit_sys(c(2,5), "twoways")       # the specification of Table 3
m_dif_tw   <- fit_dif(c(2,5), "twoways")

# ---- 2. instruments from t-3 -------------------------------
m_sys_tw_deep <- fit_sys(c(3,5), "twoways")
m_dif_tw_deep <- fit_dif(c(3,5), "twoways")

# ---- tables ------------------------------------------------
diag_tbl <- bind_rows(
  diag_row(m_sys_ind,     m_dif_ind,     "System GMM, individual effects (t-2:t-5)"),
  diag_row(m_sys_tw,      m_dif_tw,      "System GMM, two-way effects (t-2:t-5)"),
  diag_row(m_sys_tw_deep, m_dif_tw_deep, "System GMM, two-way effects (t-3:t-5)")
)

cat("\n=== AR, Hansen, Difference-in-Hansen, number of instruments ===\n")
print(diag_tbl, row.names = FALSE)

coef_tbl <- bind_rows(
  coef_row(m_sys_ind,     "individual (t-2:t-5)"),
  coef_row(m_sys_tw,      "two-way (t-2:t-5)"),
  coef_row(m_sys_tw_deep, "two-way (t-3:t-5)")
)
cat("\n=== PUB and IMM across specifications ===\n")
print(coef_tbl %>% filter(term %in% c("PUB","IMM_DPT")), row.names = FALSE)

cat("\n=== All coefficients, two-way effects, t-2:t-5 ===\n")
print(coef_row(m_sys_tw, "two-way (t-2:t-5)"), row.names = FALSE)

write.csv(diag_tbl, file.path(OUT, "gmm_robustness_diagnostics.csv"), row.names = FALSE)
write.csv(coef_tbl, file.path(OUT, "gmm_robustness_coefs.csv"),       row.names = FALSE)
cat("\nSaved: gmm_robustness_diagnostics.csv, gmm_robustness_coefs.csv\n")

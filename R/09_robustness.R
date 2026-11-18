# 09: robustness of the predictive design
#
#   1. the initial window: driver means over 2000-2003, 2000-2005 and
#      2000-2007
#   2. influence: the change in the in-sample R2 of the one-component
#      PLS on LE_shift when each country is left out
# Both feed Supplementary Figure S5 (script 10).
#
# Reads  output/panel_with_ext.rds, driver_initial.rds, le_targets.rds,
#        shape_pcs.rds; data/df_le_long.rds
# Writes output/sensitivity_window.csv, influence_pls_loo.csv

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pls)
})

panel <- readRDS(file.path(OUT, "panel_with_ext.rds"))
le_long <- readRDS(file.path(DATA, "df_le_long.rds"))
le_tg   <- readRDS(file.path(OUT, "le_targets.rds"))
pcs     <- readRDS(file.path(OUT, "shape_pcs.rds"))

# 1. initial window ----------------------------------------------

windows <- list(
  W3 = c(2000, 2003),
  W5 = c(2000, 2005),
  W7 = c(2000, 2007)
)

sens_results <- data.frame()
for (wn in names(windows)) {
  yrs <- windows[[wn]]
  d_init <- panel %>%
    filter(year >= yrs[1], year <= yrs[2]) %>%
    group_by(iso3c) %>%
    summarise(
      log_GNI = mean(log_GNI, na.rm=TRUE),
      EDU     = mean(EDU_GROSS, na.rm=TRUE),
      OOP     = mean(OOP_share, na.rm=TRUE),
      PNO     = mean(PRIV_NONOOP, na.rm=TRUE),
      PUB     = mean(PUB, na.rm=TRUE),
      EXT     = mean(EXT_share, na.rm=TRUE),
      IMM     = mean(IMM_DPT, na.rm=TRUE),
      SAN     = mean(SANITATION, na.rm=TRUE),
      .groups = "drop")
  cs <- d_init %>% inner_join(le_tg, by="iso3c") %>%
    inner_join(pcs, by="iso3c") %>% na.omit()

  X_vars <- c("log_GNI","EDU","OOP","PNO","PUB","EXT","IMM","SAN")
  X  <- as.matrix(cs[, X_vars])
  ys <- scale(cs$LE_shift)
  Ys <- scale(as.matrix(cs[, c("PC1","PC2","PC3")]))
  Xs <- scale(X)

  m1 <- plsr(ys ~ Xs, ncomp = 2, validation = "LOO", method = "oscorespls")
  m2 <- plsr(Ys ~ Xs, ncomp = 2, validation = "LOO", method = "oscorespls")
  R2_shift_w <- max(R2(m1, estimate="CV", intercept=FALSE)$val[1,1,])
  R2_shape_w <- max(R2(m2, estimate="CV", intercept=FALSE)$val[1,"PC2",])

  sens_results <- rbind(sens_results, data.frame(
    window = wn, years = paste(yrs, collapse="-"),
    n_countries = nrow(cs),
    R2_shift_max = round(R2_shift_w, 3),
    R2_PC2_max   = round(R2_shape_w, 3)
  ))
}
cat("=== Initial window ===\n")
print(sens_results)
write.csv(sens_results, file.path(OUT, "sensitivity_window.csv"), row.names = FALSE)

# 2. influence of single countries -------------------------------

d_init <- readRDS(file.path(OUT, "driver_initial.rds"))
cs <- d_init %>% inner_join(le_tg, by="iso3c") %>%
  inner_join(pcs, by="iso3c") %>% na.omit()

X_vars <- c("log_GNI_init","EDU_init","OOP_init","PRIV_NONOOP_init",
            "PUB_init","EXT_init","IMM_init","SAN_init")
X  <- as.matrix(cs[, X_vars])
y  <- cs$LE_shift

# in-sample R2 on all countries
ys <- scale(y); Xs <- scale(X)
m_full <- plsr(ys ~ Xs, ncomp = 1, method = "oscorespls")
fit_full <- drop(m_full$fitted.values[,,1])
R2_full <- 1 - sum((ys - fit_full)^2) / sum((ys - mean(ys))^2)

# the same with each country left out in turn
dR2 <- numeric(nrow(cs))
for (i in seq_len(nrow(cs))) {
  Xs_i <- scale(X[-i, ]); ys_i <- scale(y[-i])
  m_i  <- plsr(ys_i ~ Xs_i, ncomp = 1, method = "oscorespls")
  fit_i <- drop(m_i$fitted.values[,,1])
  R2_i  <- 1 - sum((ys_i - fit_i)^2) / sum((ys_i - mean(ys_i))^2)
  dR2[i] <- R2_full - R2_i   # change in R2 from removing the country
}

infl <- data.frame(iso3c = cs$iso3c, country = cs$country,
                   dR2 = round(dR2, 4)) %>% arrange(desc(abs(dR2)))

cat("\n=== Ten most influential countries (leave one out, PLS1) ===\n")
print(head(infl, 10))

write.csv(infl, file.path(OUT, "influence_pls_loo.csv"), row.names = FALSE)

cat("\nSaved: sensitivity_window.csv, influence_pls_loo.csv\n")

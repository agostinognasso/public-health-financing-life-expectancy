# 27: Supplementary Table S2, cross-sectional R2 with and without HIV and homicide
#
# OLS of LE_shift and of the shape components PC1-PC3 on the driver
# means over 2000-2005, with and without peak HIV prevalence and the
# mean homicide rate (Supplementary Section S6). The four financing shares sum to
# 100, so the external share is the omitted category. Gains are
# differences in R2, in percentage points.
#
# Reads  output/driver_initial.rds, le_targets.rds, shape_pcs.rds,
#        shocks_country.rds
# Writes output/27_table2.csv

source("00_common.R")

ds  <- readRDS(file.path(OUT, "driver_initial.rds"))
le  <- readRDS(file.path(OUT, "le_targets.rds"))
pcs <- readRDS(file.path(OUT, "shape_pcs.rds"))
shk <- readRDS(file.path(OUT, "shocks_country.rds"))
cs <- ds %>%
  inner_join(le,  by = "iso3c") %>%
  inner_join(pcs, by = "iso3c") %>%
  inner_join(shk, by = "iso3c") %>%
  na.omit()
cat(sprintf("cross-section for the HIV comparison: %d countries\n", nrow(cs)))

CTRL <- paste("log_GNI_init + EDU_init + IMM_init + SAN_init +",
              "OOP_init + PRIV_NONOOP_init + PUB_init")
r2 <- function(y) {
  f0 <- as.formula(paste(y, "~", CTRL))
  f1 <- update(f0, ~ . + HIV_peak + homicide)
  fh <- update(f0, ~ . + HIV_peak)
  m0 <- lm(f0, data = cs); m1 <- lm(f1, data = cs); mh <- lm(fh, data = cs)
  data.frame(R2_base = summary(m0)$r.squared,
             R2_full = summary(m1)$r.squared,
             R2_hiv_only = summary(mh)$r.squared,
             t_hiv = summary(m1)$coefficients["HIV_peak", "t value"],
             max_abs_t_other = max(abs(summary(m1)$coefficients[-c(1, which(rownames(summary(m1)$coefficients) == "HIV_peak")), "t value"])),
             adj_base = summary(m0)$adj.r.squared,
             adj_full = summary(m1)$adj.r.squared, n = nobs(m1))
}
d <- bind_rows(lapply(c("LE_shift","PC1","PC2","PC3"), r2))
d$target <- factor(c("Life expectancy\nshift (level)", "Shape PC1\n(volatility)",
                     "Shape PC2\n(stagnation)", "Shape PC3\n(smoothness)"),
                   levels = c("Life expectancy\nshift (level)", "Shape PC1\n(volatility)",
                              "Shape PC2\n(stagnation)", "Shape PC3\n(smoothness)"))
d <- d %>% mutate(gain = R2_full - R2_base, adj_gain = adj_full - adj_base)

cat("\n=== Supplementary Table S2 ===\n")
print(d %>% transmute(outcome = gsub("\n"," ",target),
                      n = n,
                      R2_base = round(R2_base, 2), R2_full = round(R2_full, 2),
                      gain_pp = round(100*gain, 1),
                      adj_gain_pp = round(100*adj_gain, 1),
                      gain_hiv_only_pp = round(100*(R2_hiv_only - R2_base), 1),
                      t_hiv = round(t_hiv, 2),
                      max_abs_t_other = round(max_abs_t_other, 2)),
      row.names = FALSE)

write.csv(d %>% mutate(target = gsub("\n", " ", target)),
          file.path(OUT, "27_table2.csv"), row.names = FALSE)
cat("Saved: 27_table2.csv\n")

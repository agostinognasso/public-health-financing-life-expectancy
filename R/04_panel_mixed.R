# 04: mixed model on the raw shares, between and within
#
# Random intercept by country, linear year trend, a COVID-year dummy,
# and each driver split into its country mean (between) and the
# deviation from it (within), as in Mundlak (1978). The fitted model
# is read by script 32, which compares it with the ILR model of
# script 08 by maximum likelihood.
#
# Reads  output/panel_with_ext.rds, region_map.rds; data/df_le_long.rds
# Writes output/panel_between_within.csv, panel_lmer_model.rds,
#        panel_summary.rds

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(lmerTest)
})

panel <- readRDS(file.path(OUT, "panel_with_ext.rds"))
reg   <- readRDS(file.path(OUT, "region_map.rds"))

panel <- panel %>%
  mutate(log_GNI_y = log_GNI,
         EDU_y     = EDU_GROSS,
         OOP_y     = OOP_share,
         PNO_y     = PRIV_NONOOP,
         PUB_y     = PUB,
         EXT_y     = EXT_share,
         IMM_y     = IMM_DPT,
         SAN_y     = SANITATION)

# complete cases on the drivers
panel <- panel %>%
  filter(!is.na(log_GNI_y), !is.na(EDU_y), !is.na(OOP_y),
         !is.na(PUB_y), !is.na(IMM_y), !is.na(SAN_y))

le_long <- readRDS(file.path(DATA, "df_le_long.rds"))
panel <- panel %>% inner_join(le_long %>% select(iso3c, year, LE),
                              by = c("iso3c","year")) %>%
  filter(!is.na(LE))

cat("Panel: ", nrow(panel), " rows, ",
    n_distinct(panel$iso3c), " countries\n")

# The four shares sum to 100, so only PUB and OOP enter; the private
# non-out-of-pocket and external shares together are the reference.
drv <- c("log_GNI_y","EDU_y","OOP_y","PUB_y","IMM_y","SAN_y")

mundlak <- panel %>%
  group_by(iso3c) %>%
  mutate(across(all_of(drv),
                list(bw = ~ mean(.x, na.rm = TRUE)),
                .names = "{.col}_bw")) %>%
  ungroup() %>%
  mutate(across(all_of(drv),
                ~ .x - get(paste0(cur_column(),"_bw")),
                .names = "{.col}_wn"))

mundlak$year_c <- mundlak$year - 2012

# 1. model -------------------------------------------------------

f_full <- as.formula(paste(
  "LE ~ year_c + is_covid_year +",
  paste(paste0(drv, "_bw"), collapse = " + "), " + ",
  paste(paste0(drv, "_wn"), collapse = " + "), " + ",
  "(1 | iso3c)"
))

cat("\nModel formula:\n  ", deparse(f_full), "\n")

m_full <- lmer(f_full, data = mundlak, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))

co <- coef(summary(m_full))
co <- as.data.frame(co)
co$term <- rownames(co)
rownames(co) <- NULL

bw_rows <- co %>% filter(grepl("_bw$", term)) %>%
  mutate(var = sub("_y_bw$","", term),
         type = "between (cross-country)")
wn_rows <- co %>% filter(grepl("_wn$", term)) %>%
  mutate(var = sub("_y_wn$","", term),
         type = "within (year-to-year)")

paneltable <- bind_rows(bw_rows, wn_rows) %>%
  select(var, type, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`) %>%
  rename(SE = `Std. Error`, t = `t value`, p = `Pr(>|t|)`) %>%
  mutate(across(c(Estimate, SE, t), ~ round(.x, 4)),
         p = round(p, 4),
         sig = ifelse(p < 0.001, "***",
               ifelse(p < 0.01, "**",
               ifelse(p < 0.05, "*",
               ifelse(p < 0.1, ".", "")))))

cat("\n=== Between and within coefficients ===\n")
print(paneltable %>% arrange(var, type))

write.csv(paneltable, file.path(OUT, "panel_between_within.csv"), row.names = FALSE)

# intraclass correlation
var_comp <- as.data.frame(VarCorr(m_full))
icc <- var_comp$vcov[var_comp$grp == "iso3c"] /
       sum(var_comp$vcov)
cat(sprintf("\nICC (country variance / total): %.3f\n", icc))
cat("  ", round(icc*100,1), "% of the variance of LE is between countries,\n",
    "  ", round((1-icc)*100,1), "% within countries\n")

# 2. marginal and conditional R2 ---------------------------------

# Nakagawa and Schielzeth (2013)
r2_nakagawa <- function(model) {
  v_fix <- var(predict(model, re.form = NA))
  vc <- as.data.frame(VarCorr(model))
  v_re <- sum(vc$vcov[vc$grp != "Residual"])
  v_res <- vc$vcov[vc$grp == "Residual"]
  list(
    marginal    = v_fix / (v_fix + v_re + v_res),
    conditional = (v_fix + v_re) / (v_fix + v_re + v_res)
  )
}
r2s <- r2_nakagawa(m_full)
cat(sprintf("\nR2 marginal (fixed effects only):     %.3f\n", r2s$marginal))
cat(sprintf("R2 conditional (with random effects): %.3f\n", r2s$conditional))

# 3. largest coefficients ----------------------------------------

top_bw <- paneltable %>% filter(type == "between (cross-country)") %>%
  arrange(desc(abs(Estimate))) %>% head(5)
top_wn <- paneltable %>% filter(type == "within (year-to-year)") %>%
  arrange(desc(abs(Estimate))) %>% head(5)

cat("\n--- five largest between-country coefficients ---\n")
print(top_bw)
cat("\n--- five largest within-country coefficients ---\n")
print(top_wn)

# 4. save --------------------------------------------------------

saveRDS(m_full, file.path(OUT, "panel_lmer_model.rds"))
saveRDS(list(icc = icc, r2 = r2s, paneltable = paneltable),
        file.path(OUT, "panel_summary.rds"))

cat("\nSaved: panel_between_within.csv, panel_lmer_model.rds, panel_summary.rds\n")

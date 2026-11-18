# 25: missingness, panel gaps and the inclusion rule (Section 6.2)
#
# The analysis is complete-case. The script reports the share of the
# country-by-year rectangle that is analysed, missingness by
# variable, how analysed and dropped cells differ in life expectancy
# and income, the internal gaps in each country's series and the
# System GMM on the longest contiguous run, and the inclusion rule
# re-applied on 2000-2024.
#
# Reads  data/df_panel_drivers.rds, data/df_le_long.rds, and data/
#        through build_panel()
# Writes output/25_item_missingness.csv, 25_country_gaps.csv,
#        25_gap_sensitivity.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

raw <- readRDS(file.path(DATA, "df_panel_drivers.rds"))
le  <- readRDS(file.path(DATA, "df_le_long.rds"))
dat <- build_panel(zero_rule = "multiplicative")

# ---- 1. completeness --------------------------------------
cat("=== Completeness of the analytical rectangle ===\n")
nC <- n_distinct(dat$iso3c); nY <- length(2000:2024)
cat(sprintf("countries x years          : %d x %d = %d cells\n", nC, nY, nC * nY))
cat(sprintf("observations analysed      : %d\n", nrow(dat)))
cat(sprintf("missing                    : %d (%.1f%%)\n",
            nC * nY - nrow(dat), 100 * (1 - nrow(dat) / (nC * nY))))

cat("\nobservations per country: ")
opc <- dat %>% count(iso3c)
cat(sprintf("min %d, Q1 %.0f, median %.0f, Q3 %.0f, max %d\n",
            min(opc$n), quantile(opc$n, .25), median(opc$n),
            quantile(opc$n, .75), max(opc$n)))
cat(sprintf("countries with the full 25 years: %d (%.0f%%)\n",
            sum(opc$n == 25), 100 * mean(opc$n == 25)))

# ---- 2. missingness by variable ---------------------------
cat("\n=== Item-level missingness in the 144-country raw panel ===\n")
vv <- c(LE = NA, log_GNI = "log_GNI", EDU_GROSS = "EDU_GROSS",
        PUB = "PUB", OOP_share = "OOP_share", PRIV = "PRIV",
        IMM_DPT = "IMM_DPT", SANITATION = "SANITATION")
rawle <- raw %>% left_join(le %>% dplyr::select(iso3c, year, LE),
                           by = c("iso3c", "year"))
miss <- bind_rows(lapply(c("LE", names(vv)[-1]), function(v) {
  data.frame(variable = v, n_missing = sum(is.na(rawle[[v]])),
             pct = round(100 * mean(is.na(rawle[[v]])), 1))
}))
print(miss)

# ---- 3. analysed against dropped cells ---------------------
cat("\n=== Analysed against dropped cells ===\n")
grid <- expand.grid(iso3c = unique(raw$iso3c), year = 2000:2024,
                    stringsAsFactors = FALSE)
flag <- grid %>%
  left_join(dat %>% transmute(iso3c, year, analysed = 1L),
            by = c("iso3c", "year")) %>%
  mutate(analysed = ifelse(is.na(analysed), 0L, 1L)) %>%
  left_join(le %>% dplyr::select(iso3c, year, LE), by = c("iso3c", "year")) %>%
  left_join(raw %>% dplyr::select(iso3c, year, log_GNI, region),
            by = c("iso3c", "year"))

cat("mean life expectancy and income, analysed vs dropped cells:\n")
print(flag %>% group_by(analysed) %>%
        summarise(n = n(), LE = round(mean(LE, na.rm = TRUE), 2),
                  log_GNI = round(mean(log_GNI, na.rm = TRUE), 3),
                  .groups = "drop"))
d_le <- diff(tapply(flag$LE, flag$analysed, mean, na.rm = TRUE))
cat(sprintf("\ndifference in mean LE (analysed minus dropped): %+.2f years\n", d_le))
gl <- glm(analysed ~ LE + log_GNI, data = flag, family = binomial)
cat("\nlogit of being analysed on the outcome and income:\n")
print(round(summary(gl)$coefficients, 4))
cat("\nshare of cells analysed, by region:\n")
print(as.data.frame(flag %>% filter(!is.na(region)) %>% group_by(region) %>%
        summarise(cells = n(), pct_analysed = round(100 * mean(analysed), 1),
                  .groups = "drop") %>% arrange(pct_analysed)), row.names = FALSE)

# ---- 4. internal gaps and the GMM moment conditions --------
cat("\n\n=== Internal gaps in the estimation panel ===\n")
gaps <- dat %>% arrange(iso3c, year) %>% group_by(iso3c) %>%
  summarise(n = n(), span = max(year) - min(year) + 1,
            n_gaps = sum(diff(year) > 1),
            largest_gap = ifelse(any(diff(year) > 1), max(diff(year)) - 1, 0L),
            .groups = "drop")
cat(sprintf("countries with at least one internal gap: %d of %d (%.0f%%)\n",
            sum(gaps$n_gaps > 0), nrow(gaps), 100 * mean(gaps$n_gaps > 0)))
cat(sprintf("observations in those countries          : %d (%.0f%% of the panel)\n",
            sum(gaps$n[gaps$n_gaps > 0]), 100 * sum(gaps$n[gaps$n_gaps > 0]) / nrow(dat)))
cat(sprintf("largest single gap                       : %d years\n", max(gaps$largest_gap)))
cat("\nA gap removes the lagged differences and levels that straddle it,\n")
cat("so the instruments available differ across countries; plm forms the\n")
cat("moments from the available pairs only.\n")

cat("\n--- sensitivity: contiguous runs only ---\n")
contig <- dat %>% arrange(iso3c, year) %>% group_by(iso3c) %>%
  mutate(run = cumsum(c(1, diff(year) != 1))) %>%
  group_by(iso3c, run) %>% mutate(run_len = n()) %>% ungroup() %>%
  group_by(iso3c) %>% filter(run == run[which.max(run_len)]) %>% ungroup()
cat(sprintf("longest contiguous run per country: %d obs, %d countries\n",
            nrow(contig), n_distinct(contig$iso3c)))
fg <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION | lag(LE, 2:5)
m_all <- pgmm(fg, data = pdata.frame(dat, index = c("iso3c", "year")),
              effect = "twoways", model = "twosteps",
              transformation = "ld", collapse = TRUE)
m_ctg <- tryCatch(pgmm(fg, data = pdata.frame(contig, index = c("iso3c", "year")),
                       effect = "twoways", model = "twosteps",
                       transformation = "ld", collapse = TRUE),
                  error = function(e) NULL)
cmp <- function(m, lab) {
  if (is.null(m)) return(NULL)
  b <- coef(m); V <- pgmm_vcov(m)
  data.frame(sample = lab, alpha = b[["lag(LE, 1)"]],
             PUB = b[["PUB"]], se_PUB = sqrt(V["PUB", "PUB"]),
             IMM = b[["IMM_DPT"]], se_IMM = sqrt(V["IMM_DPT", "IMM_DPT"]),
             n_inst = dim(m$W[[1]])[2], row.names = NULL)
}
out <- bind_rows(cmp(m_all, "all observations"),
                 cmp(m_ctg, "longest contiguous run only"))
out$n_obs       <- c(nrow(dat), nrow(contig))[seq_len(nrow(out))]
out$n_countries <- c(n_distinct(dat$iso3c),
                     n_distinct(contig$iso3c))[seq_len(nrow(out))]
print(out %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        mutate(p_PUB = round(2 * pnorm(-abs(PUB / se_PUB)), 3),
               p_IMM = round(2 * pnorm(-abs(IMM / se_IMM)), 3)))

# ---- 5. the inclusion rule on the estimation window --------
cat("\n\n=== The inclusion rule and the estimation window ===\n")
cat("The inclusion rule asks for at least 30 non-missing LE observations\n")
cat("over 1960-2024; here it is re-applied on 2000-2024.\n\n")
le_cnt <- le %>% group_by(iso3c) %>%
  summarise(n_6024 = sum(!is.na(LE)),
            n_0024 = sum(!is.na(LE) & year >= 2000 & year <= 2024),
            .groups = "drop")
inpanel <- le_cnt %>% filter(iso3c %in% unique(dat$iso3c))
cat(sprintf("countries in the analytical panel: %d\n", nrow(inpanel)))
cat(sprintf("of which fewer than 20 LE observations within 2000-2024: %d\n",
            sum(inpanel$n_0024 < 20)))
cat(sprintf("of which fewer than 15 within 2000-2024                : %d\n",
            sum(inpanel$n_0024 < 15)))
cat("\ncontribution of each country to the estimated panel:\n")
print(as.data.frame(opc %>% arrange(n) %>% head(12)), row.names = FALSE)
cat(sprintf("\ncountries contributing fewer than 10 country-years: %d\n",
            sum(opc$n < 10)))
cat(sprintf("countries contributing fewer than 5                : %d\n",
            sum(opc$n < 5)))

cat("\n--- sensitivity: require at least 15 observations in 2000-2024 ---\n")
keep <- opc$iso3c[opc$n >= 15]
d15 <- dat %>% filter(iso3c %in% keep)
cat(sprintf("panel: %d obs, %d countries\n", nrow(d15), n_distinct(d15$iso3c)))
m15 <- tryCatch(pgmm(fg, data = pdata.frame(d15, index = c("iso3c", "year")),
                     effect = "twoways", model = "twosteps",
                     transformation = "ld", collapse = TRUE),
                error = function(e) NULL)
print(bind_rows(cmp(m_all, "rule on 1960-2024"),
                cmp(m15, ">= 15 obs in 2000-2024")) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        mutate(p_PUB = round(2 * pnorm(-abs(PUB / se_PUB)), 3),
               p_IMM = round(2 * pnorm(-abs(IMM / se_IMM)), 3)))

write.csv(miss,  file.path(OUT, "25_item_missingness.csv"), row.names = FALSE)
write.csv(gaps,  file.path(OUT, "25_country_gaps.csv"), row.names = FALSE)
write.csv(out,   file.path(OUT, "25_gap_sensitivity.csv"), row.names = FALSE)
cat("\nSaved: 25_item_missingness.csv, 25_country_gaps.csv,\n")
cat("25_gap_sensitivity.csv\n")

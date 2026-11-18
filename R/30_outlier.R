# 30: an implausible outcome value (Section 6.2)
#
# The World Bank series gives the Central African Republic a life
# expectancy of 14.7 years in 2009, against about 48 in the years
# either side. The script lists the year-on-year changes larger
# than five years across the panel and re-estimates the static and
# System GMM models without CAF 2009 and without CAF.
#
# Reads  data/df_le_long.rds, and data/ through build_panel()
# Writes output/30_outlier_influence.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))

le <- readRDS(file.path(DATA, "df_le_long.rds"))
cat("=== Central African Republic, life expectancy series ===\n")
print(le %>% filter(iso3c == "CAF", year >= 2003, year <= 2015) %>%
        arrange(year) %>% dplyr::select(year, LE) %>% as.data.frame())

caf <- le %>% filter(iso3c == "CAF") %>% arrange(year)
cat("\nlargest year-on-year change in the CAF series:\n")
d <- data.frame(year = caf$year[-1], change = diff(caf$LE))
print(d %>% arrange(change) %>% head(3))
print(d %>% arrange(desc(change)) %>% head(3))

cat("\n=== How many implausible annual jumps are in the panel? ===\n")
jumps <- le %>% arrange(iso3c, year) %>% group_by(iso3c) %>%
  mutate(d = LE - dplyr::lag(LE)) %>% ungroup() %>% filter(!is.na(d))
cat(sprintf("year-on-year changes exceeding 5 years in size: %d\n",
            sum(abs(jumps$d) > 5)))
print(jumps %>% filter(abs(d) > 5) %>% arrange(d) %>%
        dplyr::select(country, iso3c, year, LE, d) %>% head(10) %>% as.data.frame())

# ---- influence on the headline estimates -------------------
dat  <- build_panel(zero_rule = "multiplicative")
flag <- dat$iso3c == "CAF" & dat$year == 2009
cat(sprintf("\nthe observation is in the estimation panel: %s\n", any(flag)))

fg <- LE ~ lag(LE, 1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION | lag(LE, 2:5)
run <- function(d, lab) {
  pd <- pdata.frame(d, index = c("iso3c", "year"))
  m1 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION,
            data = pd, model = "within")
  m2 <- plm(LE ~ log_GNI + EDU_GROSS + PUB + OOP_share + IMM_DPT + SANITATION,
            data = pd, model = "within", effect = "twoways")
  g  <- tryCatch(pgmm(fg, data = pd, effect = "twoways", model = "twosteps",
                      transformation = "ld", collapse = TRUE),
                 error = function(e) NULL)
  V1 <- plm::vcovHC(m1, cluster = "group", method = "arellano")
  out <- data.frame(sample = lab, n = nrow(d),
                    FE_oneway = coef(m1)["PUB"],
                    se_oneway = sqrt(V1["PUB", "PUB"]),
                    FE_twoway = coef(m2)["PUB"])
  if (!is.null(g)) {
    b <- coef(g); V <- pgmm_vcov(g)
    out$alpha <- b[["lag(LE, 1)"]]
    out$SysGMM_PUB <- b[["PUB"]]
    out$se_PUB <- sqrt(V["PUB", "PUB"])
    out$SysGMM_IMM <- b[["IMM_DPT"]]
    out$p_IMM <- 2 * pnorm(-abs(b[["IMM_DPT"]] / sqrt(V["IMM_DPT", "IMM_DPT"])))
  }
  out
}

cmp <- bind_rows(
  run(dat, "full panel"),
  run(dat %>% filter(!flag), "dropping CAF 2009"),
  run(dat %>% filter(iso3c != "CAF"), "dropping CAF entirely"))
cat("\n=== Influence of the outlier on the headline estimates ===\n")
print(cmp %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

d1 <- cmp$FE_oneway[1]; d2 <- cmp$FE_oneway[2]
cat(sprintf("\none-way FE on PUB: %+.4f with the observation, %+.4f without (%+.0f%%)\n",
            d1, d2, 100 * (d2 - d1) / abs(d1)))

write.csv(cmp, file.path(OUT, "30_outlier_influence.csv"), row.names = FALSE)
cat("\nSaved: 30_outlier_influence.csv\n")

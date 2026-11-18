# 11: HIV prevalence and the homicide rate
#
#   HIV_pct   SH.DYN.AIDS.ZS, % of the population aged 15-49
#   HOMICIDE  VC.IHR.PSRC.P5, intentional homicides per 100,000
#
# Each series is read from its cached extract in data/ and downloaded
# through WDI only if the extract is missing. The panel is kept as it
# is; the country summaries are the 2000-2024 peak and mean of HIV
# prevalence and the mean homicide rate, used by the cross-section of
# Supplementary Table S2 (script 27).
#
# Reads  data/HIV_pct_shock.rds, data/HOMICIDE_shock.rds
# Writes output/shocks_panel.rds, shocks_country.rds

source("00_common.R")

suppressPackageStartupMessages({
  library(WDI); library(dplyr); library(tidyr); library(countrycode)
})


shock_inds <- c(
  HIV_pct  = "SH.DYN.AIDS.ZS",   # HIV, % of population 15-49
  HOMICIDE = "VC.IHR.PSRC.P5"    # intentional homicides per 100,000
)

dl_one <- function(name, code) {
  cf <- file.path(DATA, paste0(name, "_shock.rds"))
  if (file.exists(cf)) return(readRDS(cf))
  x <- tryCatch(
    WDI::WDI(country = "all", indicator = setNames(code, name),
             start = 1995, end = 2024, extra = FALSE, cache = NULL),
    error = function(e) NULL
  )
  if (!is.null(x)) saveRDS(x, cf)
  x
}

dl_list <- lapply(names(shock_inds), function(n) dl_one(n, shock_inds[n]))
names(dl_list) <- names(shock_inds)
dl_list <- dl_list[!sapply(dl_list, is.null)]

add_iso <- function(d) d %>%
  mutate(iso3c = suppressWarnings(
    countrycode::countrycode(iso2c,"iso2c","iso3c"))) %>%
  filter(!is.na(iso3c), iso3c != "")

shocks <- Reduce(function(x,y) full_join(x,y,by=c("iso2c","country","year")),
                 dl_list) %>%
  add_iso() %>%
  select(iso3c, year, any_of(names(shock_inds)))

cat("HIV and homicide series for ", n_distinct(shocks$iso3c), " countries\n")
cat("Non-missing values by series:\n")
print(sapply(shocks[,3:ncol(shocks)], function(x) sum(!is.na(x))))

# country summaries over 2000-2024 -------------------------------
if (!"HOMICIDE" %in% names(shocks)) shocks$HOMICIDE <- NA_real_

shocks_country <- shocks %>%
  filter(year >= 2000, year <= 2024) %>%
  group_by(iso3c) %>%
  summarise(
    HIV_peak  = suppressWarnings(max(HIV_pct, na.rm = TRUE)),
    HIV_mean  = mean(HIV_pct, na.rm = TRUE),
    homicide  = mean(HOMICIDE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric),
                ~ ifelse(is.finite(.x), .x, NA_real_)))

cat("\nCountry summaries, medians:\n")
print(sapply(shocks_country[,-1], function(x) median(x, na.rm=TRUE)))

saveRDS(shocks_country, file.path(OUT, "shocks_country.rds"))
saveRDS(shocks, file.path(OUT, "shocks_panel.rds"))

cat("\nSaved: shocks_country.rds, shocks_panel.rds\n")

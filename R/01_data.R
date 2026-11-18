# 01: country panel and cross-sectional targets
#
# Adds the external financing share, EXT = 100 - PUB - PRIV, to the
# driver panel, and builds the cross-section used by the predictive
# design of Supplementary Section S1:
#   driver_initial  driver means over 2000-2005 (the predictors)
#   driver_full     driver means over 2000-2024 (for comparison)
#   le_targets      LE_shift = mean LE 2020-2024 minus mean LE 2000-2004
#
# Reads  data/df_le_long.rds, data/df_panel_drivers.rds
# Writes output/panel_with_ext.rds, driver_initial.rds, driver_full.rds,
#        le_targets.rds, region_map.rds

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})


le_long <- readRDS(file.path(DATA, "df_le_long.rds"))
panel   <- readRDS(file.path(DATA, "df_panel_drivers.rds"))

# 1. external share ---------------------------------------------
panel <- panel %>%
  mutate(EXT_share = pmax(100 - PUB - PRIV, 0))

cat("EXT_share summary (% CHE):\n")
print(summary(panel$EXT_share))

# 2. driver means over the initial window ------------------------
INIT_Y0 <- 2000; INIT_Y1 <- 2005
POST_Y0 <- 2005; POST_Y1 <- 2024
SHIFT_PRE_Y0  <- 2000; SHIFT_PRE_Y1  <- 2004
SHIFT_POST_Y0 <- 2020; SHIFT_POST_Y1 <- 2024

mk_summary <- function(df, y0, y1, suffix) {
  df %>%
    filter(year >= y0, year <= y1) %>%
    group_by(iso3c, country) %>%
    summarise(
      log_GNI    = mean(log_GNI,    na.rm = TRUE),
      EDU        = mean(EDU_GROSS,  na.rm = TRUE),
      OOP        = mean(OOP_share,  na.rm = TRUE),
      PRIV_NONOOP= mean(PRIV_NONOOP,na.rm = TRUE),
      PUB        = mean(PUB,        na.rm = TRUE),
      EXT        = mean(EXT_share,  na.rm = TRUE),
      IMM        = mean(IMM_DPT,    na.rm = TRUE),
      SAN        = mean(SANITATION, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric),
                  ~ ifelse(is.finite(.x), .x, NA_real_))) %>%
    rename_with(~ paste0(.x, "_", suffix),
                .cols = -c(iso3c, country))
}

ds_init <- mk_summary(panel, INIT_Y0, INIT_Y1, "init")
ds_full <- mk_summary(panel, INIT_Y0, 2024,    "full")

cat("\ndriver means 2000-2005: ", nrow(ds_init), " countries\n")
cat("driver means 2000-2024: ", nrow(ds_full), " countries\n")

# 3. outcome targets ---------------------------------------------
le_pre  <- le_long %>% filter(year >= SHIFT_PRE_Y0,  year <= SHIFT_PRE_Y1) %>%
  group_by(iso3c) %>% summarise(LE_pre  = mean(LE, na.rm=TRUE), .groups="drop")
le_post <- le_long %>% filter(year >= SHIFT_POST_Y0, year <= SHIFT_POST_Y1) %>%
  group_by(iso3c) %>% summarise(LE_post = mean(LE, na.rm=TRUE), .groups="drop")

le_targets <- le_pre %>% inner_join(le_post, by = "iso3c") %>%
  mutate(LE_shift  = LE_post - LE_pre,
         LE_level  = LE_post)   # mean LE over 2020-2024

cat("\nLE_shift summary (years):\n"); print(summary(le_targets$LE_shift))
cat("LE_level summary (years):\n"); print(summary(le_targets$LE_level))

# 4. regions, for the leave-one-region-out designs ---------------
region_map <- suppressWarnings(
  countrycode::codelist[, c("iso3c", "region", "continent")]
)

# 5. save --------------------------------------------------------
saveRDS(ds_init,    file.path(OUT, "driver_initial.rds"))
saveRDS(ds_full,    file.path(OUT, "driver_full.rds"))
saveRDS(le_targets, file.path(OUT, "le_targets.rds"))
saveRDS(panel,      file.path(OUT, "panel_with_ext.rds"))
saveRDS(region_map, file.path(OUT, "region_map.rds"))

cat("\nSaved: driver_initial.rds, driver_full.rds, le_targets.rds,\n")
cat("panel_with_ext.rds, region_map.rds\n")

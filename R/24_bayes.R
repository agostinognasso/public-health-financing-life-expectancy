# 24: Bayesian hierarchical model (Supplementary Section S2)
#
# rstanarm::stan_lmer with country intercepts and regional slopes on
# the within components of immunisation, sanitation and enrolment;
# four chains of 2000 iterations. Convergence diagnostics for every
# reported parameter, a refit at adapt_delta = 0.99 for the divergent
# transitions, and two alternative priors on the random-effect
# covariance, because the regional slopes rest on seven regions, one
# of them with three countries.
#
# Reads  output/panel_with_ext.rds, shocks_panel.rds; data/df_le_long.rds
# Writes output/24_fixed_effects.csv, 24_slope_scale.csv,
#        24_convergence.csv, 24_prior_sensitivity.csv,
#        24_region_slopes.csv, 24_adapt_delta.csv

source("00_common.R")
suppressPackageStartupMessages({
  library(rstanarm); library(posterior)
})
options(mc.cores = min(4, parallel::detectCores()))

eps <- 0.05
panel <- readRDS(file.path(OUT, "panel_with_ext.rds"))
le    <- readRDS(file.path(DATA, "df_le_long.rds"))
shk   <- readRDS(file.path(OUT, "shocks_panel.rds"))

dat <- panel %>%
  inner_join(le %>% dplyr::select(iso3c, year, LE), by = c("iso3c", "year")) %>%
  left_join(shk %>% dplyr::select(iso3c, year, HIV_pct), by = c("iso3c", "year")) %>%
  filter(!is.na(LE), !is.na(log_GNI), !is.na(EDU_GROSS), !is.na(PUB),
         !is.na(IMM_DPT), !is.na(SANITATION), !is.na(OOP_share),
         !is.na(HIV_pct), !is.na(region), region != "")

# ILR coordinates with the additive replacement of script 08
b_ilr <- function(df) {
  P <- pmax(df$PUB, eps); O <- pmax(df$OOP_share, eps)
  N <- pmax(df$PRIV_NONOOP, eps); E <- pmax(df$EXT_share, eps)
  t <- P + O + N + E; P <- P/t*100; O <- O/t*100; N <- N/t*100; E <- E/t*100
  data.frame(ILR1 = sqrt(3/4) * (log(P) - (log(O)+log(N)+log(E))/3),
             ILR2 = sqrt(2/3) * ((log(O)+log(N))/2 - log(E)),
             ILR3 = sqrt(1/2) * (log(O) - log(N)))
}
dat <- bind_cols(dat, b_ilr(dat)) %>%
  group_by(iso3c) %>%
  mutate(log_GNI_bw = mean(log_GNI), log_GNI_wn = log_GNI - log_GNI_bw,
         EDU_bw = mean(EDU_GROSS),   EDU_wn = EDU_GROSS - EDU_bw,
         IMM_bw = mean(IMM_DPT),     IMM_wn = IMM_DPT - IMM_bw,
         SAN_bw = mean(SANITATION),  SAN_wn = SANITATION - SAN_bw,
         HIV_bw = mean(HIV_pct),     HIV_wn = HIV_pct - HIV_bw,
         ILR1_bw = mean(ILR1),       ILR1_wn = ILR1 - ILR1_bw,
         ILR2_bw = mean(ILR2),       ILR2_wn = ILR2 - ILR2_bw,
         ILR3_bw = mean(ILR3),       ILR3_wn = ILR3 - ILR3_bw) %>%
  ungroup() %>% mutate(year_c = year - 2012)

P2S <- c("log_GNI_bw","EDU_bw","IMM_bw","SAN_bw","HIV_bw",
         "ILR1_bw","ILR2_bw","ILR3_bw","log_GNI_wn","EDU_wn","IMM_wn",
         "SAN_wn","HIV_wn","ILR1_wn","ILR2_wn","ILR3_wn")
for (v in P2S) dat[[paste0(v, "_s")]] <- as.numeric(scale(dat[[v]]))

cat("Bayesian panel:", nrow(dat), "obs,", n_distinct(dat$iso3c),
    "countries,", n_distinct(dat$region), "regions\n")
cat("\ncountries per region:\n")
print(dat %>% group_by(region) %>%
        summarise(countries = n_distinct(iso3c), obs = n(), .groups = "drop") %>%
        arrange(countries))

f <- as.formula(paste("LE ~ year_c + is_covid_year +",
                      paste(paste0(P2S, "_s"), collapse = " + "),
                      "+ (1 | iso3c) + (1 + IMM_wn_s + SAN_wn_s + EDU_wn_s | region)"))

fit_one <- function(pcov, label, chains = 4, iter = 2000, adapt_delta = 0.95) {
  cat(sprintf("\n--- fitting: %s (%d chains x %d iter) ---\n", label, chains, iter))
  t0 <- Sys.time()
  m <- stan_lmer(f, data = dat, chains = chains, iter = iter,
                 warmup = iter / 2, refresh = 0, seed = 42,
                 prior = normal(0, 1), prior_intercept = normal(70, 15),
                 prior_aux = exponential(1), prior_covariance = pcov,
                 adapt_delta = adapt_delta)
  cat(sprintf("    %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  m
}

# ---- main fit ----------------------------------------------
m_main <- fit_one(decov(regularization = 1, concentration = 1,
                        shape = 1, scale = 1), "main, decov(1,1,1,1)")

# ---- convergence diagnostics -------------------------------
dr  <- as_draws_array(m_main)
sm  <- summarise_draws(dr, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
show <- grepl("_s$|^year_c$|^is_covid|^Sigma\\[region|^sigma$|^\\(Intercept\\)$",
              sm$variable)
diag_tbl <- sm[show, ]
cat("\n=== Convergence diagnostics for displayed parameters ===\n")
print(as.data.frame(diag_tbl %>%
        mutate(across(c(mean, sd, rhat), ~ round(.x, 4)),
               across(c(ess_bulk, ess_tail), ~ round(.x)))), row.names = FALSE)

cat(sprintf("\nworst R-hat among displayed parameters : %.4f\n", max(diag_tbl$rhat, na.rm = TRUE)))
cat(sprintf("smallest bulk ESS                     : %.0f\n", min(diag_tbl$ess_bulk, na.rm = TRUE)))
cat(sprintf("smallest tail ESS                     : %.0f\n", min(diag_tbl$ess_tail, na.rm = TRUE)))
cat(sprintf("post-warmup draws                     : %d\n", nrow(as.data.frame(m_main))))
nd <- tryCatch(bayesplot::nuts_params(m_main), error = function(e) NULL)
if (!is.null(nd))
  cat(sprintf("divergent transitions                 : %d\n",
              sum(nd$Value[nd$Parameter == "divergent__"])))

# ---- the divergences: refit with a smaller step size --------
# Divergent transitions can signal a part of the posterior the
# sampler does not explore. A higher adapt_delta shrinks the step
# size; the refit reports the divergences left and the largest
# change in a posterior mean.
m_ad <- fit_one(decov(regularization = 1, concentration = 1,
                      shape = 1, scale = 1),
                "main, adapt_delta = 0.99", adapt_delta = 0.99)
sm_ad <- summarise_draws(as_draws_array(m_ad),
                         "mean", "sd", "rhat", "ess_bulk", "ess_tail")
sm_ad <- sm_ad[match(diag_tbl$variable, sm_ad$variable), ]
nd_ad <- bayesplot::nuts_params(m_ad)
ad_tbl <- data.frame(
  adapt_delta = c(0.95, 0.99),
  divergent = c(if (is.null(nd)) NA else sum(nd$Value[nd$Parameter == "divergent__"]),
                sum(nd_ad$Value[nd_ad$Parameter == "divergent__"])),
  rhat_max = c(max(diag_tbl$rhat, na.rm = TRUE), max(sm_ad$rhat, na.rm = TRUE)),
  ess_bulk_min = c(min(diag_tbl$ess_bulk, na.rm = TRUE),
                   min(sm_ad$ess_bulk, na.rm = TRUE)),
  ess_tail_min = c(min(diag_tbl$ess_tail, na.rm = TRUE),
                   min(sm_ad$ess_tail, na.rm = TRUE)),
  max_shift_in_sd = c(0, max(abs(sm_ad$mean - diag_tbl$mean) / diag_tbl$sd,
                             na.rm = TRUE)))
cat("\n=== Main model at adapt_delta 0.95 and 0.99 ===\n")
print(ad_tbl %>% mutate(across(c(rhat_max, max_shift_in_sd), ~ round(.x, 4)),
                        across(c(ess_bulk_min, ess_tail_min), ~ round(.x))),
      row.names = FALSE)
cat("max_shift_in_sd: largest change in a posterior mean among the displayed\n")
cat("parameters, in posterior standard deviations of the 0.95 fit.\n")

# ---- random-effect SDs and prior sensitivity ---------------
get_sds <- function(m, label) {
  vc <- as.data.frame(VarCorr(m))
  vc %>% filter(is.na(var2)) %>%
    transmute(spec = label, group = grp, term = var1, sd = sdcor)
}

m_s1 <- fit_one(decov(regularization = 2, concentration = 1,
                      shape = 1, scale = 0.5), "tighter: decov(2,1,1,0.5)")
m_s2 <- fit_one(decov(regularization = 1, concentration = 1,
                      shape = 1, scale = 2),   "looser: decov(1,1,1,2)")

sds <- bind_rows(get_sds(m_main, "main decov(1,1,1,1)"),
                 get_sds(m_s1,   "tighter decov(2,1,1,0.5)"),
                 get_sds(m_s2,   "looser decov(1,1,1,2)"))
sds_ad <- get_sds(m_ad, "main decov(1,1,1,1), adapt_delta 0.99")
cat("\n=== Regional slope SDs, main model, at adapt_delta 0.99 ===\n")
print(sds_ad %>% filter(group == "region") %>% mutate(sd = round(sd, 4)) %>%
        as.data.frame(), row.names = FALSE)
cat("\n=== Random-effect standard deviations under three scale priors ===\n")
print(sds %>% mutate(sd = round(sd, 4)) %>% as.data.frame(), row.names = FALSE)

cat("\n=== Sensitivity of the regional slope SDs to the prior ===\n")
piv <- sds %>% filter(group == "region", !is.na(term), term != "(Intercept)") %>%
  tidyr::pivot_wider(names_from = spec, values_from = sd)
print(as.data.frame(piv %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)
rng <- sds %>% filter(group == "region", term != "(Intercept)") %>%
  group_by(term) %>% summarise(lo = min(sd), hi = max(sd), .groups = "drop") %>%
  mutate(ratio = hi / lo)
cat("\nratio of largest to smallest posterior SD across priors:\n")
print(as.data.frame(rng %>% mutate(across(where(is.numeric), ~ round(.x, 3)))),
      row.names = FALSE)
cat("A ratio near 1 means the data rather than the prior set the scale.\n")

# ---- region-specific slopes with credible intervals --------
post <- as.data.frame(m_main)
slope_terms <- c("IMM_wn_s", "SAN_wn_s", "EDU_wn_s")
regs <- sort(unique(as.character(dat$region)))
rs <- lapply(slope_terms, function(tm) {
  gl <- post[[tm]]
  bind_rows(lapply(regs, function(rg) {
    cn <- sprintf("b[%s region:%s]", tm, gsub(" ", "_", rg))
    if (!cn %in% names(post)) {
      cand <- grep(paste0("^b\\[", tm, " region:"), names(post), value = TRUE)
      hit <- cand[gsub("_", " ", sub(".*region:", "", sub("\\]$", "", cand))) == rg]
      if (!length(hit)) return(NULL)
      cn <- hit[1]
    }
    tot <- gl + post[[cn]]
    data.frame(term = tm, region = rg, mean = mean(tot),
               q025 = quantile(tot, 0.025), q975 = quantile(tot, 0.975),
               dev_mean = mean(post[[cn]]),
               dev_excl0 = quantile(post[[cn]], 0.025) > 0 |
                           quantile(post[[cn]], 0.975) < 0,
               row.names = NULL)
  }))
}) %>% bind_rows()
cat("\n=== Region-specific within-country slopes (standardised) ===\n")
print(as.data.frame(rs %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)
cat(sprintf("\nregional deviations whose 95%% interval excludes zero: %d of %d\n",
            sum(rs$dev_excl0), nrow(rs)))
cat("With country intercepts on a highly persistent outcome the posterior\n")
cat("R^2 is close to one by construction; Supplementary Section S2 reports\n")
cat("these regional deviations instead.\n")

# ---- fixed effects of the main fit, with credible intervals --
# quoted in Supplementary Section S2, drawn in Supplementary Figure S8
fx <- grep("^(year_c|is_covid_year|.*_(bw|wn)_s)$", names(post), value = TRUE)
fixef_tbl <- data.frame(term = fx,
                        mean = sapply(post[fx], mean),
                        sd   = sapply(post[fx], sd),
                        q025 = sapply(post[fx], quantile, 0.025),
                        q975 = sapply(post[fx], quantile, 0.975),
                        row.names = NULL) %>%
  mutate(signif_95 = q025 > 0 | q975 < 0)
cat("\n=== Fixed effects, main fit (standardised predictors) ===\n")
print(as.data.frame(fixef_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)

# ---- regional slope SDs on the scale of the drivers ---------
# The model is fitted on standardised predictors; dividing by the
# standard deviation of each within-country component gives years
# of life expectancy per unit of the driver.
sd_wn <- c(IMM_wn_s = sd(dat$IMM_wn), SAN_wn_s = sd(dat$SAN_wn),
           EDU_wn_s = sd(dat$EDU_wn))
sd_std <- sds %>% filter(spec == "main decov(1,1,1,1)", group == "region",
                         term %in% slope_terms)
scale_tbl <- data.frame(term = slope_terms,
                        sd_within_driver = sd_wn[slope_terms],
                        slope_sd_std = sd_std$sd[match(slope_terms, sd_std$term)],
                        mean_std = colMeans(post[slope_terms]),
                        row.names = NULL) %>%
  mutate(slope_sd_unit = slope_sd_std / sd_within_driver,
         mean_unit     = mean_std / sd_within_driver,
         sd_over_mean  = slope_sd_std / abs(mean_std))
cat("\n=== Regional slope SD against the pooled within mean ===\n")
print(as.data.frame(scale_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)

# ---- posterior R^2 and PSIS-LOO of the main fit --------------
r2_main <- bayes_R2(m_main)
cat(sprintf("\nposterior R^2 of the main fit: mean %.3f\n", mean(r2_main)))
loo_main <- loo(m_main, cores = 2)
print(loo_main)
pk <- loo_main$diagnostics$pareto_k
cat(sprintf("observations with Pareto k above 0.7: %d of %d (%.1f%%)\n",
            sum(pk > 0.7), length(pk), 100 * mean(pk > 0.7)))

write.csv(fixef_tbl, file.path(OUT, "24_fixed_effects.csv"), row.names = FALSE)
write.csv(scale_tbl, file.path(OUT, "24_slope_scale.csv"), row.names = FALSE)
write.csv(as.data.frame(diag_tbl), file.path(OUT, "24_convergence.csv"), row.names = FALSE)
write.csv(sds, file.path(OUT, "24_prior_sensitivity.csv"), row.names = FALSE)
write.csv(rs,  file.path(OUT, "24_region_slopes.csv"), row.names = FALSE)
write.csv(ad_tbl, file.path(OUT, "24_adapt_delta.csv"), row.names = FALSE)
cat("\nSaved: 24_convergence.csv, 24_prior_sensitivity.csv,\n")
cat("24_region_slopes.csv, 24_adapt_delta.csv, 24_fixed_effects.csv,\n")
cat("24_slope_scale.csv\n")

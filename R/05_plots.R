# 05: Supplementary Figures S1 to S4
#
#   S1  share of variance of the shape components
#   S2  loadings of the retained components
#   S3  paired difference in VIP, level minus shape
#   S4  out-of-region R2 by held-out region
#
# Reads  output/shape_pca_variance.csv, shape_pca_loadings.csv, pls/
# Writes figures/FigureS1 to FigureS4

source("00_common.R")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

# The figures carry no legend, title or annotation: the captions in the
# Supplementary describe the symbols, and variable codes are replaced
# by their names.
DRIVER_NAMES <- c(log_GNI = "Log income per head", EDU = "Secondary enrolment",
                  IMM = "Immunisation", SAN = "Basic sanitation",
                  PUB = "Public share", OOP = "Out-of-pocket share",
                  PRIV_NONOOP = "Private non-out-of-pocket share",
                  EXT = "External share")
FEATURE_NAMES <- c(trend_slope = "Trend slope",
                   resid_sd = "Residual standard deviation",
                   worst_drop = "Worst drop", best_jump = "Best jump",
                   n_reversal = "Reversal count", pelt_n = "Changepoint count",
                   slope_range = "Segmental slope range",
                   slope_max = "Maximum segmental slope")


# S1. variance by component --------------------------------------
pca_var <- read.csv(file.path(OUT, "shape_pca_variance.csv"))
pca_var$PC <- factor(pca_var$PC, levels = pca_var$PC)
pca_var$retained <- pca_var$cum_var <= 0.85

p_var <- ggplot(pca_var, aes(x = PC, y = var_expl, fill = retained)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_line(aes(y = cum_var, group = 1), color = "black", linewidth = 0.7) +
  geom_point(aes(y = cum_var), color = "black", size = 2.2, shape = 17) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("TRUE" = "grey25", "FALSE" = "grey80"),
                    guide = "none") +
  scale_x_discrete(labels = function(x) sub("^PC", "", x)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Principal component", y = "Explained variance") +
  theme_minimal(base_size = 12)

save_fig("FigureS1", p_var, width = 8, height = 4.5)

# S2. loadings ---------------------------------------------------
load <- read.csv(file.path(OUT, "shape_pca_loadings.csv"),
                 stringsAsFactors = FALSE)
names(load)[1] <- "feature"
load_long <- pivot_longer(load, -feature, names_to = "PC", values_to = "loading")
stopifnot(all(load_long$feature %in% names(FEATURE_NAMES)))
load_long$feature <- FEATURE_NAMES[load_long$feature]
# panels A, B and C are the first, second and third components
load_long$panel <- LETTERS[match(load_long$PC, sort(unique(load_long$PC)))]

p_load <- ggplot(load_long, aes(x = loading, y = feature, fill = loading > 0)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey50") +
  facet_wrap(~ panel, ncol = 3) +
  scale_fill_manual(values = c("TRUE" = "grey30", "FALSE" = "white"),
                    guide = "none") +
  labs(x = "Loading", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12, hjust = 0))

save_fig("FigureS2", p_load, width = 9, height = 4.5)

# S3. paired difference in VIP -----------------------------------
vip <- read.csv(file.path(OUT, "pls", "VIP_predictive_paired.csv"),
                stringsAsFactors = FALSE)
vip$regressor <- sub("_init$","", vip$regressor)
stopifnot(all(vip$regressor %in% names(DRIVER_NAMES)))
vip$regressor <- DRIVER_NAMES[vip$regressor]

vip_d <- vip %>% transmute(regressor, delta_mean, delta_q025, delta_q975,
                           diff_signif)
vip_d$regressor <- factor(vip_d$regressor,
                          levels = rev(vip$regressor[order(vip$delta_mean)]))

p_delta <- ggplot(vip_d, aes(x = delta_mean, y = regressor,
                             shape = diff_signif, linetype = diff_signif)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_errorbar(aes(xmin = delta_q025, xmax = delta_q975), width = 0.2,
                colour = "black") +
  geom_point(size = 3, colour = "black", fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                     labels = c("FALSE" = "not significant", "TRUE" = "p < 0.05")) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "22"),
                        labels = c("FALSE" = "not significant", "TRUE" = "p < 0.05")) +
  labs(x = "Difference in variable importance in projection, level minus shape",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

save_fig("FigureS3", p_delta,
       width = 8, height = 5)

# S4. regional holdout -------------------------------------------
ho <- read.csv(file.path(OUT, "pls", "holdout_by_region.csv"),
               stringsAsFactors = FALSE)
ho_long <- pivot_longer(ho, cols = c(R2_shift, R2_shape),
                        names_to = "target", values_to = "R2") %>%
  mutate(target = recode(target,
                         R2_shift = "LE shift",
                         R2_shape = "Shape PCs"))

p_ho <- ggplot(ho_long, aes(x = R2, y = reorder(region, R2), fill = target)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_col(position = "dodge", colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("LE shift" = "grey30",
                               "Shape PCs" = "grey85")) +
  labs(x = "Out-of-region R²", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

save_fig("FigureS4", p_ho, width = 9, height = 5)

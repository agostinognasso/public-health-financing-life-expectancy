# 10: Supplementary Figures S5 and S6
#
#   S5  initial-window sensitivity and the ten most influential
#       countries (script 09)
#   S6  VIP for level and shape with the ILR coordinates (script 08)
#
# Reads  output/sensitivity_window.csv, influence_pls_loo.csv,
#        pls/VIP_ilr_paired.csv
# Writes figures/FigureS5, FigureS6

source("00_common.R")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(forcats); library(patchwork)
})


theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey92"),
      axis.line.x = element_line(linewidth = 0.4, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.3),
      legend.position = "none"
    )
}
# No legends, titles or annotations inside the figures: the captions in
# the Supplementary describe the symbols.
# Colour palette, checked for colour-vision deficiency; every colour is
# paired with a distinct symbol or line type, so the figures also read in
# greyscale.
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"

# S6. VIP with the ILR coordinates -------------------------------

vip <- read.csv(file.path(OUT, "pls", "VIP_ilr_paired.csv"),
                stringsAsFactors = FALSE)
vip$regressor <- sub("_init$","", vip$regressor)
vip$regressor <- recode(vip$regressor,
  ILR1_PUBvsRest    = "Balance: public vs non-public",
  ILR2_DomPrivVsEXT = "Balance: domestic private vs external",
  ILR3_OOPvsPNO    = "Balance: out-of-pocket vs prepaid",
  log_GNI          = "Log income per head",
  EDU              = "Secondary enrolment",
  IMM              = "Immunisation",
  SAN              = "Basic sanitation"
)
vip_long <- bind_rows(
  vip %>% transmute(regressor, target = "Level (LE shift)",
                    mean = VIP_shift_mean,
                    lo = VIP_shift_q025, hi = VIP_shift_q975),
  vip %>% transmute(regressor, target = "Shape (PCs)",
                    mean = VIP_shape_mean,
                    lo = VIP_shape_q025, hi = VIP_shape_q975)
)
ord <- vip$regressor[order(vip$VIP_shift_mean + vip$VIP_shape_mean,
                           decreasing = TRUE)]
vip_long$regressor <- factor(vip_long$regressor, levels = rev(ord))

p3 <- ggplot(vip_long, aes(x = mean, y = regressor,
                           shape = target, linetype = target,
                           colour = target)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2,
                position = position_dodge(0.5)) +
  geom_point(size = 2.8, position = position_dodge(0.5),
             fill = "white", stroke = 0.8) +
  scale_colour_manual(values = c(BLUE, ORANGE), name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = "Variable importance in projection", y = NULL) +
  theme_pub(11)
save_fig("FigureS6", p3,
       width = 8.5, height = 5)

# S5. initial window and influence -------------------------------

sens <- read.csv(file.path(OUT, "sensitivity_window.csv"),
                 stringsAsFactors = FALSE)
sens$window <- sens$years                 # label the windows by their years
sens_long <- pivot_longer(sens,
                          cols = c(R2_shift_max, R2_PC2_max),
                          names_to = "target", values_to = "R2") %>%
  mutate(target = recode(target,
                         R2_shift_max = "LE shift",
                         R2_PC2_max   = "Shape PC2"))

p4a <- ggplot(sens_long, aes(x = window, y = R2,
                             linetype = target, shape = target,
                             colour = target, group = target)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.8, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = c(BLUE, ORANGE), name = NULL) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_y_continuous(limits = c(0.25, 0.45)) +
  labs(x = "Initial window", y = "Cross-validated R²") +
  theme_pub(10)

inf <- read.csv(file.path(OUT, "influence_pls_loo.csv"),
                stringsAsFactors = FALSE) %>% head(10)
inf$country <- factor(inf$country, levels = rev(inf$country))

p4b <- ggplot(inf, aes(x = dR2, y = country, fill = dR2 > 0)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("TRUE" = BLUE, "FALSE" = ORANGE),
                    guide = "none") +
  labs(x = "ΔR² (removing this country)", y = NULL) +
  theme_pub(10)

p4 <- (p4a + p4b) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))
save_fig("FigureS5", p4,
       width = 11, height = 4.5)

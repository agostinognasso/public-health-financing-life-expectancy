# 32: information criteria for the raw-share and ILR mixed models,
#     and the correlations among the financing shares (Section 3.2 and
#     Supplementary Section S5)
#
# A REML likelihood is computed on contrasts that depend on the
# fixed-effect design, so REML criteria compare only models with the
# same fixed effects, which the raw-share model of script 04 and the
# ILR model of script 08 do not have; and three ILR coordinates
# against two raw shares would give the ILR block one more degree of
# freedom. The models are refitted by maximum likelihood on the same
# observations and controls, at matched degrees of freedom.
#
# Reads  output/panel_lmer_model.rds, panel_ilr_model.rds; data/
#        through build_panel()
# Writes output/32_simplex_ic.csv

source("00_common.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(lme4))
options(warn = -1)

AN <- OUT
mb <- readRDS(file.path(AN, "panel_lmer_model.rds"))
mi <- readRDS(file.path(AN, "panel_ilr_model.rds"))

cat("=== REML fits of scripts 04 and 08 ===\n")
cat(sprintf("REML fits as saved: base AIC %.2f, ILR AIC %.2f, gap %+.2f\n",
            AIC(mb), AIC(mi), AIC(mb) - AIC(mi)))
cat(sprintf("isREML: base %s, ILR %s; fixed-effect terms: %d vs %d\n",
            isREML(mb), isREML(mi), length(fixef(mb)), length(fixef(mi))))
cat("The two models differ in their fixed effects, so these REML\n")
cat("criteria are not comparable. Refitting by maximum likelihood:\n")
bm <- refitML(mb); im <- refitML(mi)
cat(sprintf("ML fits: base AIC %.2f, ILR AIC %.2f, gap %+.2f\n",
            AIC(bm), AIC(im), AIC(bm) - AIC(im)))
cat(sprintf("ML fits: base BIC %.2f, ILR BIC %.2f, gap %+.2f\n",
            BIC(bm), BIC(im), BIC(bm) - BIC(im)))

# ---- matched-df comparison on identical data ----------------
fb <- model.frame(mb); fi <- model.frame(mi)
stopifnot(nrow(fb) == nrow(fi), isTRUE(all.equal(fb$LE, fi$LE)))
D <- cbind(fb, fi[, c("ILR1_bw","ILR2_bw","ILR3_bw",
                      "ILR1_wn","ILR2_wn","ILR3_wn")])
ctl <- paste("year_c + is_covid_year +",
             "log_GNI_y_bw + EDU_y_bw + IMM_y_bw + SAN_y_bw +",
             "log_GNI_y_wn + EDU_y_wn + IMM_y_wn + SAN_y_wn")
specs <- list(
  "raw shares PUB + OOP"    = "OOP_y_bw+PUB_y_bw+OOP_y_wn+PUB_y_wn",
  "log-ratio ILR1 + ILR2"   = "ILR1_bw+ILR2_bw+ILR1_wn+ILR2_wn",
  "raw share PUB only"      = "PUB_y_bw+PUB_y_wn",
  "log-ratio ILR1 only"     = "ILR1_bw+ILR1_wn",
  "log-ratio ILR1+ILR2+ILR3"= "ILR1_bw+ILR2_bw+ILR3_bw+ILR1_wn+ILR2_wn+ILR3_wn")

cat(sprintf("\n=== Maximum likelihood, n = %d, identical controls ===\n", nrow(D)))
tab <- do.call(rbind, lapply(names(specs), function(nm) {
  m <- lmer(as.formula(paste("LE ~", ctl, "+", specs[[nm]], "+ (1|iso3c)")),
            data = D, REML = FALSE)
  data.frame(spec = nm, fin_terms = (attr(logLik(m),"df") - 15) ,
             df = attr(logLik(m),"df"),
             AIC = AIC(m), BIC = BIC(m), row.names = NULL)
}))
tab$dAIC <- tab$AIC - min(tab$AIC); tab$dBIC <- tab$BIC - min(tab$BIC)
print(tab %>% mutate(across(where(is.numeric), ~ round(.x, 2))))

g2 <- tab$AIC[tab$spec=="raw shares PUB + OOP"] -
      tab$AIC[tab$spec=="log-ratio ILR1 + ILR2"]
g1 <- tab$AIC[tab$spec=="raw share PUB only"] -
      tab$AIC[tab$spec=="log-ratio ILR1 only"]
cat(sprintf("\nAt matched degrees of freedom the log-ratio coordinates\n"))
cat(sprintf("improve the Akaike criterion by %.2f points on two financing\n", g2))
cat(sprintf("terms and by %.2f on one.\n", g1))

# ---- omitting the out-of-pocket share ----------------------
cat("\n=== Two-way FE on the panel, with and without out-of-pocket ===\n")
dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
pd  <- pdata.frame(dat, index = c("iso3c","year"))
CT  <- "log_GNI + EDU_GROSS + IMM_DPT + SANITATION"
f_with <- plm(as.formula(paste("LE ~", CT, "+ PUB + OOP_share")), data = pd,
              model = "within", effect = "twoways")
f_wo   <- plm(as.formula(paste("LE ~", CT, "+ PUB")), data = pd,
              model = "within", effect = "twoways")
V1 <- plm::vcovHC(f_with, cluster="group", method="arellano")
V0 <- plm::vcovHC(f_wo,   cluster="group", method="arellano")
cat(sprintf("PUB with OOP in the model: %+.5f (SE %.5f)\n",
            coef(f_with)["PUB"], sqrt(V1["PUB","PUB"])))
cat(sprintf("PUB with OOP omitted:      %+.5f (SE %.5f)\n",
            coef(f_wo)["PUB"], sqrt(V0["PUB","PUB"])))
cat(sprintf("shift = %+.5f, which is %.2f standard errors\n",
            coef(f_wo)["PUB"] - coef(f_with)["PUB"],
            abs(coef(f_wo)["PUB"] - coef(f_with)["PUB"])/sqrt(V1["PUB","PUB"])))

cat("\n=== Correlations among the financing shares ===\n")
cat(sprintf("corr(PUB, OOP) = %+.3f\n",
            cor(dat$PUB, dat$OOP_share, use="complete.obs")))
oop_priv <- dat$OOP_share / pmax(100 - dat$PUB, 0.05)
cat(sprintf("out-of-pocket as a share of non-public financing (100 - PUB): mean %.2f, median %.2f\n",
            mean(oop_priv, na.rm=TRUE), median(oop_priv, na.rm=TRUE)))
# domestic private financing, PRIV = OOP + PNO, as in Section 3.2
pp_ <- dat$PRIV > 0
oop_dom <- dat$OOP_share[pp_] / dat$PRIV[pp_]
cat(sprintf("out-of-pocket as a share of domestic private financing (PRIV): mean %.2f, median %.2f\n",
            mean(pmin(oop_dom, 1)), median(pmin(oop_dom, 1))))
cat(sprintf("ratio of mean out-of-pocket to mean voluntary prepaid share: %.2f\n",
            mean(dat$OOP_share) / mean(pmax(dat$PRIV - dat$OOP_share, 0))))
cat(sprintf("corr(OOP, PRIV) = %+.3f\n", cor(dat$OOP_share, dat$PRIV)))
cat(sprintf("corr(OOP, log(PRIV/PUB)) = %+.3f\n",
            cor(dat$OOP_share[pp_], log(dat$PRIV[pp_] / dat$PUB[pp_]))))

write.csv(tab, file.path(OUT, "32_simplex_ic.csv"), row.names = FALSE)
cat("\nSaved: 32_simplex_ic.csv\n")

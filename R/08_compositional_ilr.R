# 08: ILR coordinates of the financing shares
#
# The four financing shares (PUB, OOP, PNO, EXT) sum to 100 and are
# treated as a composition. Isometric log-ratio coordinates follow the
# sequential binary partition
#   ILR1  {PUB} | {OOP, PNO, EXT}   public against non-public
#   ILR2  {OOP, PNO} | {EXT}        domestic private against external
#   ILR3  {OOP} | {PNO}             out-of-pocket against voluntary
# The coordinates replace the raw shares (1) in the predictive PLS with
# the paired VIP bootstrap, and (2) in the between/within mixed model of
# script 04. Parts below 0.05 are set to 0.05 before the log-ratios.
#
# Reads  output/panel_with_ext.rds, driver_initial.rds, le_targets.rds,
#        shape_pcs.rds; data/df_le_long.rds
# Writes output/pls/VIP_ilr_paired.csv, panel_ilr_coefs.csv,
#        panel_ilr_model.rds

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(lmerTest)
  library(compositions)
})

panel <- readRDS(file.path(OUT, "panel_with_ext.rds"))
le_long <- readRDS(file.path(DATA, "df_le_long.rds"))

dat <- panel %>%
  inner_join(le_long %>% select(iso3c, year, LE),
             by = c("iso3c","year")) %>%
  filter(!is.na(LE), !is.na(log_GNI), !is.na(EDU_GROSS),
         !is.na(PUB), !is.na(OOP_share),
         !is.na(IMM_DPT), !is.na(SANITATION))

# four-part composition; parts below eps set to eps, then reclosed
eps <- 0.05
comp_mat <- dat %>%
  mutate(PUB = pmax(PUB, eps),
         OOP = pmax(OOP_share, eps),
         PNO = pmax(PRIV_NONOOP, eps),
         EXT = pmax(EXT_share, eps)) %>%
  mutate(tot = PUB + OOP + PNO + EXT,
         PUB = PUB / tot * 100,
         OOP = OOP / tot * 100,
         PNO = PNO / tot * 100,
         EXT = EXT / tot * 100) %>%
  select(PUB, OOP, PNO, EXT)

# sequential binary partition, columns PUB, OOP, PNO, EXT:
# +1 numerator, -1 denominator, 0 not involved
sbp <- matrix(c(
  # PUB OOP PNO EXT
    +1, -1, -1, -1,   # ILR1: PUB vs {OOP,PNO,EXT}
     0, +1, +1, -1,   # ILR2: {OOP,PNO} vs EXT
     0, +1, -1,  0    # ILR3: OOP vs PNO
), nrow = 3, byrow = TRUE)
colnames(sbp) <- c("PUB","OOP","PNO","EXT")
rownames(sbp) <- c("ILR1_PUBvsRest","ILR2_DomPrivVsEXT","ILR3_OOPvsPNO")

ilr_from_sbp <- function(comp, sbp) {
  comp <- as.matrix(comp)
  apply(sbp, 1, function(row) {
    pos <- which(row > 0); neg <- which(row < 0)
    r <- length(pos); s <- length(neg)
    factor_ <- sqrt(r*s / (r + s))
    log_pos <- rowMeans(log(comp[, pos, drop = FALSE]))
    log_neg <- rowMeans(log(comp[, neg, drop = FALSE]))
    factor_ * (log_pos - log_neg)
  })
}

ilr <- as.data.frame(ilr_from_sbp(comp_mat, sbp))
dat <- bind_cols(dat, ilr)

cat("ILR coordinates, first rows:\n")
print(head(round(ilr, 3)))
cat("\nCorrelations between ILR coords:\n")
print(round(cor(ilr), 3))

# 1. predictive PLS with the ILR coordinates ---------------------

ds_init <- readRDS(file.path(OUT, "driver_initial.rds"))
le_tg   <- readRDS(file.path(OUT, "le_targets.rds"))
pcs     <- readRDS(file.path(OUT, "shape_pcs.rds"))

# ILR means over the initial window, 2000-2005
init <- dat %>% filter(year >= 2000, year <= 2005) %>%
  group_by(iso3c) %>%
  summarise(across(c("ILR1_PUBvsRest","ILR2_DomPrivVsEXT","ILR3_OOPvsPNO"),
                   ~ mean(.x, na.rm=TRUE)),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.finite(.x), .x, NA_real_)))

# the other drivers, as in script 03
nc_init <- ds_init %>%
  select(iso3c, country,
         log_GNI_init, EDU_init, IMM_init, SAN_init)

cs <- nc_init %>%
  inner_join(init,  by = "iso3c") %>%
  inner_join(le_tg, by = "iso3c") %>%
  inner_join(pcs,   by = "iso3c") %>%
  na.omit()

cat("\nCross-section with ILR coordinates: ", nrow(cs), " countries\n")

X_vars <- c("log_GNI_init","EDU_init","IMM_init","SAN_init",
            "ILR1_PUBvsRest","ILR2_DomPrivVsEXT","ILR3_OOPvsPNO")
shape_targets <- c("PC1","PC2","PC3")

X <- as.matrix(cs[, X_vars])
y_shift <- cs$LE_shift
Y_pcs   <- as.matrix(cs[, shape_targets])

Xs <- scale(X); ys <- scale(y_shift); Ys <- scale(Y_pcs)

library(pls)
m_shift <- plsr(ys ~ Xs, ncomp = min(6, ncol(Xs)),
                validation = "LOO", method = "oscorespls")
ncomp_shift <- selectNcomp(m_shift, method = "onesigma", plot = FALSE)
if (is.na(ncomp_shift) || ncomp_shift < 1) ncomp_shift <- 2

m_shape <- plsr(Ys ~ Xs, ncomp = min(6, ncol(Xs)),
                validation = "LOO", method = "oscorespls")
rmsep_arr <- RMSEP(m_shape, estimate = "CV", intercept = FALSE)$val
rmsep_tot <- apply(rmsep_arr, 3, sum)
ncomp_shape <- min(which(rmsep_tot <= min(rmsep_tot) * 1.05))
if (ncomp_shape < 1) ncomp_shape <- 2

cat("ncomp PLS-shift (ILR):", ncomp_shift, "\n")
cat("ncomp PLS-shape (ILR):", ncomp_shape, "\n")

R2_shift_ilr <- R2(m_shift, estimate="CV",
                   ncomp = 1:ncomp_shift, intercept=FALSE)$val
cat("R2 CV LE_shift (ILR): ",
    paste(round(R2_shift_ilr[1,1,], 3), collapse=", "), "\n")

R2_shape_ilr <- R2(m_shape, estimate="CV",
                   ncomp = 1:ncomp_shape, intercept=FALSE)$val
cat("R2 CV shape PCs (ILR):\n")
for (i in seq_along(shape_targets))
  cat(sprintf("  %s: %s\n", shape_targets[i],
              paste(round(R2_shape_ilr[1,i,], 3), collapse=", ")))

m_shift_fit <- plsr(ys ~ Xs, ncomp = ncomp_shift, method = "oscorespls")
m_shape_fit <- plsr(Ys ~ Xs, ncomp = ncomp_shape, method = "oscorespls")

VIP <- function(model) {
  W  <- model$loading.weights; H  <- ncol(W)
  q  <- model$Yloadings; T  <- model$scores
  if (is.null(dim(q))) q <- matrix(q, ncol = H)
  SS <- sapply(1:H, function(h) sum((q[, h]^2)) * sum(T[, h]^2))
  Wnorm2 <- colSums(W^2)
  v <- sapply(1:nrow(W), function(j)
    sqrt(nrow(W) * sum(SS * (W[j, ]/sqrt(Wnorm2))^2) / sum(SS)))
  names(v) <- rownames(W); v
}

vip_shift <- VIP(m_shift_fit)
vip_shape <- VIP(m_shape_fit)

# paired bootstrap, B = 500
set.seed(42)
B <- 500; n <- nrow(cs)
boot_shift <- matrix(NA, B, length(X_vars), dimnames=list(NULL,X_vars))
boot_shape <- matrix(NA, B, length(X_vars), dimnames=list(NULL,X_vars))
for (b in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  Xb <- scale(X[idx,]); ysb <- scale(y_shift[idx]); Ysb <- scale(Y_pcs[idx,])
  tryCatch({
    boot_shift[b,] <- VIP(plsr(ysb ~ Xb, ncomp=ncomp_shift, method="oscorespls"))[X_vars]
    boot_shape[b,] <- VIP(plsr(Ysb ~ Xb, ncomp=ncomp_shape, method="oscorespls"))[X_vars]
  }, error=function(e) NULL)
}

vip_ilr_tbl <- data.frame(
  regressor = X_vars,
  VIP_shift_mean = round(colMeans(boot_shift, na.rm=TRUE), 3),
  VIP_shift_q025 = round(apply(boot_shift,2,quantile,.025,na.rm=TRUE),3),
  VIP_shift_q975 = round(apply(boot_shift,2,quantile,.975,na.rm=TRUE),3),
  VIP_shape_mean = round(colMeans(boot_shape, na.rm=TRUE), 3),
  VIP_shape_q025 = round(apply(boot_shape,2,quantile,.025,na.rm=TRUE),3),
  VIP_shape_q975 = round(apply(boot_shape,2,quantile,.975,na.rm=TRUE),3)
) %>% mutate(
  stable_shift = VIP_shift_q025 > 1,
  stable_shape = VIP_shape_q025 > 1
)

cat("\n=== VIP with ILR coordinates (bootstrap intervals) ===\n")
print(vip_ilr_tbl)

write.csv(vip_ilr_tbl, file.path(OUT, "pls", "VIP_ilr_paired.csv"), row.names = FALSE)

# 2. mixed model with the ILR coordinates ------------------------

dat <- dat %>%
  group_by(iso3c) %>%
  mutate(
    ILR1_bw = mean(ILR1_PUBvsRest, na.rm=TRUE),
    ILR1_wn = ILR1_PUBvsRest - ILR1_bw,
    ILR2_bw = mean(ILR2_DomPrivVsEXT, na.rm=TRUE),
    ILR2_wn = ILR2_DomPrivVsEXT - ILR2_bw,
    ILR3_bw = mean(ILR3_OOPvsPNO, na.rm=TRUE),
    ILR3_wn = ILR3_OOPvsPNO - ILR3_bw,
    log_GNI_bw = mean(log_GNI, na.rm=TRUE),
    log_GNI_wn = log_GNI - log_GNI_bw,
    EDU_bw     = mean(EDU_GROSS, na.rm=TRUE),
    EDU_wn     = EDU_GROSS - EDU_bw,
    IMM_bw     = mean(IMM_DPT, na.rm=TRUE),
    IMM_wn     = IMM_DPT - IMM_bw,
    SAN_bw     = mean(SANITATION, na.rm=TRUE),
    SAN_wn     = SANITATION - SAN_bw
  ) %>% ungroup() %>%
  mutate(year_c = year - 2012)

f_ilr <- LE ~ year_c + is_covid_year +
    log_GNI_bw + EDU_bw + IMM_bw + SAN_bw +
    ILR1_bw + ILR2_bw + ILR3_bw +
    log_GNI_wn + EDU_wn + IMM_wn + SAN_wn +
    ILR1_wn + ILR2_wn + ILR3_wn +
    (1 | iso3c)

m_ilr <- lmer(f_ilr, data = dat, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))

co_ilr <- as.data.frame(coef(summary(m_ilr)))
co_ilr$term <- rownames(co_ilr); rownames(co_ilr) <- NULL

panel_ilr <- co_ilr %>%
  filter(grepl("_(bw|wn)$", term)) %>%
  mutate(type = ifelse(grepl("_bw$", term), "between","within"),
         var  = sub("_(bw|wn)$","", term)) %>%
  select(var, type, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`) %>%
  rename(beta = Estimate, SE = `Std. Error`,
         t = `t value`, p = `Pr(>|t|)`) %>%
  mutate(across(c(beta,SE,t), ~ round(.x,4)),
         p = round(p,4),
         sig = ifelse(p<0.001,"***",
               ifelse(p<0.01,"**",
               ifelse(p<0.05,"*",
               ifelse(p<0.1,".","")))))

cat("\n=== Mixed model with the ILR coordinates ===\n")
print(panel_ilr %>% arrange(var, type))

r2_nakagawa <- function(model) {
  v_fix <- var(predict(model, re.form = NA))
  vc <- as.data.frame(VarCorr(model))
  v_re <- sum(vc$vcov[vc$grp != "Residual"])
  v_res <- vc$vcov[vc$grp == "Residual"]
  list(marginal = v_fix / (v_fix + v_re + v_res),
       conditional = (v_fix + v_re) / (v_fix + v_re + v_res))
}
r2 <- r2_nakagawa(m_ilr)
cat(sprintf("\nR2 marginal (ILR): %.3f\n", r2$marginal))
cat(sprintf("R2 conditional (ILR): %.3f\n", r2$conditional))

write.csv(panel_ilr, file.path(OUT, "panel_ilr_coefs.csv"), row.names = FALSE)
saveRDS(m_ilr, file.path(OUT, "panel_ilr_model.rds"))

cat("\nSaved: pls/VIP_ilr_paired.csv, panel_ilr_coefs.csv, panel_ilr_model.rds\n")

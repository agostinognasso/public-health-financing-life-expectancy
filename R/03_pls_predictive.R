# 03: predictive PLS, regional holdout and paired VIP bootstrap
#
# Predictors: driver means over 2000-2005 (script 01).
# Targets:    LE_shift, and the shape components PC1-PC3 of 2005-2024
#             (script 02).
# Validation: leave-one-out CV, and holding out one World Bank region
#             at a time.
# Inference:  paired bootstrap, B = 500, with the same resample for the
#             level and the shape model, and a 95 per cent interval for
#             the difference VIP_level - VIP_shape.
#
# Reads  output/driver_initial.rds, le_targets.rds, shape_pcs.rds,
#        region_map.rds
# Writes output/pls/

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pls)
})

set.seed(42)
dir.create(file.path(OUT, "pls"), showWarnings = FALSE, recursive = TRUE)

# 1. cross-section -----------------------------------------------

ds_init <- readRDS(file.path(OUT, "driver_initial.rds"))
le_tg   <- readRDS(file.path(OUT, "le_targets.rds"))
pcs     <- readRDS(file.path(OUT, "shape_pcs.rds"))
reg     <- readRDS(file.path(OUT, "region_map.rds"))

X_vars <- c("log_GNI_init","EDU_init","OOP_init","PRIV_NONOOP_init",
            "PUB_init","EXT_init","IMM_init","SAN_init")
shape_targets <- c("PC1","PC2","PC3")

cs <- ds_init %>%
  inner_join(le_tg, by = "iso3c") %>%
  inner_join(pcs,   by = "iso3c") %>%
  inner_join(reg,   by = "iso3c") %>%
  select(iso3c, country, region, continent,
         all_of(X_vars), LE_shift, LE_level,
         all_of(shape_targets)) %>%
  na.omit()

cat("Predictive cross-section: ", nrow(cs), " countries\n")
cat("Countries by region:\n"); print(table(cs$region))

# 2. VIP ---------------------------------------------------------

VIP <- function(model) {
  W  <- model$loading.weights
  H  <- ncol(W)
  q  <- model$Yloadings
  T  <- model$scores
  if (is.null(dim(q))) q <- matrix(q, ncol = H)
  SS <- sapply(1:H, function(h) sum((q[, h]^2)) * sum(T[, h]^2))
  Wnorm2 <- colSums(W^2)
  sapply(1:nrow(W), function(j) {
    sqrt(nrow(W) * sum(SS * (W[j, ]/sqrt(Wnorm2))^2) / sum(SS))
  }) -> v
  names(v) <- rownames(W)
  v
}

# 3. PLS1 on LE_shift, PLS2 on the shape components ---------------

X <- as.matrix(cs[, X_vars])
y_shift <- cs$LE_shift
Y_pcs   <- as.matrix(cs[, shape_targets])

Xs <- scale(X); ys_shift <- scale(y_shift); Ys_pcs <- scale(Y_pcs)

# number of components by leave-one-out CV
m_shift <- plsr(ys_shift ~ Xs, ncomp = min(6, ncol(Xs)),
                validation = "LOO", method = "oscorespls")
ncomp_shift <- selectNcomp(m_shift, method = "onesigma", plot = FALSE)
if (is.na(ncomp_shift) || ncomp_shift < 1) ncomp_shift <- 2

# multivariate target: smallest number of components within 5 per cent
# of the minimum total RMSEP
m_shape <- plsr(Ys_pcs ~ Xs, ncomp = min(6, ncol(Xs)),
                validation = "LOO", method = "oscorespls")
rmsep_arr <- RMSEP(m_shape, estimate = "CV", intercept = FALSE)$val
rmsep_tot <- apply(rmsep_arr, 3, sum)
best <- which.min(rmsep_tot)
threshold <- rmsep_tot[best] * 1.05
ncomp_shape <- min(which(rmsep_tot <= threshold))
if (ncomp_shape < 1) ncomp_shape <- 2

cat("ncomp PLS-shift: ", ncomp_shift, "\n")
cat("ncomp PLS-shape: ", ncomp_shape, "\n")

# refit at the selected number of components
m_shift_fit <- plsr(ys_shift ~ Xs, ncomp = ncomp_shift, method = "oscorespls")
m_shape_fit <- plsr(Ys_pcs   ~ Xs, ncomp = ncomp_shape, method = "oscorespls")

vip_shift <- VIP(m_shift_fit)
vip_shape <- VIP(m_shape_fit)

# cross-validated R2
R2_shift_cv <- R2(m_shift, estimate="CV", ncomp = 1:ncomp_shift, intercept=FALSE)$val
R2_shape_cv <- R2(m_shape, estimate="CV", ncomp = 1:ncomp_shape, intercept=FALSE)$val

cat("\nR2 CV LE_shift (per comp): ",
    paste(round(R2_shift_cv[1,1,], 3), collapse=", "), "\n")
cat("R2 CV shape PC1 (per comp): ",
    paste(round(R2_shape_cv[1,"PC1",], 3), collapse=", "), "\n")
cat("R2 CV shape PC2 (per comp): ",
    paste(round(R2_shape_cv[1,"PC2",], 3), collapse=", "), "\n")
cat("R2 CV shape PC3 (per comp): ",
    paste(round(R2_shape_cv[1,"PC3",], 3), collapse=", "), "\n")

# X and Y loadings
write.csv(round(loadings(m_shift_fit)[,1:ncomp_shift, drop=FALSE], 3),
          file.path(OUT, "pls", "shift_loadings_X.csv"))
write.csv(round(loadings(m_shape_fit)[,1:ncomp_shape, drop=FALSE], 3),
          file.path(OUT, "pls", "shape_loadings_X.csv"))
write.csv(round(m_shape_fit$Yloadings[,1:ncomp_shape, drop=FALSE], 3),
          file.path(OUT, "pls", "shape_loadings_Y.csv"))

# 4. one region held out at a time -------------------------------

regions <- unique(cs$region)
regions <- regions[!is.na(regions) & regions != ""]
cat("\n--- Regional holdout: ", length(regions), " regions ---\n")

holdout_results <- data.frame()
for (rg in regions) {
  test_idx  <- which(cs$region == rg)
  train_idx <- setdiff(seq_len(nrow(cs)), test_idx)
  if (length(test_idx) < 3 || length(train_idx) < 30) next

  X_tr <- scale(X[train_idx, ])
  X_te <- scale(X[test_idx, ],
                center = attr(X_tr, "scaled:center"),
                scale  = attr(X_tr, "scaled:scale"))
  ys_tr <- scale(y_shift[train_idx])
  ys_mu <- attr(ys_tr, "scaled:center"); ys_sd <- attr(ys_tr, "scaled:scale")
  Yp_tr <- scale(Y_pcs[train_idx, ])
  Yp_mu <- attr(Yp_tr, "scaled:center"); Yp_sd <- attr(Yp_tr, "scaled:scale")

  m1 <- plsr(ys_tr ~ X_tr, ncomp = ncomp_shift, method = "oscorespls")
  m2 <- plsr(Yp_tr ~ X_tr, ncomp = ncomp_shape, method = "oscorespls")

  pred_shift <- drop(predict(m1, newdata = X_te, ncomp = ncomp_shift))*ys_sd + ys_mu
  pred_shape <- predict(m2, newdata = X_te, ncomp = ncomp_shape)[,,1]
  pred_shape <- sweep(pred_shape, 2, Yp_sd, "*")
  pred_shape <- sweep(pred_shape, 2, Yp_mu, "+")

  true_shift <- y_shift[test_idx]
  true_shape <- Y_pcs[test_idx, ]

  rmse_shift <- sqrt(mean((pred_shift - true_shift)^2))
  rmse_shape <- sqrt(mean((pred_shape - true_shape)^2))
  R2_shift   <- 1 - sum((pred_shift - true_shift)^2) /
                    sum((true_shift - mean(true_shift))^2)
  R2_shape   <- 1 - sum((pred_shape - true_shape)^2) /
                    sum((true_shape - mean(true_shape))^2)

  holdout_results <- rbind(holdout_results, data.frame(
    region    = rg,
    n_test    = length(test_idx),
    n_train   = length(train_idx),
    RMSE_shift= round(rmse_shift, 3),
    R2_shift  = round(R2_shift, 3),
    RMSE_shape= round(rmse_shape, 3),
    R2_shape  = round(R2_shape, 3)
  ))
}
print(holdout_results)
write.csv(holdout_results, file.path(OUT, "pls", "holdout_by_region.csv"), row.names=FALSE)

# 5. paired bootstrap on the VIP ---------------------------------

cat("\n--- Paired bootstrap (B = 500) ---\n")
B <- 500
n <- nrow(cs)

boot_vip_shift <- matrix(NA, B, length(X_vars),
                         dimnames = list(NULL, X_vars))
boot_vip_shape <- matrix(NA, B, length(X_vars),
                         dimnames = list(NULL, X_vars))

for (b in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  Xb  <- scale(X[idx, ])
  ysb <- scale(y_shift[idx])
  Ysb <- scale(Y_pcs[idx, ])
  tryCatch({
    mb1 <- plsr(ysb ~ Xb, ncomp = ncomp_shift, method = "oscorespls")
    boot_vip_shift[b, ] <- VIP(mb1)[X_vars]
    mb2 <- plsr(Ysb ~ Xb, ncomp = ncomp_shape, method = "oscorespls")
    boot_vip_shape[b, ] <- VIP(mb2)[X_vars]
  }, error = function(e) NULL)
}

delta <- boot_vip_shift - boot_vip_shape   # same resample for both

vip_test <- data.frame(
  regressor      = X_vars,
  VIP_shift_mean = round(colMeans(boot_vip_shift, na.rm=TRUE), 3),
  VIP_shift_q025 = round(apply(boot_vip_shift, 2, quantile, .025, na.rm=TRUE), 3),
  VIP_shift_q975 = round(apply(boot_vip_shift, 2, quantile, .975, na.rm=TRUE), 3),
  VIP_shape_mean = round(colMeans(boot_vip_shape, na.rm=TRUE), 3),
  VIP_shape_q025 = round(apply(boot_vip_shape, 2, quantile, .025, na.rm=TRUE), 3),
  VIP_shape_q975 = round(apply(boot_vip_shape, 2, quantile, .975, na.rm=TRUE), 3),
  delta_mean     = round(colMeans(delta, na.rm=TRUE), 3),
  delta_q025     = round(apply(delta, 2, quantile, .025, na.rm=TRUE), 3),
  delta_q975     = round(apply(delta, 2, quantile, .975, na.rm=TRUE), 3)
)

# difference significant when its 95 per cent interval excludes zero
vip_test$diff_signif <- with(vip_test, (delta_q025 > 0) | (delta_q975 < 0))
vip_test$role_shift <- with(vip_test, VIP_shift_q025 > 1)
vip_test$role_shape <- with(vip_test, VIP_shape_q025 > 1)

cat("\n=== VIP, level and shape, with the paired difference ===\n")
print(vip_test)

write.csv(vip_test, file.path(OUT, "pls", "VIP_predictive_paired.csv"), row.names = FALSE)

saveRDS(list(
  ncomp_shift = ncomp_shift,
  ncomp_shape = ncomp_shape,
  R2_shift_cv = R2_shift_cv,
  R2_shape_cv = R2_shape_cv,
  vip_test    = vip_test,
  holdout     = holdout_results
), file.path(OUT, "pls", "pls_results.rds"))

cat("\nSaved: output/pls/\n")

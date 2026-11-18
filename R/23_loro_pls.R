# 23: leave-one-region-out, decomposed, and the paired VIP
#     bootstrap under three resampling schemes
#
# An out-of-region R^2 measures predictions against the held-out
# region's own mean, which the fitted model does not see, so it
# falls below zero whenever the regional level is missed. The
# out-of-region MSE is split into a level and a dispersion part,
#    MSE = (mean(pred) - mean(true))^2 + MSE_centred,
# and the R^2 is also reported with the held-out region given its
# own mean (Supplementary Figure S10).
#
# The paired VIP bootstrap is run with B = 2000 under three schemes:
# countries independent, regions resampled as clusters, and a
# two-stage draw of regions and then of countries within them
# (Supplementary Figure S7).
#
# Reads  output/driver_initial.rds, le_targets.rds, shape_pcs.rds,
#        region_map.rds
# Writes output/23_loro_decomposed.csv, 23_vip_bootstrap.csv

source("00_common.R")
suppressPackageStartupMessages(library(pls))
set.seed(42)

X_vars <- c("log_GNI_init", "EDU_init", "OOP_init", "PRIV_NONOOP_init",
            "PUB_init", "EXT_init", "IMM_init", "SAN_init")
shape_targets <- c("PC1", "PC2", "PC3")

cs <- readRDS(file.path(OUT, "driver_initial.rds")) %>%
  inner_join(readRDS(file.path(OUT, "le_targets.rds")), by = "iso3c") %>%
  inner_join(readRDS(file.path(OUT, "shape_pcs.rds")),  by = "iso3c") %>%
  inner_join(readRDS(file.path(OUT, "region_map.rds")), by = "iso3c") %>%
  dplyr::select(iso3c, country, region, continent, all_of(X_vars),
                LE_shift, LE_level, all_of(shape_targets)) %>%
  na.omit()

cat("Cross-section:", nrow(cs), "countries\n")
print(table(cs$region))

X <- as.matrix(cs[, X_vars]); y_shift <- cs$LE_shift
Y_pcs <- as.matrix(cs[, shape_targets])

VIP <- function(model) {
  W <- model$loading.weights; H <- ncol(W)
  q <- model$Yloadings; T <- model$scores
  if (is.null(dim(q))) q <- matrix(q, ncol = H)
  SS <- sapply(1:H, function(h) sum(q[, h]^2) * sum(T[, h]^2))
  Wn <- colSums(W^2)
  v <- sapply(seq_len(nrow(W)), function(j)
    sqrt(nrow(W) * sum(SS * (W[j, ] / sqrt(Wn))^2) / sum(SS)))
  names(v) <- rownames(W); v
}

# component counts, as in script 03
m0s <- plsr(scale(y_shift) ~ scale(X), ncomp = 6, validation = "LOO",
            method = "oscorespls")
nc_shift <- selectNcomp(m0s, method = "onesigma", plot = FALSE)
if (is.na(nc_shift) || nc_shift < 1) nc_shift <- 2
m0p <- plsr(scale(Y_pcs) ~ scale(X), ncomp = 6, validation = "LOO",
            method = "oscorespls")
rt <- apply(RMSEP(m0p, estimate = "CV", intercept = FALSE)$val, 3, sum)
nc_shape <- min(which(rt <= rt[which.min(rt)] * 1.05))
cat(sprintf("\ncomponents: shift %d, shape %d\n", nc_shift, nc_shape))

# ------------------------------------------------------------
# 1. leave-one-region-out, decomposed
# ------------------------------------------------------------
regions <- sort(unique(cs$region[!is.na(cs$region) & cs$region != ""]))
loro <- lapply(regions, function(rg) {
  te <- which(cs$region == rg); tr <- setdiff(seq_len(nrow(cs)), te)
  if (length(te) < 3 || length(tr) < 30) return(NULL)

  Xtr <- scale(X[tr, ])
  Xte <- scale(X[te, ], center = attr(Xtr, "scaled:center"),
               scale = attr(Xtr, "scaled:scale"))
  ytr <- scale(y_shift[tr])
  mu <- attr(ytr, "scaled:center"); sdv <- attr(ytr, "scaled:scale")

  m <- plsr(ytr ~ Xtr, ncomp = nc_shift, method = "oscorespls")
  pred <- drop(predict(m, newdata = Xte, ncomp = nc_shift)) * sdv + mu
  true <- y_shift[te]

  ss_tot <- sum((true - mean(true))^2)
  mse    <- mean((pred - true)^2)
  bias   <- mean(pred) - mean(true)
  # centred MSE: remove each side's own mean, leaving only the
  # part of the error that is about ordering countries within the region
  mse_c  <- mean(((pred - mean(pred)) - (true - mean(true)))^2)

  data.frame(
    region = rg, n_test = length(te), n_train = length(tr),
    mean_true = mean(true), mean_pred = mean(pred),
    sd_true = sd(true), sd_pred = sd(pred),
    R2_total = 1 - mse * length(te) / ss_tot,
    level_bias = bias,
    pct_mse_from_level = 100 * bias^2 / mse,
    R2_given_region_mean = 1 - mse_c * length(te) / ss_tot,
    cor_pred_true = suppressWarnings(cor(pred, true)),
    row.names = NULL)
}) %>% bind_rows()

cat("\n=== Leave-one-region-out, decomposed ===\n")
print(loro %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

cat("\n=== Summary ===\n")
cat(sprintf("regions with negative total R2                 : %d of %d\n",
            sum(loro$R2_total < 0), nrow(loro)))
cat(sprintf("regions still negative after granting the mean : %d of %d\n",
            sum(loro$R2_given_region_mean < 0), nrow(loro)))
cat(sprintf("median share of out-of-region MSE due to level : %.0f%%\n",
            median(loro$pct_mse_from_level)))
cat(sprintf("mean correlation between predicted and actual  : %+.3f\n",
            mean(loro$cor_pred_true, na.rm = TRUE)))

# ------------------------------------------------------------
# 2. paired bootstrap on the VIP difference
# ------------------------------------------------------------
B <- 2000
vip_diff_once <- function(idx) {
  Xb <- X[idx, ]; yb <- y_shift[idx]; Yb <- Y_pcs[idx, ]
  if (sd(yb) < 1e-8) return(NULL)
  Xs <- scale(Xb)
  if (any(!is.finite(Xs))) return(NULL)
  m1 <- tryCatch(plsr(scale(yb) ~ Xs, ncomp = nc_shift, method = "oscorespls"),
                 error = function(e) NULL)
  m2 <- tryCatch(plsr(scale(Yb) ~ Xs, ncomp = nc_shape, method = "oscorespls"),
                 error = function(e) NULL)
  if (is.null(m1) || is.null(m2)) return(NULL)
  VIP(m1) - VIP(m2)
}

n <- nrow(cs); reg_f <- cs$region; ureg <- unique(reg_f)

draw_iid     <- function() sample.int(n, n, replace = TRUE)
draw_cluster <- function() {                       # resample regions
  rs <- sample(ureg, length(ureg), replace = TRUE)
  unlist(lapply(rs, function(r) which(reg_f == r)))
}
draw_two_stage <- function() {                     # regions, then countries
  rs <- sample(ureg, length(ureg), replace = TRUE)
  unlist(lapply(rs, function(r) {
    ii <- which(reg_f == r); sample(ii, length(ii), replace = TRUE)
  }))
}

run_boot <- function(drawer, label) {
  out <- vector("list", B)
  for (b in seq_len(B)) out[[b]] <- vip_diff_once(drawer())
  ok <- !vapply(out, is.null, logical(1))
  M <- do.call(rbind, out[ok])
  colnames(M) <- X_vars
  res <- data.frame(
    scheme = label, variable = X_vars,
    delta_vip = colMeans(M),
    lo = apply(M, 2, quantile, 0.025),
    hi = apply(M, 2, quantile, 0.975),
    width = apply(M, 2, quantile, 0.975) - apply(M, 2, quantile, 0.025),
    B_ok = sum(ok), row.names = NULL)
  res$excludes_zero <- res$lo > 0 | res$hi < 0
  res
}

cat("\n\n=== Paired bootstrap on the VIP difference (level minus shape) ===\n")
cat(sprintf("B = %d for each scheme\n", B))
bt <- bind_rows(
  run_boot(draw_iid,       "A countries iid"),
  run_boot(draw_cluster,   "B cluster on the 7 regions"),
  run_boot(draw_two_stage, "C two-stage: regions then countries")
)
print(bt %>% mutate(across(c(delta_vip, lo, hi, width), ~ round(.x, 3))))

cat("\n=== Interval widths by scheme (mean across regressors) ===\n")
print(bt %>% group_by(scheme) %>%
        summarise(mean_width = round(mean(width), 3),
                  n_excluding_zero = sum(excludes_zero), .groups = "drop"))


write.csv(loro, file.path(OUT, "23_loro_decomposed.csv"), row.names = FALSE)
write.csv(bt,   file.path(OUT, "23_vip_bootstrap.csv"),   row.names = FALSE)
cat("\nSaved: 23_loro_decomposed.csv, 23_vip_bootstrap.csv\n")

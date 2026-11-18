# 02: shape features of the life expectancy trajectory, and their PCA
#
# For each country, features of the 2005-2024 trajectory: trend slope,
# residual volatility, largest annual fall and rise, number of
# reversals, PELT changepoints in the annual increments, and the range
# and maximum of the segment slopes. The tsfeatures statistics are
# computed too but not used, because they have no direct reading in
# terms of health interventions or shocks. The PCA keeps the
# components needed for 80 per cent of the variance.
#
# Reads  data/df_le_long.rds
# Writes output/shape_features_post.rds, shape_pcs.rds,
#        shape_pca_loadings.{rds,csv}, shape_pca_variance.csv

source("00_common.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(changepoint); library(tsfeatures)
})

POST_Y0 <- 2005; POST_Y1 <- 2024

le_long <- readRDS(file.path(DATA, "df_le_long.rds"))

# 1. features on 2005-2024 ---------------------------------------

extract_features_post <- function(y, years) {
  ok <- !is.na(y); y <- y[ok]; yr <- years[ok]
  n <- length(y)
  if (n < 10) return(NULL)
  dy <- diff(y); n_dy <- length(dy)

  # PELT on the annual increments
  pelt_n <- NA_integer_; slope_range <- 0; slope_max <- NA_real_
  cp <- tryCatch(
    changepoint::cpt.meanvar(dy, method = "PELT",
                             penalty = "MBIC", minseglen = 3),
    error = function(e) NULL
  )
  if (!is.null(cp)) {
    cps <- changepoint::cpts(cp)
    pelt_n <- length(cps)
    seg_break <- c(0, cps, n_dy)
    seg_slopes <- sapply(seq_along(seg_break)[-1], function(i) {
      idx <- (seg_break[i-1]+1):seg_break[i]
      if (length(idx) >= 2) mean(dy[idx], na.rm = TRUE) else NA_real_
    })
    seg_slopes <- seg_slopes[is.finite(seg_slopes)]
    if (length(seg_slopes) >= 1) {
      slope_max   <- max(seg_slopes)
      slope_range <- diff(range(seg_slopes))
    }
  }

  # tsfeatures, computed for reference and dropped below
  ts_obj <- ts(y, start = yr[1])
  tsf <- tryCatch(
    tsfeatures(list(ts_obj),
               features = c("entropy", "stability",
                            "lumpiness", "spike",
                            "linearity", "curvature")),
    error = function(e) NULL
  )

  out <- data.frame(
    n_obs       = n,
    trend_slope = coef(lm(y ~ yr))[2],
    resid_sd    = sd(resid(lm(y ~ yr))),
    worst_drop  = min(dy),
    best_jump   = max(dy),
    n_reversal  = sum(dy < 0),
    pelt_n      = pelt_n,
    slope_range = slope_range,
    slope_max   = slope_max
  )
  if (!is.null(tsf)) out <- cbind(out, as.data.frame(tsf))
  out
}

shape_post <- le_long %>%
  filter(year >= POST_Y0, year <= POST_Y1) %>%
  arrange(iso3c, year) %>%
  group_by(iso3c) %>%
  do(extract_features_post(.$LE, .$year)) %>%
  ungroup() %>%
  filter(!is.na(trend_slope))

cat("Shape features, 2005-2024: ", nrow(shape_post), " countries\n")
cat("Features computed:\n"); print(names(shape_post))

# 2. features retained -------------------------------------------

selected_features <- c(
  "trend_slope",   # average annual gain
  "resid_sd",      # volatility around the trend
  "worst_drop",    # largest annual fall
  "best_jump",     # largest annual rise
  "n_reversal",    # number of annual falls
  "pelt_n",        # number of changepoints
  "slope_range",   # spread of the segment slopes
  "slope_max"      # steepest segment
)
dropped_features <- c("entropy","stability","lumpiness",
                      "spike","linearity","curvature")

cat("\nFeatures retained (", length(selected_features), "):\n  ",
    paste(selected_features, collapse=", "), "\n")
cat("Features dropped (", length(dropped_features), "):\n  ",
    paste(dropped_features, collapse=", "), "\n")

X_shape <- shape_post[, c("iso3c", selected_features)] %>% na.omit()

# 3. PCA ---------------------------------------------------------

mat <- as.matrix(X_shape[, -1])
pca <- prcomp(mat, center = TRUE, scale. = TRUE)

eig <- pca$sdev^2
var_expl <- eig / sum(eig)
cum_var  <- cumsum(var_expl)

cat("\nPCA, share of variance by component:\n")
for (i in seq_along(var_expl)) {
  cat(sprintf("  PC%d: %.3f   (cum: %.3f)\n", i, var_expl[i], cum_var[i]))
}

# components up to 80 per cent of the variance
n_pcs <- min(which(cum_var >= 0.80))
cat("\nComponents retained (>= 80% of variance):", n_pcs, "\n")

# 4. loadings and scores -----------------------------------------
loadings_PCs <- pca$rotation[, 1:n_pcs, drop = FALSE]
cat("\nLoadings of the retained components:\n")
print(round(loadings_PCs, 3))

scores <- as.data.frame(pca$x[, 1:n_pcs, drop = FALSE])
names(scores) <- paste0("PC", 1:n_pcs)
shape_pcs <- cbind(iso3c = X_shape$iso3c, scores)

# 5. save --------------------------------------------------------
saveRDS(shape_post,    file.path(OUT, "shape_features_post.rds"))
saveRDS(shape_pcs,     file.path(OUT, "shape_pcs.rds"))
saveRDS(loadings_PCs,  file.path(OUT, "shape_pca_loadings.rds"))
write.csv(round(loadings_PCs, 3),
          file.path(OUT, "shape_pca_loadings.csv"), row.names = TRUE)
write.csv(data.frame(PC = paste0("PC", seq_along(var_expl)),
                     var_expl = round(var_expl, 3),
                     cum_var  = round(cum_var, 3)),
          file.path(OUT, "shape_pca_variance.csv"), row.names = FALSE)

cat("\nSaved: shape_features_post.rds, shape_pcs.rds,\n")
cat("shape_pca_loadings.{rds,csv}, shape_pca_variance.csv\n")

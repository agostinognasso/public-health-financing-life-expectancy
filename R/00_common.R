# 00: setup, sourced by the other scripts (run from R/)
#
# Paths, the estimation panel used from script 16 on, and the shared
# helpers:
#   build_panel()     four-part financing composition, zero
#                     replacement and ILR coordinates
#   ilr_from_sbp()    ILR coordinates from a sequential binary partition
#   longrun()         long-run multipliers with delta-method and
#                     Fieller intervals
#   fieller_levels()  confidence levels at which a Fieller set changes
#                     shape
#   save_fig()        writes a figure under its printed number

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

# plm is attached by the scripts that fit panel models, not here: it
# masks dplyr's lag, lead and between, which the others use.

DATA <- file.path("..", "data")
OUT  <- file.path("..", "output")
FIG  <- file.path("..", "figures")
LOG  <- file.path("..", "log")

for (d in c(OUT, file.path(OUT, "pls"), FIG, LOG))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Figures are saved under the number they carry in the article
# (Figure1, Figure2, FigureS1 to FigureS12), as vector PDF for the
# typesetter and as PNG at 300 dpi.
save_fig <- function(name, plot, width = 7, height = 5) {
  ggplot2::ggsave(file.path(FIG, paste0(name, ".pdf")), plot,
                  width = width, height = height,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(FIG, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = 300)
  cat("  wrote", name, "\n")
  invisible(name)
}

# regressors of the main specification
KEY_X <- c("log_GNI", "EDU_GROSS", "PUB", "OOP_share", "IMM_DPT", "SANITATION")

# ------------------------------------------------------------
# build_panel(): the analytical panel
#   year_max : 2024 (full) or 2019 (pre-pandemic window)
#   delta    : detection limit for zero replacement, in percentage
#              points of current health expenditure
#   zero_rule: "multiplicative" (Martin-Fernandez et al. 2003),
#              "additive" (parts below delta set to delta, then
#              reclosed), or "none" (zeros kept, so the ILR
#              coordinates are NA where a part is zero)
# ------------------------------------------------------------
build_panel <- function(year_min = 2000, year_max = 2024,
                        delta = 0.05, zero_rule = "multiplicative") {

  panel   <- readRDS(file.path(OUT, "panel_with_ext.rds"))
  le_long <- readRDS(file.path(DATA, "df_le_long.rds"))

  dat <- panel %>%
    inner_join(le_long %>% dplyr::select(iso3c, year, LE), by = c("iso3c", "year")) %>%
    filter(year >= year_min, year <= year_max) %>%
    filter(!is.na(LE), !is.na(log_GNI), !is.na(EDU_GROSS),
           !is.na(PUB), !is.na(IMM_DPT), !is.na(SANITATION),
           !is.na(OOP_share)) %>%
    arrange(iso3c, year)

  # ---- four-part composition, before replacement ------------
  dat <- dat %>%
    dplyr::mutate(c_PUB = PUB,
           c_OOP = OOP_share,
           c_PNO = pmax(PRIV - OOP_share, 0),
           c_EXT = pmax(100 - PUB - PRIV, 0))

  parts <- c("c_PUB", "c_OOP", "c_PNO", "c_EXT")
  dat$n_zero  <- rowSums(dat[, parts] <= 0)
  # A part recorded below delta is treated as below detection and
  # replaced like a zero. This includes the floating-point residues
  # of order 1e-15 that the subtraction for c_EXT leaves where
  # PUB + PRIV = 100, which would otherwise enter the log-ratios.
  dat$n_below <- rowSums(dat[, parts] < delta)

  # ---- zero replacement -------------------------------------
  cm <- as.matrix(dat[, parts])
  if (zero_rule == "additive") {
    cm <- pmax(cm, delta)
    cm <- cm / rowSums(cm) * 100
  } else if (zero_rule == "multiplicative") {
    # Martin-Fernandez, Barcelo-Vidal & Pawlowsky-Glahn (2003):
    # parts below delta -> delta; the others rescaled by
    # (1 - sum(delta)/total) so that their ratios are preserved.
    tot <- rowSums(cm)
    for (i in seq_len(nrow(cm))) {
      z <- cm[i, ] < delta
      if (any(z)) {
        k <- sum(z)
        cm[i, z]  <- delta
        cm[i, !z] <- cm[i, !z] * (1 - k * delta / tot[i])
      }
    }
    cm <- cm / rowSums(cm) * 100
  } else if (zero_rule == "none") {
    cm[cm <= 0] <- NA_real_
    cm <- cm / rowSums(cm) * 100
  } else stop("unknown zero_rule")

  colnames(cm) <- c("z_PUB", "z_OOP", "z_PNO", "z_EXT")
  dat <- bind_cols(dat, as.data.frame(cm))

  # ---- ILR coordinates ---------------------------------------
  sbp <- matrix(c(
    +1, -1, -1, -1,   # ILR1: PUB vs {OOP,PNO,EXT}
     0, +1, +1, -1,   # ILR2: {OOP,PNO} vs EXT
     0, +1, -1,  0    # ILR3: OOP vs PNO
  ), nrow = 3, byrow = TRUE,
     dimnames = list(c("ILR1", "ILR2", "ILR3"),
                     c("z_PUB", "z_OOP", "z_PNO", "z_EXT")))

  ilr <- ilr_from_sbp(dat[, colnames(sbp)], sbp)
  dat <- bind_cols(dat, as.data.frame(ilr))

  # ---- public spending as a level rather than a share -------
  # HEALTH_PC is current health expenditure per capita (current US$);
  # the public component is its share times that total.
  dat <- dat %>%
    mutate(PUB_PC     = HEALTH_PC * PUB / 100,
           log_PUB_PC = ifelse(PUB_PC > 0, log(PUB_PC), NA_real_))

  dat
}

# ILR coordinates from a sequential binary partition (as in script 08)
ilr_from_sbp <- function(comp, sbp) {
  comp <- as.matrix(comp)
  out <- apply(sbp, 1, function(row) {
    pos <- which(row > 0); neg <- which(row < 0)
    r <- length(pos); s <- length(neg)
    sqrt(r * s / (r + s)) *
      (rowMeans(log(comp[, pos, drop = FALSE])) -
       rowMeans(log(comp[, neg, drop = FALSE])))
  })
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(comp))
  colnames(out) <- rownames(sbp)
  out
}

# ------------------------------------------------------------
# Long-run multiplier beta/(1-alpha) with a delta-method SE and a
# Fieller (1954) set, which keeps its coverage when 1-alpha is
# imprecisely estimated.
# ------------------------------------------------------------
longrun <- function(b, V, alpha_name, level = 0.95) {
  za <- qnorm(1 - (1 - level) / 2)
  a  <- b[alpha_name]
  va <- V[alpha_name, alpha_name]
  den <- 1 - a

  terms <- setdiff(names(b), alpha_name)
  res <- lapply(terms, function(tn) {
    bb  <- b[tn]; vb <- V[tn, tn]; cab <- V[tn, alpha_name]
    lr  <- bb / den
    # delta method: d/db = 1/(1-a); d/da = b/(1-a)^2
    g1 <- 1 / den; g2 <- bb / den^2
    v  <- g1^2 * vb + g2^2 * va + 2 * g1 * g2 * cab
    se <- sqrt(pmax(v, 0))

    # Fieller set for the ratio bb/den, where
    # den = 1 - a has variance va and Cov(bb, den) = -cab
    cbd <- -cab
    A <- den^2 - za^2 * va
    B <- -2 * (bb * den - za^2 * cbd)
    C <- bb^2 - za^2 * vb
    disc <- B^2 - 4 * A * C
    if (A > 0 && disc >= 0) {
      r <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
      fl <- r[1]; fu <- r[2]; ftype <- "bounded"
    } else if (disc >= 0) {
      # 1 - alpha not significantly different from zero at this
      # level: the set is the complement of an interval
      r <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
      fl <- r[1]; fu <- r[2]; ftype <- "exclusive/unbounded"
    } else {
      fl <- -Inf; fu <- Inf; ftype <- "whole line"
    }

    data.frame(term = tn, beta_SR = bb, se_SR = sqrt(vb),
               beta_LR = lr, se_LR = se,
               lo_delta = lr - za * se, hi_delta = lr + za * se,
               fieller_lo = fl, fieller_hi = fu, fieller_type = ftype,
               row.names = NULL)
  })
  out <- do.call(rbind, res)
  attr(out, "alpha")     <- a
  attr(out, "se_alpha")  <- sqrt(va)
  attr(out, "one_minus") <- den
  out
}

# ------------------------------------------------------------
# Confidence levels at which the Fieller set for b/(1-alpha)
# changes shape (Dufour 1997). Below level_bounded it is an
# interval; this level depends on alpha alone, being the level at
# which 1-alpha stops being significantly different from zero.
# Between level_bounded and level_wholeline it is the union of
# two half-lines; above level_wholeline it is the whole real line.
# ------------------------------------------------------------
fieller_levels <- function(b, V, alpha_name) {
  a   <- b[alpha_name]
  va  <- V[alpha_name, alpha_name]
  den <- 1 - a

  terms <- setdiff(names(b), alpha_name)
  res <- lapply(terms, function(tn) {
    bb  <- b[tn]; vb <- V[tn, tn]
    cbd <- -V[tn, alpha_name]                 # Cov(bb, 1 - a)
    # the discriminant of the Fieller quadratic is non-negative
    # for z^2 up to q_whole
    det     <- va * vb - cbd^2
    q_whole <- (den^2 * vb - 2 * bb * den * cbd + bb^2 * va) / det
    data.frame(term = tn,
               level_bounded   = 2 * pnorm(abs(den) / sqrt(va)) - 1,
               level_wholeline = if (det > 0) 2 * pnorm(sqrt(q_whole)) - 1
                                 else NA_real_,
               row.names = NULL)
  })
  do.call(rbind, res)
}

# a Fieller set written out according to its shape
fmt_fieller <- function(lo, hi, type, digits = 4) {
  f <- function(x) formatC(x, format = "f", digits = digits, flag = "+")
  switch(type,
         "bounded"             = sprintf("[%s, %s]", f(lo), f(hi)),
         "exclusive/unbounded" = sprintf("(-Inf, %s] U [%s, +Inf)", f(lo), f(hi)),
         "whole line"          = "the whole real line")
}

# robust vcov for a pgmm fit, matching summary(., robust = TRUE)
pgmm_vcov <- function(m) {
  V <- tryCatch(plm::vcovHC(m), error = function(e) NULL)
  if (is.null(V)) V <- vcov(m)
  V
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

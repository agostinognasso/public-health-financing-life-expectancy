# 26: Kiviet-Bruno bias-corrected LSDV (Table 2, rung 4)
#
# The analytical bias correction of Kiviet (1995), extended to
# unbalanced panels by Bruno (2005). Stata computes it with xtlsdvc;
# we found no CRAN implementation, so the approximation formulas are
# coded here and checked by simulation before they are applied.
#
#   y_it = gamma y_i,t-1 + x_it'beta + eta_i + eps_it
#
#   c1(Tbar^-1)      = s2 tr(Pi) q1
#   c2(N^-1 Tbar^-1) = -s2 { Q Wb'Pi Ms Wb + tr(Q Wb'Pi Ms Wb) I_k
#                            + 2 s2 q11 tr(Pi'PiPi) I_k } q1
#   c3(N^-1 Tbar^-2) = s4 tr(Pi) [ 2 q11 Q Wb'Pi Pi'Wb q1
#                        + { q1'Wb'Pi Pi'Wb q1 + q11 tr(Q Wb'Pi Pi'Wb)
#                            + 2 tr(Pi'Pi Pi'Pi) q11^2 } q1 ]
#   Q = { Wb'Ms Wb + s2 tr(Pi'Pi) e1 e1' }^-1,  q1 = Q e1, q11 = e1'q1
#   Ms = S{I - D(D'SD)^-1 D'}S,  Pi = Ms L Gamma,
#   L = I_N kron L_T,  Gamma = I_N kron (I_T - gamma L_T)^-1
#
#   B1 = c1 ; B2 = B1 + c2 ; B3 = B2 + c3 ; LSDVC_i = LSDV - B_i
#
# Standard errors follow the parametric bootstrap of Kiviet and
# Bun (2001): draw eps* ~ N(0, s2hat), regenerate y* from the
# recursion, and re-run the whole procedure, initial estimator
# included, on each bootstrap sample.
#
# All NT x NT objects are block diagonal by country, so every trace
# and cross-product is accumulated over N blocks of T x T.
#
# Reads  data/ through build_panel()
# Writes output/26_lsdvc.csv, 26_lsdvc_bootstrap.csv
# The number of bootstrap replications is 500 unless the environment
# variable LSDVC_BREP sets another value.

source("00_common.R")
suppressPackageStartupMessages(library(plm))
set.seed(20260921)

YEARS <- 2000:2024
TT    <- length(YEARS)
XV    <- c("log_GNI", "EDU_GROSS", "PUB", "OOP_share", "IMM_DPT", "SANITATION")
FOC   <- "PUB"

# ---- 1. The Bruno (2005) machinery -----------------------

make_LT <- function(T) { M <- matrix(0,T,T); if (T>1) M[cbind(2:T,1:(T-1))] <- 1; M }

make_GammaT <- function(g, T) {
  G <- matrix(0,T,T); for (t in 1:T) for (s in 1:t) G[t,s] <- g^(t-s); G
}

make_Msi <- function(s_i) {
  Ti <- sum(s_i)
  if (Ti == 0) return(matrix(0, length(s_i), length(s_i)))
  M <- -outer(s_i, s_i)/Ti; diag(M) <- diag(M) + s_i; M
}

lsdvc_bias <- function(S, Wbar, gamma, sigma2, gamma_cap = 0.995) {
  if (!is.finite(gamma)) stop("non-finite gamma")
  capped <- abs(gamma) >= gamma_cap
  if (capped) gamma <- sign(gamma)*gamma_cap
  N <- nrow(S); T <- ncol(S); k <- dim(Wbar)[3]
  LG <- make_LT(T) %*% make_GammaT(gamma, T)

  tr_Pi <- 0; tr_PtP <- 0; tr_PtPP <- 0; tr_PtPPtP <- 0
  WMW <- matrix(0,k,k); WPMW <- matrix(0,k,k); WPPW <- matrix(0,k,k)
  for (i in 1:N) {
    s_i <- S[i,]; if (sum(s_i) == 0) next
    Ms <- make_Msi(s_i); Pi <- Ms %*% LG; PtP <- crossprod(Pi)
    tr_Pi     <- tr_Pi     + sum(diag(Pi))
    tr_PtP    <- tr_PtP    + sum(diag(PtP))
    tr_PtPP   <- tr_PtPP   + sum(diag(PtP %*% Pi))
    tr_PtPPtP <- tr_PtPPtP + sum(diag(PtP %*% PtP))
    MW <- Ms %*% matrix(Wbar[i,,], nrow=T, ncol=k)
    WMW  <- WMW  + crossprod(MW)
    WPMW <- WPMW + crossprod(MW, LG %*% MW)
    WPPW <- WPPW + crossprod(LG %*% MW)
  }
  e1 <- c(1, rep(0,k-1))
  Qi <- WMW + sigma2*tr_PtP*tcrossprod(e1)
  Q  <- tryCatch(solve(Qi), error=function(e) MASS::ginv(Qi))
  q1 <- as.vector(Q %*% e1); q11 <- q1[1]; Ik <- diag(k)

  c1 <- sigma2*tr_Pi*q1
  A2 <- Q %*% WPMW
  c2 <- -sigma2*as.vector((A2 + sum(diag(A2))*Ik +
                           2*sigma2*q11*tr_PtPP*Ik) %*% q1)
  A3 <- Q %*% WPPW
  sc <- as.numeric(t(q1) %*% WPPW %*% q1) + q11*sum(diag(A3)) +
        2*tr_PtPPtP*q11^2
  c3 <- sigma2^2*tr_Pi*as.vector(2*q11*(A3 %*% q1) + sc*q1)

  list(B1=c1, B2=c1+c2, B3=c1+c2+c3, c1=c1, c2=c2, c3=c3, capped=capped)
}

# E(W) conditional on the initial value and the exogenous block.
#   prev tracks E(y_{i,t-1}); the recursion starts at the country's
#   first observed outcome, which is taken as given (Kiviet's y_0).
#   anchor = "first"  : one recursion per country, gaps propagated
#   anchor = "regap"  : the recursion restarts at every observed
#                       outcome, so no gap is propagated
make_Wbar <- function(X, y0, eta, gamma, beta, first_t,
                      anchor = "first", Yobs = NULL, obs = NULL) {
  N <- dim(X)[1]; T <- dim(X)[2]; kx <- dim(X)[3]
  Wb <- array(0, c(N, T, kx + 1))
  for (i in 1:N) {
    prev <- y0[i]                       # E(y at period first_t[i])
    ft <- first_t[i]
    if (ft < T) for (t in (ft + 1):T) {
      if (anchor == "regap" && !is.null(obs) && obs[i, t - 1])
        prev <- Yobs[i, t - 1]
      Wb[i, t, 1] <- prev
      prev <- gamma * prev + sum(X[i, t, ] * beta) + eta[i]
    }
    if (ft >= 1) Wb[i, 1:ft, 1] <- y0[i]
    Wb[i, , 2:(kx + 1)] <- X[i, , ]
  }
  Wb
}

lsdv_fit <- function(Y, W, S) {
  N <- nrow(S); T <- ncol(S); k <- dim(W)[3]
  XtX <- matrix(0,k,k); Xty <- numeric(k)
  for (i in 1:N) {
    idx <- which(S[i,]==1); if (length(idx) <= 1) next
    Wi <- matrix(W[i,idx,], nrow=length(idx), ncol=k); yi <- Y[i,idx]
    Wi <- sweep(Wi,2,colMeans(Wi)); yi <- yi - mean(yi)
    XtX <- XtX + crossprod(Wi); Xty <- Xty + crossprod(Wi,yi)
  }
  as.vector(qr.solve(XtX, Xty))
}

# ---- 2. Validation by simulation -------------------------
cat("=== Validation by simulation: bias of LSDV and of the corrections ===\n")

sim_panel <- function(N, T, g, b, s2=1, miss=0, burn=100) {
  eta <- rnorm(N)
  X <- array(rnorm(N*(T+1)), c(N,T+1,1))
  yl <- rnorm(N, sd=sqrt(s2/(1-g^2))) + eta/(1-g)
  for (u in 1:burn) yl <- g*yl + b*rnorm(N) + eta + rnorm(N, sd=sqrt(s2))
  y0 <- g*yl + b*X[,1,1] + eta + rnorm(N, sd=sqrt(s2))
  Y <- matrix(0,N,T); prev <- y0
  for (t in 1:T) { prev <- g*prev + b*X[,t+1,1] + eta + rnorm(N,sd=sqrt(s2))
                   Y[,t] <- prev }
  R <- matrix(TRUE,N,T)
  if (miss>0) R[sample(length(R), floor(miss*length(R)))] <- FALSE
  S <- matrix(0,N,T)
  for (i in 1:N) for (t in 1:T)
    S[i,t] <- as.integer(R[i,t] && (t==1 || R[i,t-1]))
  W <- array(0,c(N,T,2)); W[,1,1] <- y0
  if (T>1) W[,2:T,1] <- Y[,1:(T-1)]
  W[,,2] <- array(X[,2:(T+1),,drop=FALSE],c(N,T,1))
  list(Y=Y, X=array(X[,2:(T+1),,drop=FALSE],c(N,T,1)), W=W,
       y0=y0, eta=eta, S=S)
}

g0 <- 0.87; b0 <- 0.5; Nv <- 138; Tv <- 25; Rv <- 200
for (mv in c(0, 0.32)) {
  out <- matrix(NA, Rv, 4, dimnames=list(NULL,c("lsdv","B1","B2","B3")))
  for (r in 1:Rv) {
    d <- sim_panel(Nv,Tv,g0,b0,miss=mv)
    dl <- lsdv_fit(d$Y, d$W, d$S)
    Wb <- make_Wbar(d$X, d$y0, d$eta, g0, b0, rep(0, Nv))
    bb <- lsdvc_bias(d$S, Wb, g0, 1)
    out[r,] <- c(dl[1], dl[1]-bb$B1[1], dl[1]-bb$B2[1], dl[1]-bb$B3[1])
  }
  cat(sprintf("\n-- N=%d T=%d gamma=%.2f, %.0f%% of cells missing, %d reps --\n",
              Nv,Tv,g0,100*mv,Rv))
  cat(sprintf("%-6s %10s %10s %10s\n","", "mean","bias","RMSE"))
  for (nm in colnames(out)) {
    v <- out[,nm]
    cat(sprintf("%-6s %10.5f %10.5f %10.5f\n", nm, mean(v), mean(v)-g0,
                sqrt(mean((v-g0)^2))))
  }
}
cat("\nThe bias approximation is evaluated at the true parameters here,\n")
cat("so this isolates the formulas from the initial estimator.\n")

# ---- 3. The estimation panel on a common 2000-2024 grid ---
cat("\n=== Application to the estimation panel ===\n")

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
iso <- sort(unique(dat$iso3c)); N <- length(iso)
dat$ii <- match(dat$iso3c, iso); dat$tt <- match(dat$year, YEARS)

Y   <- matrix(NA_real_, N, TT)
obs <- matrix(FALSE, N, TT)
Xr  <- array(NA_real_, c(N, TT, length(XV)))
Y[cbind(dat$ii, dat$tt)]  <- dat$LE
obs[cbind(dat$ii, dat$tt)] <- TRUE
for (j in seq_along(XV)) Xr[cbind(dat$ii, dat$tt, j)] <- dat[[XV[j]]]

# usable: outcome and drivers at t, outcome at t-1
S <- matrix(0L, N, TT)
for (i in 1:N) for (t in 2:TT) S[i,t] <- as.integer(obs[i,t] && obs[i,t-1])

Ti <- rowSums(S)
cat(sprintf("countries %d, grid T = %d, usable (t,t-1) pairs = %d\n",
            N, TT, sum(S)))
cat(sprintf("countries with Ti >= 2: %d; mean Ti = %.1f\n",
            sum(Ti>=2), mean(Ti[Ti>0])))

# X on the grid, country-mean imputation for the recursion only
Xf <- Xr
for (i in 1:N) for (j in seq_along(XV)) {
  v <- Xf[i,,j]
  if (all(is.na(v))) { Xf[i,,j] <- 0 } else {
    v[is.na(v)] <- mean(v, na.rm=TRUE); Xf[i,,j] <- v }
}

# Year dummies. Period 1 (2000) never enters the dynamic sample, since
# a usable observation needs the outcome at t-1, so the base category is
# 2001 and the dummies run 2002..2024 -- the same set pgmm reports.
YD <- YEARS[-(1:2)]
Xall <- array(0, c(N, TT, length(XV)+length(YD)))
Xall[,,seq_along(XV)] <- Xf
for (j in seq_along(YD)) Xall[, match(YD[j], YEARS), length(XV)+j] <- 1
xnames <- c(XV, as.character(YD))
kx <- dim(Xall)[3]; k <- kx + 1

# W = (y_{-1}, X) using the observed outcome
W <- array(0, c(N, TT, k))
W[,2:TT,1] <- Y[,1:(TT-1)]
W[,,2:k] <- Xall
W[is.na(W)] <- 0            # only at cells with S = 0

first_t <- apply(obs, 1, function(z) if (any(z)) which(z)[1] else TT)
y0 <- sapply(1:N, function(i) if (any(obs[i,])) Y[i, first_t[i]] else 0)

# ---- LSDV on the grid, checked against plm ------------------
dl <- lsdv_fit(Y, W, S); names(dl) <- c("alpha", xnames)
pd <- pdata.frame(dat, index = c("iso3c","year"))
rhs <- paste(XV, collapse = " + ")
m_tw <- plm(as.formula(paste("LE ~ lag(LE,1) +", rhs)), data = pd,
            model = "within", effect = "twoways")
cat(sprintf("\nLSDV alpha: this code %.5f | plm twoways %.5f | diff %.2e\n",
            dl["alpha"], coef(m_tw)["lag(LE, 1)"],
            dl["alpha"] - coef(m_tw)["lag(LE, 1)"]))
cat(sprintf("LSDV beta(PUB): this code %+.5f | plm %+.5f\n",
            dl[FOC], coef(m_tw)[FOC]))

# ---- initial N-consistent estimators ------------------------
fg <- LE ~ lag(LE,1) + log_GNI + EDU_GROSS + PUB + OOP_share +
           IMM_DPT + SANITATION | lag(LE, 2:5)
m_ab <- pgmm(fg, data=pd, effect="twoways", model="twosteps",
             transformation="d",  collapse=TRUE)
m_bb <- pgmm(fg, data=pd, effect="twoways", model="twosteps",
             transformation="ld", collapse=TRUE)

# map a pgmm coefficient vector onto (alpha, xnames)
map_coef <- function(m) {
  cf <- coef(m); out <- setNames(numeric(k), c("alpha", xnames))
  out["alpha"] <- cf["lag(LE, 1)"]
  for (nm in xnames) if (nm %in% names(cf)) out[nm] <- cf[nm]
  out
}

# eta_i and sigma^2 given a coefficient vector
eta_sig <- function(Y, W, S, d) {
  eta <- numeric(N); num <- 0; n <- 0; Nu <- 0
  for (i in 1:N) {
    idx <- which(S[i,]==1); if (!length(idx)) next
    r <- Y[i,idx] - as.vector(matrix(W[i,idx,],nrow=length(idx),ncol=k) %*% d)
    eta[i] <- mean(r); Nu <- Nu + 1; n <- n + length(idx)
    num <- num + sum((r - eta[i])^2)
  }
  list(eta=eta, s2=num/(n - Nu - k), n=n, Nu=Nu)
}

# ---- the correction, for each initializer -------------------
run_lsdvc <- function(d_init, Ymat, Wmat, y0v, anchor="first") {
  es <- eta_sig(Ymat, Wmat, S, d_init)
  Wb <- make_Wbar(Xall, y0v, es$eta, d_init["alpha"], d_init[-1],
                  first_t, anchor=anchor, Yobs=Ymat, obs=obs)
  bb <- lsdvc_bias(S, Wb, d_init["alpha"], es$s2)
  dlv <- lsdv_fit(Ymat, Wmat, S)
  list(lsdv=dlv, s2=es$s2, bias=bb, capped=bb$capped,
       B1=dlv-bb$B1, B2=dlv-bb$B2, B3=dlv-bb$B3)
}

inits <- list(AB = map_coef(m_ab), BB = map_coef(m_bb))
cat("\n=== Initial N-consistent estimators ===\n")
for (nm in names(inits))
  cat(sprintf("%-3s: alpha = %.4f, beta(PUB) = %+.5f\n",
              nm, inits[[nm]]["alpha"], inits[[nm]][FOC]))

res <- lapply(inits, function(di) run_lsdvc(di, Y, W, y0))

cat("\n=== Bias-corrected LSDV (Kiviet 1995 / Bruno 2005) ===\n")
cat(sprintf("%-22s %10s %12s\n", "", "alpha", "beta(PUB)"))
cat(sprintf("%-22s %10.4f %12.5f\n", "LSDV (uncorrected)",
            res$BB$lsdv[1], res$BB$lsdv[which(xnames==FOC)+1]))
for (nm in names(res)) for (bo in c("B1","B2","B3")) {
  v <- res[[nm]][[bo]]
  cat(sprintf("%-22s %10.4f %12.5f\n",
              sprintf("LSDVC-%s init %s", bo, nm), v[1], v[which(xnames==FOC)+1]))
}
cat(sprintf("\nsigma^2: AB %.4f, BB %.4f\n", res$AB$s2, res$BB$s2))
cat(sprintf("bias term on alpha (BB): c1 %+.5f, c2 %+.5f, c3 %+.5f\n",
            res$BB$bias$c1[1], res$BB$bias$c2[1], res$BB$bias$c3[1]))

# ---- sensitivity to how E(y_{-1}) is anchored ---------------
cat("\n=== Sensitivity: construction of E(W) ===\n")
res_re <- lapply(inits, function(di) run_lsdvc(di, Y, W, y0, anchor="regap"))
cat(sprintf("%-6s %11s %11s %13s %13s\n","init","a (first)","a (regap)",
            "PUB (first)","PUB (regap)"))
for (nm in names(res)) {
  a1 <- res[[nm]]$B1; a2 <- res_re[[nm]]$B1
  cat(sprintf("%-6s %11.4f %11.4f %13.5f %13.5f\n", nm, a1[1], a2[1],
              a1[which(xnames==FOC)+1], a2[which(xnames==FOC)+1]))
}
cat("'first' propagates the recursion through gaps; 'regap' restarts it\n")
cat("at every observed outcome. The correction is read off 'first'.\n")

# ---- 4. Parametric bootstrap of Bun and Kiviet (2001) ----
BREP <- as.integer(Sys.getenv("LSDVC_BREP", "500"))
cat(sprintf("\n=== Parametric bootstrap, %d replications ===\n", BREP))

gen_star <- function(d_hat, eta, s2) {
  Ys <- matrix(NA_real_, N, TT)
  for (i in 1:N) {
    if (!any(obs[i,])) next
    prev <- y0[i]; Ys[i, first_t[i]] <- prev
    if (first_t[i] < TT) for (t in (first_t[i]+1):TT) {
      base <- if (obs[i,t-1] && !is.na(Ys[i,t-1])) Ys[i,t-1] else prev
      v <- d_hat["alpha"]*base + sum(Xall[i,t,]*d_hat[-1]) + eta[i] +
           rnorm(1, sd=sqrt(s2))
      if (obs[i,t]) Ys[i,t] <- v
      prev <- v
    }
  }
  Ys
}

boot_one <- function(d_hat, eta, s2) {
  Ys <- gen_star(d_hat, eta, s2)
  Ws <- W; Ws[,2:TT,1] <- Ys[,1:(TT-1)]; Ws[is.na(Ws)] <- 0
  db <- dat; db$LE <- Ys[cbind(dat$ii, dat$tt)]
  if (any(is.na(db$LE))) return(NULL)
  pdb <- pdata.frame(db, index=c("iso3c","year"))
  mb <- tryCatch(pgmm(fg, data=pdb, effect="twoways", model="twosteps",
                      transformation="ld", collapse=TRUE),
                 error=function(e) NULL)
  if (is.null(mb)) return(NULL)
  di <- map_coef(mb)
  if (!is.finite(di["alpha"])) return(NULL)
  r <- tryCatch(run_lsdvc(di, Ys, Ws, y0), error=function(e) NULL)
  if (is.null(r)) return(NULL)
  c(r$B1[1], r$B1[which(xnames==FOC)+1])
}

d_star <- res$BB$B1; names(d_star) <- c("alpha", xnames)
es_star <- eta_sig(Y, W, S, d_star)
bm <- matrix(NA_real_, BREP, 2)
t0 <- Sys.time()
for (b in 1:BREP) {
  v <- boot_one(d_star, es_star$eta, res$BB$s2)
  if (!is.null(v)) bm[b,] <- v
}
bm <- bm[complete.cases(bm),,drop=FALSE]
cat(sprintf("%d/%d replications usable, %.0f s\n", nrow(bm), BREP,
            as.numeric(Sys.time()-t0, units="secs")))

se_a <- sd(bm[,1]); se_b <- sd(bm[,2])
a_hat <- unname(res$BB$B1[1]); b_hat <- unname(res$BB$B1[which(xnames==FOC)+1])
cat(sprintf("\nalpha      = %.4f  (bootstrap SE %.4f)  95%% CI [%.4f, %.4f]\n",
            a_hat, se_a, a_hat-1.96*se_a, a_hat+1.96*se_a))
cat(sprintf("beta(PUB)  = %+.5f (bootstrap SE %.5f)  95%% CI [%+.5f, %+.5f]  p = %.3f\n",
            b_hat, se_b, b_hat-1.96*se_b, b_hat+1.96*se_b,
            2*pnorm(-abs(b_hat/se_b))))
cat(sprintf("1 - alpha  = %.4f  (SE %.4f)  95%% CI [%.4f, %.4f]\n",
            1-a_hat, se_a, 1-a_hat-1.96*se_a, 1-a_hat+1.96*se_a))
cat(sprintf("long-run multiplier beta/(1-alpha) = %+.4f\n", b_hat/(1-a_hat)))
cat(sprintf("MDE on the long-run multiplier at 80%% power = %.4f\n",
            2.802 * (se_b/(1-a_hat))))   # 2.802 = qnorm(0.975) + qnorm(0.80)

# Fieller set for beta/(1-alpha), from the bootstrap covariance
Vb <- cov(bm); colnames(Vb) <- rownames(Vb) <- c("alpha","beta")
cab <- Vb["beta","alpha"]; cbd <- -cab
den <- 1 - a_hat; za <- 1.96
A <- den^2 - za^2*Vb["alpha","alpha"]
B <- -2*(b_hat*den - za^2*cbd)
Cc <- b_hat^2 - za^2*Vb["beta","beta"]
disc <- B^2 - 4*A*Cc
cat(sprintf("\nbootstrap corr(alpha, beta) = %+.3f\n",
            cab/sqrt(Vb[1,1]*Vb[2,2])))
if (A > 0 && disc >= 0) {
  r <- sort((-B + c(-1,1)*sqrt(disc))/(2*A))
  cat(sprintf("Fieller 95%% set: bounded, [%+.4f, %+.4f]\n", r[1], r[2]))
} else if (disc >= 0) {
  r <- sort((-B + c(-1,1)*sqrt(disc))/(2*A))
  cat(sprintf("Fieller 95%% set: unbounded, the complement of [%+.4f, %+.4f]\n",
              r[1], r[2]))
} else cat("Fieller 95% set: the whole real line\n")
cat(sprintf("t on (1 - alpha) = %.2f\n", den/se_a))

# ---- the bias-correction rung, three ways -------------------
cat("\n=== The bias-correction rung, three ways ===\n")
yrs <- sort(unique(dat$year)); cut <- yrs[floor(length(yrs)/2)]
fit_half <- function(d) {
  p <- pdata.frame(d, index=c("iso3c","year"))
  tryCatch(coef(plm(as.formula(paste("LE ~ lag(LE,1) +", rhs)), data=p,
                    model="within", effect="twoways")), error=function(e) NULL)
}
h1 <- fit_half(dplyr::filter(dat, year <= cut))
h2 <- fit_half(dplyr::filter(dat, year >  cut))
cf <- coef(m_tw); nmv <- intersect(names(cf), intersect(names(h1), names(h2)))
hpj <- 2*cf[nmv] - (h1[nmv] + h2[nmv])/2
cat(sprintf("%-34s %9s %12s\n","method","alpha","beta(PUB)"))
cat(sprintf("%-34s %9.4f %12.5f\n","LSDV, uncorrected",
            cf["lag(LE, 1)"], cf[FOC]))
cat(sprintf("%-34s %9.4f %12.5f\n","half-panel jackknife (DJ 2015)",
            hpj["lag(LE, 1)"], hpj[FOC]))
cat(sprintf("%-34s %9.4f %12.5f\n","Kiviet-Bruno LSDVC (init BB)", a_hat, b_hat))

out <- data.frame(
  method = c("LSDV (uncorrected)","HPJ jackknife (Dhaene-Jochmans)",
             "LSDVC B1 (init AB)","LSDVC B2 (init AB)","LSDVC B3 (init AB)",
             "LSDVC B1 (init BB)","LSDVC B2 (init BB)","LSDVC B3 (init BB)"),
  alpha = c(cf["lag(LE, 1)"], hpj["lag(LE, 1)"],
            res$AB$B1[1], res$AB$B2[1], res$AB$B3[1],
            res$BB$B1[1], res$BB$B2[1], res$BB$B3[1]),
  beta_PUB = c(cf[FOC], hpj[FOC],
               res$AB$B1[which(xnames==FOC)+1], res$AB$B2[which(xnames==FOC)+1],
               res$AB$B3[which(xnames==FOC)+1],
               res$BB$B1[which(xnames==FOC)+1], res$BB$B2[which(xnames==FOC)+1],
               res$BB$B3[which(xnames==FOC)+1]),
  se_alpha = c(NA,NA,NA,NA,NA,NA,NA,se_a),
  se_beta  = c(NA,NA,NA,NA,NA,NA,NA,se_b),
  row.names = NULL)
write.csv(out, file.path(OUT, "26_lsdvc.csv"), row.names = FALSE)
write.csv(as.data.frame(bm), file.path(OUT, "26_lsdvc_bootstrap.csv"),
          row.names = FALSE)
cat("\nSaved: 26_lsdvc.csv, 26_lsdvc_bootstrap.csv\n")

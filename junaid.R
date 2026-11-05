## =====================================================================
## ATLAS-Ridge: computational timing harness  (self-contained)
## Implements the estimator and the ten competitors from the equations
## in Sections 3 and 4 of the manuscript. No external project code.
## =====================================================================
suppressPackageStartupMessages({ library(robustbase); library(MASS) })


}

## ------------------------------------------------- the ATLAS loss family
w_cq <- function(u, c, q) exp(-(abs(u) / c)^q)

hill_trimmed <- function(r) {                        # Equation (6)
  a <- sort(abs(r[r != 0]), decreasing = TRUE); ns <- length(a)
  if (ns < 10) return(4)
  m <- min(max(2L, floor(ns^0.6)), ns - 1L)
  1 / mean(log(a[1:m] / a[m + 1L]))
}
e_target <- function(ahat) 0.85 + 0.10 * min(1, ahat / 10)     # Eq (8)
q_shape  <- function(ahat) 2 + 2 * (1 - min(1, 3 / ahat))      # Eq (9)

c0_from_eff <- function(eff, q) {                              # Eq (7)
  f <- function(c) {
    num <- integrate(function(z) (w_cq(z, c, q) *
             (1 - q * (abs(z) / c)^q)) * dnorm(z), -8, 8)$value
    den <- integrate(function(z) z^2 * w_cq(z, c, q)^2 * dnorm(z),
             -8, 8)$value
    num^2 / den - eff
  }
  if (f(0.5) * f(12) > 0) return(2.5)
  uniroot(f, c(0.5, 12), tol = 1e-8)$root
}
calibrate_loss <- function(ahat) {
  q <- q_shape(ahat)
  list(ahat = ahat, q = q, c0 = c0_from_eff(e_target(ahat), q))
}
attenuation <- function(X, S, gamma = 1.5) {            # Section 3.4
  rd2 <- mahalanobis(X, S$center, S$cov)
  pmin(1, (qchisq(0.975, ncol(X)) / pmax(rd2, 1e-12))^(gamma / 2))
}

## ------------------------------------------- penalized M-fit at fixed k
irls_ridge <- function(X, y, ci, sig, q, k, b0, tol = 1e-8, maxit = 200) {
  b <- b0
  for (it in seq_len(maxit)) {
    r <- as.vector(y - X %*% b); wt <- w_cq(r / sig, ci, q)
    bn <- tryCatch(solve(crossprod(X, X * wt) + k * diag(ncol(X)),
                         crossprod(X, wt * y)), error = function(e) b)
    if (sqrt(sum((bn - b)^2)) / max(1e-12, sqrt(sum(b^2))) < tol) {
      b <- bn; break }
    b <- bn
  }
  list(beta = as.vector(b), w = wt)
}

KGRID <- c(0, 10^seq(-4, 4, length.out = 40))          # |K| = 41

gcv_grid <- function(X, y, ci, sig, q, b0, wref) {     # Eq (16)
  n <- nrow(X)
  ev <- eigen(crossprod(X * sqrt(wref)), symmetric = TRUE,
              only.values = TRUE)$values
  best <- list(gcv = Inf, k = 0, beta = b0)
  for (k in KGRID) {
    fit <- irls_ridge(X, y, ci, sig, q, k, b0)
    r <- as.vector(y - X %*% fit$beta)
    den <- (1 - sum(ev / (ev + k)) / n)^2
    g <- if (den > 1e-8) (sum(wref * r^2) / n) / den else Inf
    if (is.finite(g) && g < best$gcv)
      best <- list(gcv = g, k = k, beta = fit$beta)
  }
  best
}

## ------------------------------------------------------- the estimators
fit_ols <- function(X, y) qr.solve(X, y)
fit_ridge_hk <- function(X, y) {
  b <- qr.solve(X, y); r <- y - X %*% b
  k <- ncol(X) * (sum(r^2) / (nrow(X) - ncol(X))) / sum(b^2)
  as.vector(solve(crossprod(X) + k * diag(ncol(X)), crossprod(X, y)))
}
fit_ridge_gcv <- function(X, y) {
  n <- nrow(X); A <- crossprod(X)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  best <- list(g = Inf, b = NULL)
  for (k in KGRID) {
    b <- solve(A + k * diag(ncol(X)), crossprod(X, y)); r <- y - X %*% b
    g <- (sum(r^2) / n) / (1 - sum(ev / (ev + k)) / n)^2
    if (is.finite(g) && g < best$g) best <- list(g = g, b = b)
  }
  as.vector(best$b)
}
fit_m_huber <- function(X, y) as.vector(coef(rlm(X, y, psi = psi.huber,    maxit = 100)))
fit_m_bisq  <- function(X, y) as.vector(coef(rlm(X, y, psi = psi.bisquare, maxit = 100)))
fit_mm <- function(X, y) as.vector(coef(lmrob.fit(X, y,
  control = lmrob.control(setting = "KS2014", max.it = 100))))

rrm <- function(X, y, b0, w) {
  s2 <- sum(w * (y - X %*% b0)^2) / max(1, sum(w) - ncol(X))
  k  <- ncol(X) * s2 / sum(b0^2)
  as.vector(solve(crossprod(X, X * w) + k * diag(ncol(X)),
                  crossprod(X, w * y)))
}
fit_rrm_huber <- function(X, y) { f <- rlm(X, y, psi = psi.huber,    maxit = 100); rrm(X, y, coef(f), f$w) }
fit_rrm_bisq  <- function(X, y) { f <- rlm(X, y, psi = psi.bisquare, maxit = 100); rrm(X, y, coef(f), f$w) }
fit_rrmm <- function(X, y) {
  f <- lmrob.fit(X, y, control = lmrob.control(setting = "KS2014", max.it = 100))
  rrm(X, y, f$coefficients, f$rweights)
}
fit_gm_ridge <- function(X, y) {
  S <- covMcd(X, alpha = 0.75, nsamp = "deterministic")
  rrm(X, y, qr.solve(X, y), attenuation(X, S))
}

## ---------------------------------------------------------- ATLAS-Ridge
fit_atlas <- function(X, y, parts = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  S  <- covMcd(X, alpha = 0.75, nsamp = "deterministic")
  t1 <- proc.time()[["elapsed"]]
  b0 <- lmrob.fit(X, y, control = lmrob.control(setting = "KS2014", max.it = 100))
  bmm <- b0$coefficients; sig <- b0$scale
  t2 <- proc.time()[["elapsed"]]
  ah <- hill_trimmed(as.vector(y - X %*% bmm) / sig)
  t3 <- proc.time()[["elapsed"]]
  th <- calibrate_loss(ah)
  t4 <- proc.time()[["elapsed"]]
  ci <- th$c0 * attenuation(X, S)
  wref <- w_cq(as.vector(y - X %*% bmm) / sig, ci, th$q)
  best <- gcv_grid(X, y, ci, sig, th$q, bmm, wref)
  t5 <- proc.time()[["elapsed"]]
  if (!parts) return(best$beta)
  list(beta = best$beta, k = best$k, ahat = ah, q = th$q, c0 = th$c0,
       stage = c(detmcd = t1 - t0, mm_init = t2 - t1, tailidx = t3 - t2,
                 calib = t4 - t3, gcv = t5 - t4))
}

}
bench2 <- function(n, p, reps = 10L) {
  out <- setNames(numeric(length(FITS)), names(FITS))
  for (nm in names(FITS)) {
    v <- numeric(reps)
    for (r in seq_len(reps)) {
      d <- gen_data(n, p, 0.95, "t3", 1, "bad10")
      v[r] <- tryCatch(time_one(FITS[[nm]], d$X, d$y),
                       error = function(e) NA_real_)
    }
    out[nm] <- median(v, na.rm = TRUE)
  }
  cat("  done n =", n, " p =", p, "\n"); flush.console()
  c(n = n, p = p, out)
}

set.seed(20261015)
cat("\n===== 1. timing (this is the slow part)\n")
timing <- as.data.frame(do.call(rbind, list(
  bench2(100, 8), bench2(200, 12), bench2(500, 20))))

cat("\n--- milliseconds per fit\n")
ms <- timing; ms[, -(1:2)] <- round(1000 * ms[, -(1:2)], 3)
print(format(ms, scientific = FALSE))
cat("\n--- ATLAS / MM (its own initializer)\n"); print(round(timing$ATLAS / timing$MM, 2))
cat("\n--- ATLAS / RRMM\n");                     print(round(timing$ATLAS / timing$RRMM, 2))
cat("\n--- ATLAS / GM_Ridge\n");                 print(round(timing$ATLAS / timing$GM_Ridge, 2))
cat("\n--- ATLAS / OLS\n");                      print(round(timing$ATLAS / timing$OLS, 0))

## =====================================================================
## 2. ATLAS-Ridge decomposed by pipeline stage
## =====================================================================
set.seed(20261015); acc <- NULL
for (r in seq_len(20L)) {
  d <- gen_data(200, 12, 0.95, "t3", 1, "bad10")
  acc <- rbind(acc, fit_atlas(d$X, d$y, parts = TRUE)$stage)
}
stage <- apply(acc, 2, median)
cat("\n===== 2. stage decomposition, n = 200, p = 12\n"); print(round(stage, 4))
cat("\npercent of total:\n"); print(round(100 * stage / sum(stage), 1))

## =====================================================================
## 3. agreement check against the manuscript
## =====================================================================
set.seed(20261015); res <- NULL
for (r in seq_len(30L)) {
  d <- gen_data(100, 8, 0.95, "t3", 1, "bad10")
  res <- rbind(res, sapply(FITS, function(f) tryCatch(
    sum((f(d$X, d$y) - d$beta)^2), error = function(e) NA_real_)))
}
cat("\n===== 3. median MSE, 30 reps, n=100 p=8 rho=0.95 t3 bad-lev 10%\n")
print(round(apply(res, 2, median, na.rm = TRUE), 4))

## =====================================================================
## 4. environment
## =====================================================================
cat("\n===== 4. environment\n"); print(sessionInfo())
cat("\ncores:", parallel::detectCores(), "\n")
if (.Platform$OS.type == "windows") {
  cat("CPU :", system("wmic cpu get name", intern = TRUE)[2], "\n")
  cat("RAM :", round(as.numeric(system(
    "wmic computersystem get TotalPhysicalMemory", intern = TRUE)[2])/2^30, 1), "GB\n")
}

# Penalised maximum-entropy estimator ------------------------------------------

# Damped Newton solver for the penalised maximum-entropy dual at one barrier d.
# k: strikes / Fp (first element 0 = forward constraint), c: prices / Fp,
# sig: price tolerances, DF: discount factor, support of x = [0, xmax] with a
# uniform prior; the mass below d is the PoD. Equity Y = DF * (x - d).
solve_me <- function(k, c, sig, DF, d, xmax, n, theta0 = NULL, tol = 1e-9, maxit = 100) {
  B <- length(k)
  x <- sort(unique(c(seq(d, xmax, length.out = n), (d + k)[d + k < xmax])))
  dx <- diff(x)
  h <- c(dx, 0) / 2 + c(0, dx) / 2
  G <- DF * pmax(outer(x - d, k, "-"), 0)
  s2 <- sig^2
  th <- if (is.null(theta0)) rep(0, B) else theta0
  evalL <- function(th) {
    e <- as.vector(G %*% th)
    M <- max(0, e)
    w <- h * exp(e - M)
    Zs <- d * exp(-M) + sum(w)
    list(L = M + log(Zs) - log(xmax) - sum(th * c) + 0.5 * sum(s2 * th^2), w = w, Zs = Zs, M = M)
  }
  cur <- evalL(th)
  it <- 0L; conv <- FALSE
  repeat {
    p <- cur$w / cur$Zs
    m <- as.vector(crossprod(G, p))
    g <- m - c + s2 * th
    if (max(abs(g)) < tol) { conv <- TRUE; break }
    if (it >= maxit) break
    H <- crossprod(G * p, G) - tcrossprod(m)
    diag(H) <- diag(H) + s2
    step <- tryCatch(-solve(H, g), error = function(e) -g / pmax(diag(H), 1e-12))
    dec <- sum(g * step)
    if (dec >= 0) { step <- -g; dec <- -sum(g^2) }
    t <- 1
    repeat {
      new <- evalL(th + t * step)
      if (is.finite(new$L) && new$L <= cur$L + 1e-4 * t * dec) break
      t <- t / 2
      if (t < 1e-12) break
    }
    if (t < 1e-12) { conv <- max(abs(g)) < 1e3 * tol; break }
    th <- th + t * step
    cur <- new
    it <- it + 1L
  }
  p <- cur$w / cur$Zs
  list(theta = th, L = cur$L, log_pod = log(d) - cur$M - log(cur$Zs),
       fitted = as.vector(crossprod(G, p)), conv = conv, iter = it, x = x, h = h, p = p, grad = max(abs(g)))
}

# Solver inputs from a chain: quotes with a positive bid, tolerances from the
# half-spread (missing -> median), floored; forward row first.
me_inputs <- function(chain, sigma_floor, sigma_kappa, sigma_stock) {
  if (!inherits(chain, "ipod_chain")) stop("`chain` must be an `ipod_chain`.", call. = FALSE)
  if (!isTRUE(chain$ok)) stop("chain is not eligible: ", chain$flag, call. = FALSE)
  q <- chain$quotes[chain$quotes$bid > 0, , drop = FALSE]
  hs <- (q$ca - q$cb) / 2
  sig <- ifelse(is.finite(hs) & hs > 0, hs, NA)
  sig[is.na(sig)] <- stats::median(sig, na.rm = TRUE)
  sig[is.na(sig)] <- sigma_floor
  sig <- pmax(sigma_kappa * sig, sigma_floor)
  list(k = c(0, q$k / chain$DF), c = c(1, q$cm), sig = c(sigma_stock, sig), n = nrow(q))
}

default_sigma_stock <- function(chain) if (identical(chain$forward, "spot")) 0.02 else 0.002

#' Fit the penalised maximum-entropy risk-neutral distribution
#'
#' Fits the option-implied distribution of Vilsmeier (2016) in a penalised,
#' scale-free form. Asset value \eqn{x} (relative to the prepaid forward) has a
#' uniform prior on \eqn{[0, x_{max}]}; equity is \eqn{Y = \mathrm{DF}\,(x-d)^+}
#' on the \eqn{S_T/F} scale, and the prior mass below the barrier \eqn{d} is
#' the probability of zero equity (the iPoD). The prior is exponentially tilted
#' so that the forward and the call prices (puts converted by parity) are
#' matched within tolerances \eqn{\sigma_i}, by minimising the penalised dual
#' \deqn{\log Z(\lambda) - \lambda' c + \tfrac12 \sum_i \sigma_i^2 \lambda_i^2}
#' with a damped Newton method.
#'
#' Option prices bound but do not identify the probability of zero equity: any
#' value in the set returned by `tail_bounds(chain, 0)` is consistent with the
#' quotes, and the fitted PoD depends on `d` and the prior. Tail probabilities
#' \eqn{Q(Y < a)} with `a` inside the quoted strike range (`a >= chain$kmin`)
#' are much less sensitive to these choices; see [tail_prob()].
#'
#' @param chain An eligible [ipod_chain()].
#' @param d Barrier relative to the prepaid forward.
#' @param xmax Upper end of the asset-value support (relative to the prepaid
#'   forward).
#' @param n_grid Number of grid points of the asset-value support.
#' @param sigma_floor Lower bound of the price tolerances (call units).
#' @param sigma_kappa Multiplier of the half-spread tolerances.
#' @param sigma_stock Tolerance of the forward constraint. Defaults to 0.002
#'   for a parity forward and 0.02 for a spot forward, where it absorbs unknown
#'   payouts.
#' @param theta0 Optional starting values of the multipliers.
#' @param tol,maxit Convergence tolerance (maximum absolute gradient) and
#'   maximum number of Newton steps.
#' @return An object of class `ipod_fit`: a list with `d`, `pod`, `log_pod`,
#'   the multipliers `theta`, `converged`, `iter`, `grad`, `cost` (minus the
#'   dual objective), the observed and fitted call-unit prices `observed`,
#'   `fitted`, their tolerances `sigma`, `stock_error`, the grid `x` with
#'   probabilities `p`, and chain information `T`, `DF`, `F`, `Fp`, `kmin`,
#'   `xmax`.
#' @references Vilsmeier, J. (2016). Updating the option implied probability
#'   of default methodology. *Journal of Computational Finance* 19(3), 1--27.
#'   \doi{10.21314/jcf.2015.206}
#'
#'   Capuano, C. (2008). The option-iPoD: The probability of default implied
#'   by option prices based on entropy. IMF Working Paper 08/194.
#'   \doi{10.2139/ssrn.1339927}
#' @seealso [ipod_profile()], [tail_prob()], [tail_bounds()].
#' @examples
#' ch <- ipod_chain(sim_quotes(sigma = 0.4))
#' fit <- ipod_fit(ch, d = 0.1, n_grid = 500)
#' fit
#' tail_prob(fit, c(0.5, 0.8))
#' @export
ipod_fit <- function(chain, d = 0.1, xmax = 20, n_grid = 1500, sigma_floor = 5e-4,
                     sigma_kappa = 1, sigma_stock = NULL, theta0 = NULL, tol = 1e-9, maxit = 100) {
  if (is.null(sigma_stock)) sigma_stock <- default_sigma_stock(chain)
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d <= 0 || d >= xmax)
    stop("`d` must be a single number in (0, xmax).", call. = FALSE)
  inp <- me_inputs(chain, sigma_floor, sigma_kappa, sigma_stock)
  s <- solve_me(inp$k, inp$c, inp$sig, chain$DF, d, xmax, n_grid, theta0, tol, maxit)
  structure(list(d = d, pod = exp(s$log_pod), log_pod = s$log_pod, theta = s$theta,
                 converged = s$conv, iter = s$iter, grad = s$grad, cost = -s$L,
                 observed = inp$c[-1], fitted = s$fitted[-1], sigma = inp$sig[-1],
                 stock_error = s$fitted[1] - 1, x = s$x, h = s$h, p = s$p,
                 T = chain$T, DF = chain$DF, F = chain$F, Fp = chain$Fp, kmin = chain$kmin, xmax = xmax),
            class = "ipod_fit")
}

#' @export
print.ipod_fit <- function(x, ...) {
  z <- (x$fitted - x$observed) / x$sigma
  cat(sprintf("<ipod_fit> d = %.4g, PoD = %.3g (converged: %s, %d iterations)\n", x$d, x$pod, x$converged, x$iter))
  cat(sprintf("  %d prices, %.0f%% within tolerance, forward error %.2g\n", length(z), 100 * mean(abs(z) <= 1), x$stock_error))
  cat(sprintf("  kmin = %.3f (tail probabilities below kmin are not identified)\n", x$kmin))
  invisible(x)
}

#' Risk-neutral tail probabilities and quantiles of a fitted distribution
#'
#' `tail_prob()` returns \eqn{Q(Y < a)} with \eqn{Y = S_T/F}, including the
#' mass at zero (the PoD). `rnd_quantile()` returns quantiles of \eqn{Y}.
#' `rnd_summary()` collects the tail probabilities and shape measures used in
#' the accompanying paper.
#'
#' Only tail probabilities with `a` inside the quoted strike range
#' (`a >= kmin`) are close to quote-determined; below `kmin` they reflect the
#' prior and the barrier. Use [tail_bounds()] to see the range consistent with
#' the quotes.
#'
#' @param fit An [ipod_fit()].
#' @param a Thresholds on the \eqn{S_T/F} scale.
#' @param prob Probabilities.
#' @return `tail_prob()`: numeric vector of probabilities. `rnd_quantile()`:
#'   numeric vector of quantiles. `rnd_summary()`: named numeric vector with
#'   `pod`, `p_eq80`, `p_eq70`, `p_eq50`, `p_eq25`, `q05`, `iqr`, `bowley`,
#'   `sd`, `kl` (Kullback-Leibler divergence from the prior) and `kmin`.
#' @examples
#' fit <- ipod_fit(ipod_chain(sim_quotes()), n_grid = 500)
#' tail_prob(fit, 0.8)
#' rnd_quantile(fit, c(0.05, 0.5))
#' rnd_summary(fit)
#' @export
tail_prob <- function(fit, a) {
  stopifnot(inherits(fit, "ipod_fit"))
  y <- fit$DF * (fit$x - fit$d)
  cp <- cumsum(fit$p)
  vapply(a, function(ai) {
    if (!is.finite(ai)) return(NA_real_)
    j <- sum(y < ai)
    fit$pod + if (j > 0) cp[j] else 0
  }, numeric(1))
}

#' @rdname tail_prob
#' @export
rnd_quantile <- function(fit, prob) {
  stopifnot(inherits(fit, "ipod_fit"))
  y <- fit$DF * (fit$x - fit$d)
  cp <- fit$pod + cumsum(fit$p)
  vapply(prob, function(a) {
    if (!is.finite(a)) return(NA_real_)
    if (fit$pod >= a) return(0)
    j <- which(cp >= a)[1]
    if (is.na(j)) max(y) else y[j]
  }, numeric(1))
}

#' @rdname tail_prob
#' @export
rnd_summary <- function(fit) {
  stopifnot(inherits(fit, "ipod_fit"))
  y <- fit$DF * (fit$x - fit$d)
  p <- fit$p; pod <- fit$pod
  pe <- tail_prob(fit, c(0.8, 0.7, 0.5, 0.25))
  q <- rnd_quantile(fit, c(0.05, 0.25, 0.5, 0.75))
  dens <- p / pmax(fit$h, 1e-300)
  kl <- sum(p[p > 0] * log(dens[p > 0] * fit$xmax)) + if (pod > 0) pod * log(pod * fit$xmax / fit$d) else 0
  m1 <- sum(p * y)
  c(pod = pod, p_eq80 = pe[1], p_eq70 = pe[2], p_eq50 = pe[3], p_eq25 = pe[4], q05 = q[1],
    iqr = q[4] - q[2], bowley = if (q[4] > q[2]) (q[4] + q[2] - 2 * q[3]) / (q[4] - q[2]) else NA_real_,
    sd = sqrt(max(sum(p * y^2) - m1^2, 0)), kl = kl, kmin = fit$kmin)
}

#' Profile of the estimator over barriers
#'
#' Fits [ipod_fit()] over a grid of relative barriers, warm-starting each fit
#' from the previous converged solution. The profile shows how strongly the
#' PoD depends on the unidentified barrier, and is the input of [pod_select()].
#'
#' @inheritParams ipod_fit
#' @param d_grid Increasing grid of barriers relative to the prepaid forward.
#' @param ... Further arguments passed to [ipod_fit()] (not `d`, `theta0`).
#' @return A data frame with one row per barrier: `d`, `log_pod`, `pod`,
#'   `cost`, `converged`, `iter`, `fit_in` (share of prices within tolerance),
#'   `fit_rms_z`, `stock_error` and the elements of [rnd_summary()]. Attributes
#'   `Fp` and `kmin` are kept for [pod_select()].
#' @examples
#' ch <- ipod_chain(sim_quotes())
#' pr <- ipod_profile(ch, d_grid = c(0.02, 0.05, 0.1, 0.2, 0.4), n_grid = 500)
#' pr[, c("d", "pod", "p_eq80")]
#' @export
ipod_profile <- function(chain, d_grid = exp(seq(log(0.002), log(1.5), length.out = 28)), ...) {
  d_grid <- sort(d_grid)
  th <- NULL
  rows <- lapply(d_grid, function(d) {
    f <- ipod_fit(chain, d = d, theta0 = th, ...)
    th <<- if (f$converged) f$theta else NULL
    z <- (f$fitted - f$observed) / f$sigma
    data.frame(d = d, log_pod = f$log_pod, pod = f$pod, cost = f$cost, converged = f$converged,
               iter = f$iter, fit_in = mean(abs(z) <= 1), fit_rms_z = sqrt(mean(z^2)),
               stock_error = f$stock_error, t(rnd_summary(f)[-1]))
  })
  out <- do.call(rbind, rows)
  attr(out, "Fp") <- chain$Fp
  attr(out, "kmin") <- chain$kmin
  out
}

#' Select a PoD from a barrier profile
#'
#' Applies a barrier-selection rule to a profile from [ipod_profile()]. These
#' rules pick one point of the set of PoDs consistent with the quotes by
#' assumption; none of them identifies the PoD (compare [tail_bounds()] at
#' `a = 0`).
#'
#' * `"mean"`, `"median"`: among barriers in `range`, the PoD closest to the
#'   mean (median) PoD.
#' * `"absolute"`: the averaging rule of Vilsmeier (2016): barriers
#'   \eqn{D = 1, \dots, \lfloor\min(20, S)\rfloor} in currency units, i.e.
#'   \eqn{d = D/S} with \eqn{S} = `scale`; log-PoD is interpolated in log
#'   \eqn{d} and the PoD closest to the mean is returned. It depends on the
#'   nominal share price.
#' * `"maxent"`: the barrier with the smallest cost (least informative fit).
#' * `"fixed"`: the PoD interpolated at barrier `d`.
#'
#' @param profile Output of [ipod_profile()].
#' @param rule Selection rule, see Details.
#' @param range Barrier range for `"mean"` and `"median"`.
#' @param d Barrier for `"fixed"`.
#' @param scale Currency scale for `"absolute"`; defaults to the prepaid
#'   forward stored with the profile.
#' @return A one-row data frame with `rule`, `d` and `pod`.
#' @examples
#' ch <- ipod_chain(sim_quotes())
#' pr <- ipod_profile(ch, d_grid = exp(seq(log(0.01), log(1), length.out = 8)), n_grid = 500)
#' rbind(pod_select(pr, "mean"), pod_select(pr, "absolute"), pod_select(pr, "fixed", d = 0.1))
#' @export
pod_select <- function(profile, rule = c("mean", "median", "absolute", "maxent", "fixed"),
                       range = c(0.05, 0.5), d = 0.1, scale = attr(profile, "Fp")) {
  rule <- match.arg(rule)
  z <- profile[profile$converged & is.finite(profile$log_pod), , drop = FALSE]
  z <- z[order(z$d), , drop = FALSE]
  na <- data.frame(rule = rule, d = NA_real_, pod = NA_real_)
  if (nrow(z) < 2) return(na)
  ilp <- function(dt) stats::approx(log(z$d), z$log_pod, xout = log(dt), rule = 2)$y
  closest <- function(dd, pods, target) { i <- which.min(abs(pods - target)); data.frame(rule = rule, d = dd[i], pod = pods[i]) }
  switch(rule,
    mean = , median = {
      r <- z$d >= range[1] & z$d <= range[2]
      if (!any(r)) return(na)
      pods <- z$pod[r]
      closest(z$d[r], pods, if (rule == "mean") mean(pods) else stats::median(pods))
    },
    absolute = {
      if (is.null(scale) || !is.finite(scale) || scale <= 0) stop("`scale` must be a positive price.", call. = FALSE)
      dd <- seq_len(max(1, floor(min(20, scale)))) / scale
      pods <- exp(ilp(dd))
      closest(dd, pods, mean(pods))
    },
    maxent = { i <- which.min(z$cost); data.frame(rule = rule, d = z$d[i], pod = z$pod[i]) },
    fixed = data.frame(rule = rule, d = d, pod = exp(ilp(d))))
}

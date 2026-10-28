# Benchmarks: Black implied volatility, lognormal tails, equity-based measures --

#' Black (1976) prices and implied volatility on the forward scale
#'
#' `black_call()` returns the normalised call price \eqn{E^Q[(Y-k)^+]} of a
#' lognormal \eqn{Y = S_T/F} with mean one and volatility `sigma`;
#' `implied_vol()` inverts it by bisection. Prices in currency units are
#' converted by \eqn{c = C/(\mathrm{DF}\,F)}, and put prices by put-call parity.
#'
#' @param k Strike divided by the forward.
#' @param sigma Annualised volatility.
#' @param T Time to expiry in years.
#' @param price Option prices in currency units.
#' @param strike Strikes in currency units.
#' @param forward Forward price \eqn{F}.
#' @param DF Discount factor to expiry.
#' @param type `"call"` or `"put"` (recycled).
#' @return `black_call()`: normalised call prices. `implied_vol()`: implied
#'   volatilities (`NA` outside the no-arbitrage range).
#' @examples
#' c <- black_call(k = 0.9, sigma = 0.35, T = 0.5)
#' implied_vol(c * 0.99 * 100, strike = 90, forward = 100, T = 0.5, DF = 0.99)
#' @export
black_call <- function(k, sigma, T) {
  sT <- sigma * sqrt(T)
  d1 <- (-log(k) + sT^2 / 2) / sT
  stats::pnorm(d1) - k * stats::pnorm(d1 - sT)
}

#' @rdname black_call
#' @export
implied_vol <- function(price, strike, forward, T, DF = 1, type = "call") {
  n <- max(length(price), length(strike), length(forward), length(T), length(DF), length(type))
  price <- rep_len(price, n); k <- rep_len(strike / forward, n); T <- rep_len(T, n)
  call <- tolower(substr(rep_len(as.character(type), n), 1, 1)) == "c"
  c <- price / (rep_len(DF, n) * rep_len(forward, n))
  c[!call] <- c[!call] + 1 - k[!call]
  iv_norm(c, k, T)
}

# vectorised bisection on the normalised call price
iv_norm <- function(c, k, T, lo = 1e-4, hi = 8, iter = 100) {
  ok <- is.finite(c) & is.finite(k) & is.finite(T) & k > 0 & T > 0 & c > pmax(1 - k, 0) + 1e-12 & c < 1
  a <- rep(lo, length(c)); b <- rep(hi, length(c))
  ok <- ok & black_call(k, hi, T) > c
  for (i in seq_len(iter)) {
    m <- (a + b) / 2
    up <- black_call(k, m, T) < c
    a <- ifelse(up, m, a); b <- ifelse(up, b, m)
  }
  v <- (a + b) / 2
  v[!(ok %in% TRUE)] <- NA_real_
  v
}

#' At-the-money implied volatility of a chain
#'
#' Black implied volatility at the forward (\eqn{k = 1}), interpolated
#' linearly in \eqn{k} between the implied volatilities of the two
#' out-of-the-money mid quotes with positive bids closest to the forward. In
#' the accompanying paper this is the benchmark against which tail and
#' distress measures are compared.
#'
#' @param chain An [ipod_chain()].
#' @return The at-the-money implied volatility (annualised), or `NA`.
#' @examples
#' atm_iv(ipod_chain(sim_quotes(sigma = 0.3)))
#' @export
atm_iv <- function(chain) {
  if (!inherits(chain, "ipod_chain")) stop("`chain` must be an `ipod_chain`.", call. = FALSE)
  q <- chain$quotes
  if (is.null(q)) return(NA_real_)
  q <- q[q$bid > 0, , drop = FALSE]
  lo <- utils::tail(which(q$k < 1), 1); hi <- utils::head(which(q$k >= 1), 1)
  if (!length(lo) || !length(hi)) return(NA_real_)
  s <- iv_norm(q$cm[c(lo, hi)], q$k[c(lo, hi)], chain$T)
  if (q$k[hi] == q$k[lo]) return(s[2])
  s[1] + (s[2] - s[1]) * (1 - q$k[lo]) / (q$k[hi] - q$k[lo])
}

#' Lognormal tail probability
#'
#' \eqn{Q(Y \le a)} for a lognormal \eqn{Y = S_T/F} with mean one and
#' volatility `sigma`: \eqn{\Phi((\log a + \sigma^2 T/2)/(\sigma\sqrt T))}.
#' With `sigma` equal to at-the-money implied volatility this is the
#' benchmark that richer tail measures have to beat.
#'
#' @param a Thresholds on the \eqn{S_T/F} scale.
#' @inheritParams black_call
#' @return Probabilities.
#' @examples
#' lognormal_tail(c(0.5, 0.8), sigma = 0.4, T = 0.5)
#' @export
lognormal_tail <- function(a, sigma, T) {
  sT <- sigma * sqrt(T)
  stats::pnorm((log(a) + sT^2 / 2) / sT)
}

#' Constant-maturity interpolation
#'
#' Linear interpolation in maturity between the two observations bracketing
#' `target`. If all maturities lie on one side, the closest one is used only if
#' it is within a factor 2/3 (below) or 3/2 (above) of the target.
#'
#' @param ttm Maturities (any unit, same as `target`).
#' @param value Values at `ttm`.
#' @param target Target maturity.
#' @return The interpolated value, or `NA`.
#' @examples
#' const_maturity(c(90, 270), c(0.40, 0.30), target = 182)
#' @export
const_maturity <- function(ttm, value, target = 182) {
  ok <- is.finite(value) & is.finite(ttm); t <- ttm[ok]; v <- value[ok]
  lo <- which(t <= target); hi <- which(t >= target)
  if (!length(hi)) return(if (length(lo) && max(t[lo]) >= target * 2 / 3) v[lo][which.max(t[lo])] else NA_real_)
  if (!length(lo)) return(if (min(t[hi]) <= target * 3 / 2) v[hi][which.min(t[hi])] else NA_real_)
  a <- lo[which.max(t[lo])]; b <- hi[which.min(t[hi])]
  if (t[b] == t[a]) return(v[a])
  v[a] + (v[b] - v[a]) * (target - t[a]) / (t[b] - t[a])
}

#' Drawdown and distance to an equity-loss threshold
#'
#' `drawdown()` returns \eqn{1 - P_t / \max(P_{t-w+1}, \dots, P_t)} over a
#' trailing window of `window` observations (use split-adjusted prices).
#' `dtd()` converts a drawdown into a distance, in volatility units, of the
#' current price from the threshold \eqn{\theta} times the trailing peak:
#' \deqn{\frac{\log((1-dd)/\theta) - \sigma^2 T/2}{\sigma\sqrt{T}}.}
#'
#' @param price Price series, ordered in time.
#' @param window Window length in observations.
#' @param dd Drawdown(s) in \eqn{[0, 1)}.
#' @param sigma Volatility (for example at-the-money implied volatility).
#' @param T Horizon in years.
#' @param theta Threshold as a fraction of the trailing peak.
#' @return Numeric vector.
#' @examples
#' p <- c(100, 105, 98, 90, 95, 80)
#' drawdown(p, window = 4)
#' dtd(0.3, sigma = 0.5)
#' @export
drawdown <- function(price, window = 252) {
  n <- length(price)
  vapply(seq_len(n), function(i) {
    if (!is.finite(price[i])) return(NA_real_)
    pk <- max(price[max(1, i - window + 1):i], na.rm = TRUE)
    1 - price[i] / pk
  }, numeric(1))
}

#' @rdname drawdown
#' @export
dtd <- function(dd, sigma, T = 0.5, theta = 0.3) {
  (log((1 - dd) / theta) - sigma^2 * T / 2) / (sigma * sqrt(T))
}

#' Anchored distress probability
#'
#' The anchored distress probability (aDP) is the risk-neutral probability
#' that the share price ends below a fixed fraction `theta` of its trailing
#' peak, \eqn{Q(S_T < \theta M)}, i.e. \eqn{Q(Y < a)} with
#' \eqn{a = \theta M / F}. It combines the bank's equity history (through the
#' anchor) with the option-implied distribution, is scale-free, and has no
#' fitted weights. It is a terminal event: the path between now and expiry
#' does not matter.
#'
#' @param peak Trailing peak of split-adjusted prices (for example 252 days).
#' @param forward Forward price \eqn{F}.
#' @param theta Fraction of the peak.
#' @param fit An [ipod_fit()]; if `NULL`, `sigma` and `T` give the lognormal
#'   version.
#' @param sigma,T Volatility and maturity for the lognormal version.
#' @return `adp_threshold()`: the threshold \eqn{a}. `adp()`: the probability.
#' @examples
#' ch <- ipod_chain(sim_quotes(sigma = 0.5))
#' fit <- ipod_fit(ch, n_grid = 500)
#' adp(peak = 1.6 * ch$F, forward = ch$F, fit = fit)
#' adp(peak = 1.6 * ch$F, forward = ch$F, sigma = atm_iv(ch), T = ch$T)
#' @export
adp_threshold <- function(peak, forward, theta = 0.5) theta * peak / forward

#' @rdname adp_threshold
#' @export
adp <- function(peak, forward, theta = 0.5, fit = NULL, sigma = NULL, T = NULL) {
  a <- adp_threshold(peak, forward, theta)
  if (!is.null(fit)) return(tail_prob(fit, a))
  if (is.null(sigma) || is.null(T)) stop("Provide `fit`, or `sigma` and `T`.", call. = FALSE)
  lognormal_tail(a, sigma, T)
}

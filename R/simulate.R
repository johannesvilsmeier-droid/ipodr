# Simulated option quotes with a known risk-neutral distribution ---------------

#' Simulate an option chain with known tail probabilities
#'
#' Generates call and put quotes for one expiry from a distribution with an
#' atom of mass `pod` at zero and a lognormal remainder:
#' \eqn{Y = S_T/F = 0} with probability `pod`, otherwise
#' \eqn{Y = Z/(1-\mathrm{pod})} with \eqn{Z} lognormal with mean one and
#' volatility `sigma`. The true tail probability is
#' \eqn{Q(Y \le a) = \mathrm{pod} + (1-\mathrm{pod})\,\Phi((\log(a(1-\mathrm{pod})) + \sigma^2T/2)/(\sigma\sqrt T))},
#' available as `attr(x, "tail")(a)`. Bid and ask quotes are the model price
#' minus and plus half a spread of `max(spread_abs, spread_rel * price)`,
#' rounded outward to `tick`.
#'
#' @param ttm Time to expiry in calendar days.
#' @param sigma Volatility of the lognormal part.
#' @param pod Probability of zero equity.
#' @param forward Forward price.
#' @param rate Continuously compounded annual rate.
#' @param strikes Strikes as fractions of the forward.
#' @param spread_abs,spread_rel Absolute and relative quoted spread.
#' @param tick Tick size.
#' @param oi Open interest assigned to every quote.
#' @return A data frame in the long format expected by [ipod_chain()], with
#'   columns `ttm`, `strike`, `type`, `bid`, `ask`, `oi`, `rate`, `spot`, and
#'   an attribute `tail`, a function returning the true \eqn{Q(Y \le a)}.
#' @examples
#' q <- sim_quotes(sigma = 0.4, pod = 0.02)
#' head(q)
#' attr(q, "tail")(c(0, 0.5, 0.8))
#' @export
sim_quotes <- function(ttm = 182, sigma = 0.35, pod = 0, forward = 100, rate = 0.03,
                       strikes = seq(0.3, 1.6, by = 0.05), spread_abs = 0.05, spread_rel = 0.05,
                       tick = 0.01, oi = 100) {
  stopifnot(pod >= 0, pod < 1, sigma > 0, ttm > 0)
  T <- ttm / 365; DF <- exp(-rate * T)
  K <- strikes * forward
  cn <- black_call(strikes * (1 - pod), sigma, T)
  C <- DF * forward * cn
  P <- C - DF * (forward - K)
  quote <- function(v) {
    hs <- pmax(spread_abs, spread_rel * v) / 2
    list(bid = pmax(floor((v - hs) / tick + 1e-9) * tick, 0), ask = ceiling((v + hs) / tick - 1e-9) * tick)
  }
  qc <- quote(C); qp <- quote(pmax(P, 0))
  out <- data.frame(ttm = ttm, strike = rep(K, 2), type = rep(c("call", "put"), each = length(K)),
                    bid = c(qc$bid, qp$bid), ask = c(qc$ask, qp$ask), oi = oi, rate = rate, spot = DF * forward)
  out <- out[out$ask > 0, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "tail") <- function(a) ifelse(a <= 0, pod,
    pod + (1 - pod) * stats::pnorm((log(pmax(a, 1e-300) * (1 - pod)) + sigma^2 * T / 2) / (sigma * sqrt(T))))
  out
}

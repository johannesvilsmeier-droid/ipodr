# Model-free bounds from bid and ask quotes (single maturity) -------------------
# Closed forms of Dillschneider (2026), Propositions 3.1 and 3.3, for the
# normalised call price c(z) = E[(Y - z)^+], E[Y] = 1, Y >= 0.
# kk: strikes incl. 0 (sorted, unique); BB, AA: call-unit bids / asks with
# BB[1] = AA[1] = 1 (the forward).

ub_call <- function(z, kk, AA) {
  i <- kk <= z; u <- min(1, AA[i])
  if (any(i) && any(!i)) {
    w <- outer(kk[i], kk[!i], function(a, b) (b - z) / (b - a))
    u <- min(u, w * AA[i] + (1 - w) * rep(AA[!i], each = sum(i)))
  }
  u
}

lb_call <- function(z, kk, BB, AA) {
  l <- max(1 - z, 0, BB[kk >= z])
  g <- kk >= z
  if (sum(g) >= 2) { k1 <- kk[g]; W <- outer(k1, k1, function(a, b) (b - z) / (b - a))
    V <- W * BB[g] + (1 - W) * rep(AA[g], each = sum(g)); l <- max(l, V[upper.tri(V)]) }
  s <- kk <= z
  if (sum(s) >= 2) { k2 <- kk[s]; W <- outer(k2, k2, function(a, b) (b - z) / (b - a))
    V <- W * AA[s] + (1 - W) * rep(BB[s], each = sum(s)); l <- max(l, V[upper.tri(V)]) }
  l
}

ltp_bounds <- function(a, kk, BB, AA, eta = 1e-9) {
  UP <- function(z) ub_call(z, kk, AA) + z - 1
  ku <- unique(c(kk[kk > a], if (1 > a) 1))
  la <- lb_call(a, kk, BB, AA) + a - 1
  U <- min(1, vapply(ku, function(k) (UP(k) - la) / (k - a), 0))
  kl <- unique(c(kk[kk <= a], if (1 <= a) 1, a))
  le <- lb_call(a + eta, kk, BB, AA) + a + eta - 1
  L <- max(0, vapply(kl, function(k) (le - UP(k)) / (a + eta - k), 0))
  c(L = L, U = U)
}

quotes_feasible <- function(kk, BB, AA) {
  zs <- sort(unique(c(kk, (kk[-1] + kk[-length(kk)]) / 2, 1)))
  all(vapply(zs, function(z) lb_call(z, kk, BB, AA) <= ub_call(z, kk, AA) + 1e-10, TRUE))
}

norm_quotes <- function(chain, w) {
  q <- chain$quotes
  B <- pmax(q$bid - w, 0) / chain$Fp; A <- (q$ask + w) / chain$Fp
  B <- ifelse(q$call, B, B + 1 - q$k); A <- ifelse(q$call, A, A + 1 - q$k)
  list(kk = c(0, q$k), BB = c(1, pmax(B, pmax(1 - q$k, 0))), AA = c(1, pmin(A, 1)))
}

#' Normalised quotes consistent with no arbitrage
#'
#' Converts a chain into forward-normalised call bids and asks and widens all
#' quotes by `0, tick, 2 * tick, ...` currency units until they admit an
#' arbitrage-free distribution. Used internally by the bound functions.
#'
#' @param chain An [ipod_chain()].
#' @param tick Tick size in currency units used for widening.
#' @param max_widen Maximum number of ticks.
#' @return A list with `k` (strikes on the \eqn{K/F} scale, including 0),
#'   `bid`, `ask` (call units, first element the forward) and `widen` (ticks
#'   used; `NA` if no widening up to `max_widen` is feasible).
#' @examples
#' arbitrage_free_quotes(ipod_chain(sim_quotes()))$widen
#' @export
arbitrage_free_quotes <- function(chain, tick = 0.01, max_widen = 2) {
  if (!inherits(chain, "ipod_chain")) stop("`chain` must be an `ipod_chain`.", call. = FALSE)
  if (is.null(chain$quotes) || !nrow(chain$quotes)) return(list(k = NA, bid = NA, ask = NA, widen = NA_integer_))
  for (j in 0:max_widen) {
    q <- norm_quotes(chain, j * tick)
    if (quotes_feasible(q$kk, q$BB, q$AA)) return(list(k = q$kk, bid = q$BB, ask = q$AA, widen = j))
  }
  list(k = q$kk, bid = q$BB, ask = q$AA, widen = NA_integer_)
}

#' Bounds on tail probabilities and call prices consistent with the quotes
#'
#' Sharp single-maturity bounds implied by bid and ask quotes, a known forward
#' and no arbitrage, using the closed forms of Dillschneider (2026,
#' Propositions 3.1 and 3.3). `tail_bounds()` bounds the left-tail probability
#' \eqn{Q(Y \le a)} with \eqn{Y = S_T/F}; at `a = 0` it returns the set of
#' probabilities of zero equity (the PoD) consistent with the quotes, whose
#' lower bound is zero whenever the lowest quoted strike is positive.
#' `call_bounds()` bounds the normalised call price \eqn{E^Q[(Y-k)^+]} and
#' `atm_iv_bounds()` translates the bounds at \eqn{k=1} into bounds on
#' at-the-money (Black) implied volatility.
#'
#' Quotes are treated as European; with American options and a parity
#' forward the bounds are approximate.
#'
#' @param chain An [ipod_chain()]. Bounds use all quotes, including those with
#'   zero bid; build the chain with `drop_zero_oi = FALSE` and `min_contracts = 0`
#'   to use every quote.
#' @param a Thresholds on the \eqn{S_T/F} scale.
#' @param k Strikes on the \eqn{K/F} scale.
#' @inheritParams arbitrage_free_quotes
#' @return A data frame with the threshold (`a` or `k`), `lower`, `upper` and
#'   `widen` (ticks used). `atm_iv_bounds()` returns a named vector `lower`,
#'   `upper`, `widen`. Bounds are `NA` if the quotes are not arbitrage-free
#'   after `max_widen` ticks.
#' @references Dillschneider, Y. (2026). Option prices and option-implied
#'   probabilities in incomplete markets. SSRN working paper 6976498.
#'   \doi{10.2139/ssrn.6976498}
#' @seealso [lp_bounds()] for general payoffs.
#' @examples
#' ch <- ipod_chain(sim_quotes(sigma = 0.4, pod = 0.01))
#' tail_bounds(ch, c(0, 0.5, 0.8))
#' call_bounds(ch, 1)
#' atm_iv_bounds(ch)
#' @export
tail_bounds <- function(chain, a, tick = 0.01, max_widen = 2) {
  q <- arbitrage_free_quotes(chain, tick, max_widen)
  b <- if (is.na(q$widen)) matrix(NA_real_, length(a), 2)
       else t(vapply(a, function(ai) ltp_bounds(ai, q$k, q$bid, q$ask), c(L = 0, U = 0)))
  data.frame(a = a, lower = b[, 1], upper = b[, 2], widen = q$widen)
}

#' @rdname tail_bounds
#' @export
call_bounds <- function(chain, k, tick = 0.01, max_widen = 2) {
  q <- arbitrage_free_quotes(chain, tick, max_widen)
  if (is.na(q$widen)) return(data.frame(k = k, lower = NA_real_, upper = NA_real_, widen = NA_integer_))
  data.frame(k = k, lower = vapply(k, lb_call, 0, kk = q$k, BB = q$bid, AA = q$ask),
             upper = vapply(k, ub_call, 0, kk = q$k, AA = q$ask), widen = q$widen)
}

#' @rdname tail_bounds
#' @export
atm_iv_bounds <- function(chain, tick = 0.01, max_widen = 2) {
  b <- call_bounds(chain, 1, tick, max_widen)
  c(lower = atm_iv_from_price(b$lower, chain$T), upper = atm_iv_from_price(b$upper, chain$T), widen = b$widen)
}

# ATM Black IV from the normalised at-the-money-forward call c = 2 Phi(s sqrt(T) / 2) - 1
atm_iv_from_price <- function(c, T) {
  ifelse(is.finite(c) & c > 0 & c < 1, 2 * stats::qnorm((c + 1) / 2) / sqrt(T),
         ifelse(is.finite(c) & c <= 0, 0, NA_real_))
}

#' Bounds on general expectations by linear programming
#'
#' Lower and upper bounds of \eqn{E^Q[f(Y)]} over all distributions of
#' \eqn{Y = S_T/F \ge 0} on a discrete support that have mean one and call
#' prices within the (arbitrage-free) quotes. For tail probabilities and call
#' prices the closed forms in [tail_bounds()] are exact and faster; the linear
#' program serves other payoffs and as a check. Requires the \pkg{lpSolve}
#' package.
#'
#' @param chain An [ipod_chain()].
#' @param payoff A function of \eqn{Y} (vectorised).
#' @param support Support points of \eqn{Y}; the strikes, 0 and 1 are added.
#'   Bounds for payoffs that grow quickly in \eqn{Y} depend on its upper end.
#' @inheritParams arbitrage_free_quotes
#' @return A named vector with `lower`, `upper` and `widen`.
#' @examples
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   ch <- ipod_chain(sim_quotes())
#'   lp_bounds(ch, function(y) as.numeric(y <= 0.8))
#'   tail_bounds(ch, 0.8)
#' }
#' @export
lp_bounds <- function(chain, payoff, support = c(seq(0, 3, by = 0.0025), seq(3.05, 10, by = 0.05), seq(10.5, 60, by = 0.5)),
                      tick = 0.01, max_widen = 2) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) stop("`lp_bounds()` requires the 'lpSolve' package.", call. = FALSE)
  q <- arbitrage_free_quotes(chain, tick, max_widen)
  if (is.na(q$widen)) return(c(lower = NA_real_, upper = NA_real_, widen = NA_real_))
  y <- sort(unique(round(c(0, support, q$k[-1], 1), 10)))
  pay <- t(pmax(outer(y, q$k[-1], "-"), 0))
  A <- rbind(rep(1, length(y)), y, pay, pay)
  dir <- c("=", "=", rep(">=", nrow(pay)), rep("<=", nrow(pay)))
  rhs <- c(1, 1, q$bid[-1], q$ask[-1])
  obj <- as.numeric(payoff(y))
  f <- function(s) { z <- lpSolve::lp(s, obj, A, dir, rhs, scale = 4); if (z$status == 0) z$objval else NA_real_ }
  c(lower = f("min"), upper = f("max"), widen = q$widen)
}

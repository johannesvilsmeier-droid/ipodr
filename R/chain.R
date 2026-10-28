# Option chains: input checks, put-call-parity forward, normalisation ---------

.check_quotes <- function(quotes) {
  if (!is.data.frame(quotes)) stop("`quotes` must be a data frame.", call. = FALSE)
  need <- c("strike", "type", "bid", "ask", "rate")
  miss <- setdiff(need, names(quotes))
  if (length(miss)) stop("`quotes` lacks column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (!"ttm" %in% names(quotes)) {
    if (!all(c("date", "expiry") %in% names(quotes)))
      stop("`quotes` needs `ttm` (calendar days) or both `date` and `expiry`.", call. = FALSE)
    quotes$ttm <- as.numeric(as.Date(quotes$expiry) - as.Date(quotes$date))
  }
  ty <- tolower(substr(as.character(quotes$type), 1, 1))
  if (anyNA(ty) || any(!ty %in% c("c", "p")))
    stop("`type` must be \"call\" or \"put\" (or \"C\"/\"P\").", call. = FALSE)
  quotes$call <- ty == "c"
  if (length(unique(quotes$ttm)) != 1L || length(unique(quotes$rate)) != 1L)
    stop("`quotes` must contain a single expiry with one rate; split the data first.", call. = FALSE)
  if ("date" %in% names(quotes) && length(unique(quotes$date)) != 1L)
    stop("`quotes` must contain a single quote date; split the data first.", call. = FALSE)
  if (!is.finite(quotes$ttm[1]) || quotes$ttm[1] <= 0) stop("`ttm` must be positive.", call. = FALSE)
  quotes
}

#' Forward price from put-call parity
#'
#' Estimates the prepaid forward \eqn{F_p = \mathrm{DF}\,F} from call and put
#' mid quotes at common strikes. Parity \eqn{C - P = F_p - \mathrm{DF}\,K}
#' is applied at the `n` strikes with the smallest absolute call-put
#' difference, where early-exercise premia of American options are smallest,
#' and the median is returned. The estimate nets out dividends paid before
#' expiry.
#'
#' @param strike Strikes common to calls and puts.
#' @param call_mid,put_mid Call and put mid quotes at `strike`.
#' @param DF Discount factor to expiry.
#' @param n Number of strikes used.
#' @return The prepaid forward \eqn{F_p}; `NA` if fewer than two strikes are
#'   available. Divide by `DF` for the forward.
#' @examples
#' K <- c(90, 100, 110)
#' parity_forward(K, call_mid = c(14, 7.5, 3.2), put_mid = c(3.6, 7.0, 12.6), DF = 0.99)
#' @export
parity_forward <- function(strike, call_mid, put_mid, DF, n = 3) {
  ok <- is.finite(strike) & is.finite(call_mid) & is.finite(put_mid)
  if (sum(ok) < 2) return(NA_real_)
  strike <- strike[ok]; cm <- call_mid[ok]; pm <- put_mid[ok]
  i <- order(abs(cm - pm))[seq_len(min(n, length(cm)))]
  Fp <- stats::median(cm[i] - pm[i] + DF * strike[i])
  if (!is.finite(Fp) || Fp <= 0) NA_real_ else Fp
}

#' Prepare and normalise one option chain
#'
#' Converts the quotes of one underlying, one quote date and one expiry into
#' the normalised form used by all estimators and bounds in the package. Prices
#' are expressed relative to the prepaid forward \eqn{F_p} and strikes relative
#' to the forward \eqn{F = F_p/\mathrm{DF}}, so that \eqn{Y = S_T/F} has
#' risk-neutral mean one and a call with strike \eqn{K} is worth
#' \eqn{F_p\,E^Q[(Y-k)^+]} with \eqn{k = K/F}.
#'
#' Two forward conventions are available:
#' * `"parity"` (default): the forward is estimated by [parity_forward()];
#'   out-of-the-money puts (\eqn{K<F}) are converted to calls by parity and
#'   combined with out-of-the-money calls (\eqn{K \ge F}).
#' * `"spot"`: calls only, with the zero-payout forward \eqn{F^0 = S/\mathrm{DF}}.
#'   Payouts before expiry are then absorbed by the tolerance on the forward
#'   constraint in [ipod_fit()].
#'
#' Quotes are treated as European. Exchange-listed single-stock options are
#' American; using out-of-the-money quotes and the parity forward limits, but
#' does not remove, the resulting bias (put-based left-tail probabilities may
#' be slightly too high).
#'
#' @param quotes A data frame in long format with one row per quote and columns
#'   `strike`, `type` (`"call"`/`"put"` or `"C"`/`"P"`), `bid`, `ask`, `rate`
#'   (continuously compounded annual rate) and either `ttm` (calendar days to
#'   expiry) or `date` and `expiry`. Optional: `oi` (open interest) and `spot`.
#'   All rows must refer to one quote date and one expiry.
#' @param forward `"parity"` or `"spot"`, see Details.
#' @param spot Spot price; required for `forward = "spot"` unless `quotes`
#'   has a `spot` column.
#' @param drop_zero_oi Drop quotes with zero open interest (missing open
#'   interest is kept). The parity forward is computed before this filter,
#'   from all quotes with a positive bid.
#' @param oi_min,min_contracts The chain is eligible (`ok = TRUE`) if at least
#'   `min_contracts` retained quotes have open interest above `oi_min`. Without
#'   an `oi` column, quotes with a positive bid are counted.
#' @return An object of class `ipod_chain`: a list with `T` (years), `DF`,
#'   `F`, `Fp`, `forward`, `kmin` and `kmax` (lowest and highest strike with a
#'   positive bid, on the \eqn{k=K/F} scale), `ok`, `flag` and `quotes`, a data
#'   frame with `strike`, `call`, `bid`, `ask`, `oi`, `k` and the call-unit
#'   quotes `cb`, `ca`, `cm` (bid, ask, mid divided by \eqn{F_p}, puts
#'   converted by parity).
#' @seealso [ipod_fit()], [tail_bounds()], [sim_quotes()].
#' @examples
#' q <- sim_quotes(ttm = 182, sigma = 0.4)
#' ch <- ipod_chain(q)
#' ch
#' @export
ipod_chain <- function(quotes, forward = c("parity", "spot"), spot = NULL,
                       drop_zero_oi = TRUE, oi_min = 5, min_contracts = 4) {
  forward <- match.arg(forward)
  q <- .check_quotes(quotes)
  q <- q[is.finite(q$strike) & q$strike > 0 & is.finite(q$bid) & is.finite(q$ask) &
           q$bid >= 0 & q$ask > 0 & q$ask >= q$bid, , drop = FALSE]
  q <- q[!duplicated(q[, c("strike", "call")]), , drop = FALSE]
  has_oi <- "oi" %in% names(q)
  if (!has_oi) q$oi <- NA_real_
  T <- q$ttm[1] / 365
  DF <- exp(-q$rate[1] * T)
  out <- list(T = T, DF = DF, F = NA_real_, Fp = NA_real_, forward = forward,
              kmin = NA_real_, kmax = NA_real_, ok = FALSE, flag = "", quotes = NULL)
  class(out) <- "ipod_chain"
  q$mid <- (q$bid + q$ask) / 2
  if (forward == "parity") {
    cc <- q[q$call & q$bid > 0, c("strike", "mid")]
    pp <- q[!q$call & q$bid > 0, c("strike", "mid")]
    m <- merge(cc, pp, by = "strike", suffixes = c("_c", "_p"))
    Fp <- if (nrow(m) >= 2) parity_forward(m$strike, m$mid_c, m$mid_p, DF) else NA_real_
    if (!is.finite(Fp)) { out$flag <- "no parity forward"; return(out) }
    F <- Fp / DF
    o <- q[(!q$call & q$strike < F) | (q$call & q$strike >= F), , drop = FALSE]
  } else {
    S <- if (!is.null(spot)) spot else if ("spot" %in% names(q)) q$spot[1] else NA_real_
    if (!is.finite(S) || S <= 0) stop("`forward = \"spot\"` needs a positive `spot`.", call. = FALSE)
    Fp <- S; F <- S / DF
    o <- q[q$call, , drop = FALSE]
  }
  if (drop_zero_oi) o <- o[is.na(o$oi) | o$oi > 0, , drop = FALSE]
  o <- o[order(o$strike), , drop = FALSE]
  out$F <- F; out$Fp <- Fp
  if (!nrow(o)) { out$flag <- "no quotes"; return(out) }
  k <- o$strike / F
  conv <- function(v) ifelse(o$call, v / Fp, v / Fp + 1 - k)
  qd <- data.frame(strike = o$strike, call = o$call, bid = o$bid, ask = o$ask, oi = o$oi, k = k,
                   cb = conv(o$bid), ca = conv(o$ask), cm = conv(o$mid))
  out$quotes <- qd
  pos <- qd$bid > 0
  if (any(pos)) { out$kmin <- min(qd$k[pos]); out$kmax <- max(qd$k[pos]) }
  n_liq <- if (has_oi) sum(qd$oi > oi_min, na.rm = TRUE) else sum(pos)
  if (n_liq < min_contracts) { out$flag <- "too few contracts"; return(out) }
  if (max(qd$k) < 0.3 || min(qd$k) > 3) { out$flag <- "strike scale mismatch"; return(out) }
  out$ok <- TRUE
  out
}

#' @export
print.ipod_chain <- function(x, ...) {
  cat(sprintf("<ipod_chain> %s forward, T = %.3f years, DF = %.4f, F = %.4g\n", x$forward, x$T, x$DF, x$F))
  if (!is.null(x$quotes))
    cat(sprintf("  %d quotes (%d calls, %d puts), k in [%.3f, %.3f]\n", nrow(x$quotes),
                sum(x$quotes$call), sum(!x$quotes$call), x$kmin, x$kmax))
  cat(sprintf("  eligible: %s%s\n", x$ok, if (nzchar(x$flag)) paste0(" (", x$flag, ")") else ""))
  invisible(x)
}

#' Select the expiry closest to a target maturity
#'
#' @param quotes Quotes of one underlying and one quote date (long format, see
#'   [ipod_chain()]), possibly with several expiries.
#' @param target Target maturity in calendar days.
#' @param range Admissible maturities in calendar days.
#' @param both If `TRUE`, only expiries with both calls and puts are eligible.
#' @return The rows of `quotes` belonging to the selected expiry (zero rows if
#'   none is eligible).
#' @examples
#' q <- rbind(sim_quotes(ttm = 120), sim_quotes(ttm = 200))
#' unique(select_expiry(q, target = 182)$ttm)
#' @export
select_expiry <- function(quotes, target = 182, range = c(121, 273), both = TRUE) {
  if (!"ttm" %in% names(quotes)) quotes$ttm <- as.numeric(as.Date(quotes$expiry) - as.Date(quotes$date))
  call <- tolower(substr(as.character(quotes$type), 1, 1)) == "c"
  ex <- sort(unique(quotes$ttm[quotes$ttm >= range[1] & quotes$ttm <= range[2]]))
  if (both) ex <- ex[ex %in% quotes$ttm[call] & ex %in% quotes$ttm[!call]]
  if (!length(ex)) return(quotes[0, , drop = FALSE])
  quotes[quotes$ttm == ex[which.min(abs(ex - target))], , drop = FALSE]
}

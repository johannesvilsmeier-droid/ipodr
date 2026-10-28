#' ipodr: Option-Implied Distress Measures and Their Identified Sets
#'
#' Estimates option-implied risk-neutral distributions of a single stock and
#' the distress measures derived from them, and reports which of these
#' quantities the option quotes actually determine.
#'
#' Workflow:
#' 1. [ipod_chain()] prepares one quote date and expiry from a long data frame
#'    of quotes (parity forward, out-of-the-money calls and puts, normalised
#'    to \eqn{Y = S_T/F}).
#' 2. [ipod_fit()] fits the penalised maximum-entropy distribution with a
#'    default barrier; [tail_prob()], [rnd_summary()] and [adp()] read tail
#'    probabilities from it; [ipod_profile()] and [pod_select()] show how the
#'    PoD depends on the barrier.
#' 3. [tail_bounds()], [call_bounds()], [atm_iv_bounds()] and [lp_bounds()]
#'    give the ranges consistent with the quotes and no arbitrage.
#' 4. [atm_iv()], [lognormal_tail()], [drawdown()] and [dtd()] provide the
#'    benchmarks; [label_events()], [auc_test()], [auc_bounds()] and
#'    [first_alarm()] evaluate indicators against distress events.
#'
#' All probabilities are risk-neutral. They are not physical default
#' probabilities unless calibrated.
#'
#' @keywords internal
#' @importFrom utils head tail
"_PACKAGE"

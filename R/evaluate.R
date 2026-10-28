# Evaluation: event windows, AUC, cluster bootstrap, AUC bounds, alarms ---------

#' Area under the ROC and precision-recall curves
#'
#' `auc()` is the rank-based (Mann-Whitney) area under the ROC curve, with ties
#' counted one half; `pr_auc()` is average precision. Higher scores are taken
#' to indicate the positive class; negate scores where low values signal risk.
#'
#' @param score Numeric scores.
#' @param label Logical (or 0/1) labels; `TRUE` = positive (event window).
#' @return A number (`NA` if one class is empty).
#' @examples
#' auc(c(0.1, 0.4, 0.35, 0.8), c(FALSE, FALSE, TRUE, TRUE))
#' pr_auc(c(0.1, 0.4, 0.35, 0.8), c(FALSE, FALSE, TRUE, TRUE))
#' @export
auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label); s <- score[ok]; y <- as.logical(label[ok])
  n1 <- sum(y); n0 <- sum(!y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(s)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @rdname auc
#' @export
pr_auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label); s <- score[ok]; y <- as.logical(label[ok])
  if (!any(y)) return(NA_real_)
  y <- y[order(-s)]
  sum((cumsum(y) / seq_along(y))[y]) / sum(y)
}

#' Label event windows
#'
#' Marks the observations of event banks in the `horizon` days before their
#' event as positives and keeps, as negatives, observations of control banks
#' (no event) within the date range spanned by the positives. Earlier
#' observations of event banks and all observations on or after an event are
#' dropped.
#'
#' @param data A data frame with one row per bank and date.
#' @param horizon Horizon in days.
#' @param bank,date,event Column names: bank identifier, observation date
#'   (`Date`) and event date (`Date`, `NA` for control banks).
#' @return The retained rows of `data` with a logical column `pos`.
#' @examples
#' d <- data.frame(bank = rep(c("A", "B"), each = 4),
#'                 date = rep(as.Date("2023-01-01") + c(0, 30, 60, 90), 2),
#'                 event = rep(as.Date(c("2023-03-15", NA)), each = 4))
#' label_events(d, horizon = 60)
#' @export
label_events <- function(data, horizon, bank = "bank", date = "date", event = "event") {
  ev <- as.Date(data[[event]]); dt <- as.Date(data[[date]])
  is_ev <- !is.na(ev)
  pos <- is_ev & dt < ev & dt >= ev - horizon
  if (!any(pos)) stop("No observation falls in an event window.", call. = FALSE)
  rng <- range(dt[pos])
  keep <- pos | (!is_ev & dt >= rng[1] & dt <= rng[2])
  out <- data[keep, , drop = FALSE]
  out$pos <- pos[keep]
  out
}

#' Compare indicators by AUC with a bank-cluster bootstrap
#'
#' For each score column, reports the pooled AUC, average precision and a
#' bank-level AUC (bank means), with bootstrap intervals that resample event
#' banks and control banks separately (clusters of all their observations).
#' If `ref` is given, the AUC difference to the reference indicator is computed
#' on the observations where both are available.
#'
#' @param data Output of [label_events()] (needs `pos` and a bank column).
#' @param scores Names of score columns (higher = riskier).
#' @param ref Optional name of the reference score column.
#' @param bank Name of the bank column.
#' @param B Number of bootstrap replications.
#' @param level Confidence level of the percentile intervals.
#' @return A data frame with one row per score: `score`, `n_pos`, `n_neg`,
#'   `auc`, `auc_lower`, `auc_upper`, `pr_auc`, `base_rate`, `bank_auc` and,
#'   with `ref`, `diff`, `diff_lower`, `diff_upper` (bootstrap median and
#'   percentile interval).
#' @examples
#' set.seed(1)
#' d <- data.frame(bank = rep(letters[1:8], each = 20), pos = rep(c(TRUE, FALSE), c(40, 120)))
#' d$x1 <- rnorm(160, mean = d$pos); d$x2 <- rnorm(160, mean = 0.5 * d$pos)
#' auc_test(d, c("x1", "x2"), ref = "x2", B = 200)
#' @export
auc_test <- function(data, scores, ref = NULL, bank = "bank", B = 500, level = 0.90) {
  if (!"pos" %in% names(data)) stop("`data` needs a logical column `pos`; see label_events().", call. = FALSE)
  y <- as.logical(data$pos); bk <- data[[bank]]
  banks1 <- unique(bk[y]); banks0 <- unique(bk[!y])
  idx <- split(seq_len(nrow(data)), bk)
  boot <- replicate(B, {
    b <- c(banks1[sample.int(length(banks1), replace = TRUE)], banks0[sample.int(length(banks0), replace = TRUE)])
    unlist(idx[as.character(b)], use.names = FALSE)
  }, simplify = FALSE)
  al <- (1 - level) / 2
  qq <- function(v, p) if (all(is.na(v))) NA_real_ else stats::quantile(v, p, na.rm = TRUE, names = FALSE)
  do.call(rbind, lapply(scores, function(cl) {
    s <- data[[cl]]; ok <- is.finite(s)
    bm <- tapply(s[ok], bk[ok], mean); by <- tapply(y[ok], bk[ok], any)
    bs <- vapply(boot, function(i) { i <- i[ok[i]]; auc(s[i], y[i]) }, numeric(1))
    out <- data.frame(score = cl, n_pos = sum(y[ok]), n_neg = sum(!y[ok]), auc = auc(s[ok], y[ok]),
                      auc_lower = qq(bs, al), auc_upper = qq(bs, 1 - al), pr_auc = pr_auc(s[ok], y[ok]),
                      base_rate = mean(y[ok]), bank_auc = auc(as.numeric(bm), as.logical(by)))
    if (!is.null(ref)) {
      r <- data[[ref]]
      bd <- if (cl == ref) NA_real_ else vapply(boot, function(i) {
        i <- i[ok[i] & is.finite(r[i])]; auc(s[i], y[i]) - auc(r[i], y[i]) }, numeric(1))
      both <- ok & is.finite(r)
      out$diff <- if (cl == ref) NA_real_ else auc(s[both], y[both]) - auc(r[both], y[both])
      out$diff_lower <- qq(bd, al); out$diff_upper <- qq(bd, 1 - al)
    }
    out
  }))
}

#' AUC bounds consistent with option prices
#'
#' If each observation's score is only known to lie in an interval
#' \eqn{[L_i, U_i]} (for example the tail-probability bounds from
#' [tail_bounds()]) and the intervals vary independently, the sharp bounds on
#' the AUC are attained with events at their lower and controls at their upper
#' bounds (lower AUC) and the reverse (upper AUC).
#'
#' @param lower,upper Lower and upper bounds of the scores.
#' @param label Logical labels (`TRUE` = positive).
#' @return A named vector with `lower` and `upper`.
#' @references Bondarenko, O., Dillschneider, Y., Schneider, P. and Trojani, F.
#'   (2026). What can you really tell from option prices? Working paper.
#' @examples
#' auc_bounds(lower = c(0.1, 0.2, 0.5, 0.6), upper = c(0.3, 0.4, 0.7, 0.9),
#'            label = c(FALSE, FALSE, TRUE, TRUE))
#' @export
auc_bounds <- function(lower, upper, label) {
  label <- as.logical(label)
  c(lower = auc(ifelse(label, lower, upper), label), upper = auc(ifelse(label, upper, lower), label))
}

#' First alarm date
#'
#' The first date on which a score has exceeded `threshold` on `run`
#' consecutive observations.
#'
#' @param date Observation dates (sorted).
#' @param score Scores.
#' @param threshold Alarm threshold.
#' @param run Number of consecutive observations required.
#' @return The alarm date (`NA` if there is none).
#' @examples
#' d <- as.Date("2023-01-01") + 0:9
#' first_alarm(d, c(1, 3, 3, 1, 3, 3, 3, 3, 3, 1), threshold = 2, run = 5)
#' @export
first_alarm <- function(date, score, threshold, run = 5) {
  ok <- is.finite(score); date <- date[ok]; score <- score[ok]
  if (length(score) < run) return(date[NA_integer_][1])
  r <- stats::filter(as.numeric(score > threshold), rep(1, run), sides = 1)
  w <- which(r >= run)
  if (length(w)) date[w[1]] else date[NA_integer_][1]
}

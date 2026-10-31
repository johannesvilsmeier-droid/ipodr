test_that("auc and pr_auc", {
  expect_equal(auc(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(auc(c(1, 1, 1, 1), c(FALSE, FALSE, TRUE, TRUE)), 0.5)
  expect_equal(auc(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, TRUE)), 0)
  expect_true(is.na(auc(1:3, c(TRUE, TRUE, TRUE))))
  expect_equal(pr_auc(c(4, 3, 2, 1), c(TRUE, TRUE, FALSE, FALSE)), 1)
})

test_that("event labels exclude post-event and early event-bank rows", {
  d <- data.frame(bank = rep(c("A", "B"), each = 6),
                  date = rep(as.Date("2023-01-01") + c(0, 30, 60, 90, 120, 150), 2),
                  event = rep(as.Date(c("2023-04-15", NA)), each = 6))
  L <- label_events(d, horizon = 60)
  expect_true(all(L$date[L$pos] < as.Date("2023-04-15")))
  expect_equal(sum(L$pos), 2)
  expect_true(all(L$bank[!L$pos] == "B"))
  expect_true(all(L$date[!L$pos] >= min(L$date[L$pos]) & L$date[!L$pos] <= max(L$date[L$pos])))
})

test_that("auc_test reports differences and intervals", {
  set.seed(2)
  d <- data.frame(bank = rep(letters[1:10], each = 10), pos = rep(c(TRUE, FALSE), c(30, 70)))
  d$x1 <- rnorm(100, mean = 2 * d$pos); d$x2 <- rnorm(100)
  r <- auc_test(d, c("x1", "x2"), ref = "x2", B = 100)
  expect_equal(r$score, c("x1", "x2"))
  expect_true(r$auc[1] > 0.8 && r$auc_lower[1] <= r$auc[1] && r$auc[1] <= r$auc_upper[1])
  expect_true(r$diff_lower[1] > 0)
  expect_true(is.na(r$diff[2]))
})

test_that("auc_bounds corners", {
  b <- auc_bounds(c(0, 0, 1, 1), c(0.5, 0.5, 2, 2), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(unname(b), c(1, 1))
  b <- auc_bounds(c(0, 0, 0, 0), c(1, 1, 1, 1), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(unname(b), c(0, 1))
})

test_that("first_alarm, drawdown, dtd, adp, const_maturity", {
  d <- as.Date("2023-01-01") + 0:9
  expect_equal(first_alarm(d, c(1, 3, 3, 1, 3, 3, 3, 3, 3, 1), 2), d[9])
  expect_true(is.na(first_alarm(d, rep(1, 10), 2)))
  expect_equal(drawdown(c(100, 50, 75), window = 3), c(0, 0.5, 0.25))
  expect_equal(dtd(0, 0.2, T = 1, theta = 1), -0.1)
  expect_equal(adp(peak = 200, forward = 100, theta = 0.5, sigma = 0.3, T = 0.5), lognormal_tail(1, 0.3, 0.5))
  expect_equal(const_maturity(c(90, 270), c(1, 3), 180), 2)
  expect_true(is.na(const_maturity(c(30, 60), c(1, 2), 182)))
})

test_that("realised_vol, padp, peer_percentile", {
  set.seed(4)
  p <- 100 * exp(cumsum(rnorm(400, sd = 0.02)))
  rv <- realised_vol(p)
  expect_true(all(is.na(rv[1:251])))
  expect_equal(rv[400], sd(diff(log(p[148:400]))) * sqrt(252))
  q <- p; q[300] <- NA
  expect_true(is.na(realised_vol(q)[300]))
  expect_equal(realised_vol(q)[301], sd(diff(log(q[c(49:299, 301)]))) * sqrt(252))
  s <- p; s[100:300] <- s[100]
  expect_true(is.na(realised_vol(s)[300]))
  pd <- padp(p)
  pk <- vapply(seq_along(p), function(i) max(p[max(1, i - 251):i]), 0)
  expect_equal(pd[400], adp(peak = pk[400], forward = p[400], theta = 0.5, sigma = rv[400], T = 0.5))
  expect_equal(pd[400], pnorm(-dtd(drawdown(p)[400], rv[400], T = 0.5, theta = 0.5)))
  expect_true(all(is.na(pd[1:251])) && all(pd[252:400] > 0 & pd[252:400] < 1))
  expect_true(padp(p, theta = 0.3)[400] < pd[400])
  expect_equal(padp(p, sigma = rep(0.3, 400))[10], pnorm(-dtd(drawdown(p)[10], 0.3, 0.5, 0.5)))
  expect_equal(peer_percentile(c(1, 2, 3, NA, 5, 5), c(1, 1, 1, 1, 2, 2), min_n = 2),
               c(1, 3, 5, NA, 2, 2) / c(6, 6, 6, 1, 4, 4))
  expect_true(all(is.na(peer_percentile(1:3, 1:3, min_n = 2))))
})

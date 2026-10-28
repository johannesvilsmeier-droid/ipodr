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

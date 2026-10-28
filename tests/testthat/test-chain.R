test_that("parity forward and normalisation recover the simulated forward", {
  q <- sim_quotes(forward = 50, rate = 0.04, sigma = 0.3)
  ch <- ipod_chain(q)
  expect_true(ch$ok)
  expect_equal(ch$F, 50, tolerance = 0.002)
  expect_true(all(ch$quotes$call == (ch$quotes$k >= 1)))
  # spot convention: calls only, zero-payout forward
  cs <- ipod_chain(q, forward = "spot")
  expect_true(all(cs$quotes$call))
  expect_equal(cs$F, 50, tolerance = 1e-10)
})

test_that("input checks", {
  q <- sim_quotes()
  expect_error(ipod_chain(q[, -which(names(q) == "bid")]), "bid")
  expect_error(ipod_chain(rbind(q, transform(q, ttm = 90))), "single expiry")
  q2 <- q; q2$oi <- 0
  expect_false(ipod_chain(q2)$ok)
  expect_identical(nrow(select_expiry(rbind(sim_quotes(ttm = 30), q), target = 182)), nrow(q))
})

test_that("Black implied volatility round trip", {
  k <- c(0.6, 0.9, 1, 1.2)
  c <- black_call(k, 0.45, 0.5)
  expect_equal(implied_vol(c * 0.98 * 80, k * 80, 80, 0.5, DF = 0.98), rep(0.45, 4), tolerance = 1e-8)
  p <- c * 0.98 * 80 - 0.98 * (80 - k * 80)
  expect_equal(implied_vol(p, k * 80, 80, 0.5, DF = 0.98, type = "put"), rep(0.45, 4), tolerance = 1e-8)
  expect_true(is.na(implied_vol(0, 120, 100, 0.5)))
  expect_equal(atm_iv(ipod_chain(sim_quotes(sigma = 0.3, spread_abs = 0, spread_rel = 0, tick = 1e-6))), 0.3, tolerance = 1e-3)
})

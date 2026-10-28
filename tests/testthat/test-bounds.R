test_that("bounds at the strikes equal exact prices (Dillschneider example)", {
  kx <- c(0.5, 0.75, 1, 1.25, 1.5); cx <- c(0.55, 0.35, 0.2, 0.1, 0.05)
  kk <- c(0, kx); BB <- c(1, cx); AA <- c(1, cx)
  expect_equal(vapply(kx, ipodr:::lb_call, 0, kk = kk, BB = BB, AA = AA), cx)
  expect_equal(vapply(kx, ipodr:::ub_call, 0, kk = kk, AA = AA), cx)
  # PoD bounds: [0, 1 + slope of the first segment] = [0, 1 - (1 - 0.55)/0.5]
  b0 <- ipodr:::ltp_bounds(0, kk, BB, AA)
  expect_equal(unname(b0), c(0, 1 - 0.45 / 0.5), tolerance = 1e-8)
})

test_that("tail bounds contain the truth and the PoD lower bound is zero", {
  q <- sim_quotes(sigma = 0.5, pod = 0.02)
  ch <- ipod_chain(q, drop_zero_oi = FALSE, min_contracts = 0)
  a <- c(0, 0.3, 0.5, 0.8, 1)
  b <- tail_bounds(ch, a)
  tr <- attr(q, "tail")(a)
  expect_true(all(b$lower <= tr + 1e-8 & tr <= b$upper + 1e-8))
  expect_equal(b$lower[1], 0)
  expect_true(all(b$lower <= b$upper))
  ivb <- atm_iv_bounds(ch)
  iv_true <- ipodr:::iv_norm(black_call(0.98, 0.5, ch$T), 1, ch$T)
  expect_true(ivb[["lower"]] <= iv_true && iv_true <= ivb[["upper"]])
  cb <- call_bounds(ch, c(0.8, 1, 1.2))
  expect_true(all(cb$lower <= black_call(c(0.8, 1, 1.2) * 0.98, 0.5, 0.5) + 1e-8))
})

test_that("closed-form bounds equal the linear program", {
  skip_if_not_installed("lpSolve")
  ch <- ipod_chain(sim_quotes(sigma = 0.4, spread_rel = 0.1))
  for (a in c(0.5, 0.8)) {
    sup <- sort(c(seq(0, 3, by = 0.0025), seq(3.05, 10, by = 0.05), seq(10.5, 60, by = 0.5), a, a + 1e-6))
    lp <- lp_bounds(ch, function(y) as.numeric(y <= a), support = sup)
    cf <- tail_bounds(ch, a)
    expect_equal(lp[["lower"]], cf$lower, tolerance = 2e-3)
    expect_equal(lp[["upper"]], cf$upper, tolerance = 2e-3)
  }
  lc <- lp_bounds(ch, function(y) pmax(y - 1, 0))
  cc <- call_bounds(ch, 1)
  expect_equal(lc[["lower"]], cc$lower, tolerance = 1e-4)
  expect_equal(lc[["upper"]], cc$upper, tolerance = 1e-4)
})

test_that("ME tail probabilities lie inside the bounds", {
  ch <- ipod_chain(sim_quotes(sigma = 0.4))
  fit <- ipod_fit(ch, n_grid = 800)
  b <- tail_bounds(ch, c(0.6, 0.8))
  p <- tail_prob(fit, c(0.6, 0.8))
  expect_true(all(p >= b$lower - 0.01 & p <= b$upper + 0.01))
})

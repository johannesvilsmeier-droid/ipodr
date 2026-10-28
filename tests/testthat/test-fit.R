test_that("solver converges and fits within tolerance", {
  ch <- ipod_chain(sim_quotes(sigma = 0.4))
  fit <- ipod_fit(ch, d = 0.1, n_grid = 800)
  expect_true(fit$converged)
  expect_gt(mean(abs(fit$fitted - fit$observed) <= fit$sigma + 1e-8), 0.9)
  expect_equal(sum(fit$p) + fit$pod, 1, tolerance = 1e-10)
  expect_lt(abs(fit$stock_error), 0.01)
})

test_that("tail probabilities inside the strike range are close to the truth", {
  for (pod in c(0, 0.03)) {
    q <- sim_quotes(sigma = 0.45, pod = pod)
    fit <- ipod_fit(ipod_chain(q), d = 0.1, n_grid = 1000)
    a <- c(0.5, 0.7, 0.8)
    expect_lt(max(abs(tail_prob(fit, a) - attr(q, "tail")(a))), 0.03)
  }
})

test_that("features, quantiles and profile are consistent", {
  ch <- ipod_chain(sim_quotes())
  fit <- ipod_fit(ch, n_grid = 500)
  s <- rnd_summary(fit)
  expect_named(s, c("pod", "p_eq80", "p_eq70", "p_eq50", "p_eq25", "q05", "iqr", "bowley", "sd", "kl", "kmin"))
  expect_true(all(diff(tail_prob(fit, c(0, 0.25, 0.5, 0.8, 1))) >= 0))
  q50 <- rnd_quantile(fit, 0.5)
  expect_gte(tail_prob(fit, q50 + 1e-6), 0.5)
  pr <- ipod_profile(ch, d_grid = c(0.02, 0.05, 0.1, 0.2, 0.4), n_grid = 500)
  expect_equal(nrow(pr), 5)
  expect_true(all(pr$converged))
  # tail probabilities inside the strike range barely move with the barrier
  expect_lt(diff(range(pr$p_eq80)), 0.02)
  for (r in c("mean", "median", "absolute", "maxent", "fixed")) {
    s <- pod_select(pr, r)
    expect_true(is.finite(s$pod) && s$pod >= 0 && s$pod <= 1)
  }
  expect_equal(pod_select(pr, "fixed", d = 0.1)$pod, pr$pod[3], tolerance = 1e-8)
})

test_that("ineligible chains are rejected", {
  q <- sim_quotes(); q$oi <- 1
  expect_error(ipod_fit(ipod_chain(q)), "not eligible")
})

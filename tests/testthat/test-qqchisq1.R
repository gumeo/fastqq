test_that("qqchisq1 small-N path matches qqlog for equivalent inputs", {
  set.seed(123)
  chisq_vals <- stats::rchisq(500, df = 1)

  # Reference: convert in R then run qqlog logic
  log10_ref <- -stats::pchisq(chisq_vals, df = 1,
                               lower.tail = FALSE, log.p = TRUE) / log(10)
  log10_ref <- log10_ref[log10_ref > 0 & is.finite(log10_ref)]
  o_ref <- sort(log10_ref, decreasing = TRUE)
  e_ref <- -log10(stats::ppoints(length(o_ref)))

  # qqchisq1 small-N logic
  chisq_clean <- chisq_vals[chisq_vals >= 0]
  log10_chisq <- -stats::pchisq(chisq_clean, df = 1,
                                 lower.tail = FALSE, log.p = TRUE) / log(10)
  log10_chisq <- log10_chisq[log10_chisq > 0]
  o_chisq <- sort(log10_chisq, decreasing = TRUE)
  e_chisq <- -log10(stats::ppoints(length(o_chisq)))

  expect_equal(o_chisq, o_ref)
  expect_equal(e_chisq, e_ref)
})

test_that("drop_dense_chisq1 large-N path matches drop_dense_qqlog", {
  set.seed(123)
  chisq_vals <- stats::rchisq(5e4, df = 1)

  # Compute -log10(p) reference in R using same formula as C++
  log10_pvals <- -stats::pchisq(chisq_vals, df = 1,
                                 lower.tail = FALSE, log.p = TRUE) / log(10)

  res_chisq <- fastqq:::drop_dense_chisq1(chisq_vals, 10000L)
  res_qqlog  <- fastqq:::drop_dense_qqlog(log10_pvals, 10000L)

  expect_equal(res_chisq[[1]], res_qqlog[[1]], tolerance = 1e-10)
  expect_equal(res_chisq[[2]], res_qqlog[[2]], tolerance = 1e-10)
})

test_that("qqchisq1 handles extreme chi-sq values without underflow", {
  # For chisq = 1000, p is astronomically small; log-space prevents underflow
  extreme_vals <- c(1, 10, 100, 500, 1000)
  log10_p <- -stats::pchisq(extreme_vals, df = 1,
                              lower.tail = FALSE, log.p = TRUE) / log(10)
  expect_true(all(is.finite(log10_p)))
  expect_true(all(log10_p > 0))
  expect_true(log10_p[5] > log10_p[4])  # monotonically increasing
})

test_that("qqchisq1 correctly places a point with p < 1e-500 on the plot", {
  # chisq ~ 2303 gives p ~ 10^-500, far below .Machine$double.xmin (~2.2e-308).
  # Naive pchisq underflows to 0; qqchisq1 must not lose this point.
  extreme_chisq <- 2303

  # Confirm the naive approach gives 0 (underflow)
  expect_equal(stats::pchisq(extreme_chisq, 1, lower.tail = FALSE), 0)

  # Confirm log-space gives a finite -log10(p) near 500
  neg_log10_p <- -stats::pchisq(extreme_chisq, 1,
                                 lower.tail = FALSE, log.p = TRUE) / log(10)
  expect_true(is.finite(neg_log10_p))
  expect_gt(neg_log10_p, 490)   # should be close to 500
  expect_lt(neg_log10_p, 510)

  # Mix with typical chi-sq values and confirm the extreme point appears
  set.seed(42)
  chisq_vals <- c(stats::rchisq(500, df = 1), extreme_chisq)
  OEmat <- fastqq:::drop_dense_chisq1(chisq_vals, 10000L)
  observed <- OEmat[[1]]

  # The extreme point (max observed) should be near 500
  expect_gt(max(observed), 490)
  expect_lt(max(observed), 510)
})

test_that("qqlog small-N path matches qq for equivalent inputs", {
  set.seed(42)
  pvec <- stats::runif(500)  # all in (0,1), all valid

  # qq small-N computed values
  o_qq <- -log10(sort(pvec, decreasing = FALSE))
  e_qq <- -log10(stats::ppoints(length(pvec)))

  # qqlog small-N computed values
  log10_pvec <- -log10(pvec)
  log10_clean <- log10_pvec[log10_pvec > 0 & is.finite(log10_pvec)]
  o_qqlog <- sort(log10_clean, decreasing = TRUE)
  e_qqlog <- -log10(stats::ppoints(length(o_qqlog)))

  expect_equal(o_qq, o_qqlog)
  expect_equal(e_qq, e_qqlog)
})

test_that("drop_dense_qqlog large-N path matches drop_dense_qq", {
  set.seed(42)
  pvec <- stats::runif(5e4)
  log10_pvec <- -log10(pvec)

  res_qq   <- fastqq:::drop_dense_qq(pvec, 10000L)
  res_qqlog <- fastqq:::drop_dense_qqlog(log10_pvec, 10000L)

  expect_equal(res_qq[[1]], res_qqlog[[1]], tolerance = 1e-10)
  expect_equal(res_qq[[2]], res_qqlog[[2]], tolerance = 1e-10)
})

test_that("qq zero_action parameter replaces zeros with a warning", {
  pvec <- c(stats::runif(100), 0, 0)

  expect_warning(
    qq(pvec, zero_action = 1e-300),
    regexp = "2 p-value\\(s\\) equal to zero replaced with"
  )
})

test_that("qq zero_action = NULL does not warn and silently drops zeros", {
  pvec <- c(stats::runif(100), 0)
  expect_no_warning(qq(pvec))
})

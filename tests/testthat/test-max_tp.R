
test_that("max_tp", {

  # setup -------------------------------------------------------------------

  s <- round(runif(1, min = 15, max = 20))
  b <- round(runif(1, min = 1, max = round(s / 2)))
  l <- s - b
  patch <- round(runif(1, min = 1, max = 10))
  theta <- runif(1, min = 0, max = 10)
  A <- ppm(s, b, l, theta)
  alpha0 <- alpha <- to_alpha(A,
                              attack = list(min = 1, max = 1),
                              convert = list(min = 1, max = 1))
  alpha0[lower.tri(alpha0, diag = TRUE)] <- 0

  x <- n2 <- n1 <- n0 <- rep(0, s)
  n1[attr(alpha, "tp") == 1] <- 1
  n2[attr(alpha, "tp") <= 2] <- 1

  # should be >= n2
  n <- c(rep(1, b), rbinom(s - b, 1, 0.5))

  # should be = 1
  x[attr(alpha, "tp") == 1] <- 1
  x[attr(alpha, "tp") == 3] <- 1

  # test --------------------------------------------------------------------

  expect_equal(max_tp(s, n = n0, alpha = alpha), 0)
  expect_equal(max_tp(s, n = n1, alpha = alpha), 1)
  expect_equal(max_tp(s, n = n2, alpha = alpha), 2)
  expect_equal(max_tp(s, n = x, alpha = alpha), 1)

  expect_gte(max_tp(s, n = n, alpha = alpha),
             max_tp(s, n = n2, alpha = alpha))
})


test_that("ppm/to_alpha", {

  # setup -------------------------------------------------------------------

  s <- round(runif(1, min = 5, max = 20))
  b <- round(runif(1, min = 1, max = s - 1))
  l <- round(runif(1, min = s - b, max = sum(b:(s - 1))))
  theta <- runif(1, min = 0, max = 10)
  A <- ppm(s, b, l, theta)
  alpha <- to_alpha(A)

  mu_l <- mean(sapply(1:1000, function(x) sum(ppm(s, b, l, theta))))

  # test --------------------------------------------------------------------

  # check # of basal species
  expect_equal(length(which(colSums(A) == 0)), b)
  expect_equal(length(which(colSums(alpha > 0) == 0)), b)

  # check # of consumer species
  expect_equal(length(which(colSums(A) != 0)), s - b)
  expect_equal(length(which(colSums(alpha > 0) > 0)), s - b)

  # compare expected and realized number of links
  expect_equal(l, mu_l, tolerance = 0.1)
})

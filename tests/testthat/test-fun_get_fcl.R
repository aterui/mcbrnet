
test_that("get_fcl", {

  # setup -------------------------------------------------------------------

  n_patch <- round(runif(1, 8, 50))

  # sample matrix
  m <- matrix(rpois(n_patch * 3, 2),
              nrow = 3,
              ncol = n_patch)
  m[, sample(n_patch, 3)] <- 0

  ## all possible combinations
  m[, 1] <- c(1, 0, 0)
  m[, 2] <- c(0, 1, 0)
  m[, 3] <- c(0, 0, 1)
  m[, 4] <- c(1, 1, 0)
  m[, 5] <- c(1, 0, 1)
  m[, 6] <- c(0, 1, 1)
  m[, 7] <- c(1, 1, 1)

  # preference
  delta <- runif(n_patch)

  # fcl
  fcl <- fun_get_fcl(x = m, delta = delta)
  fcl0 <- fcl[which(m[1, ] == 0)]
  fcl1 <- fcl[which(m[1, ] > 0 & m[2, ] == 0 & m[3, ] == 0)]
  fcl2 <- fcl[which(m[1, ] > 0 & m[2, ] > 0 & m[3, ] == 0)]
  fcl3 <- fcl[which(m[1, ] > 0 & m[2, ] > 0 & m[3, ] > 0)]

  # test --------------------------------------------------------------------

  expect_equal(unique(fcl0), 0)
  expect_equal(unique(fcl1), 1)
  expect_equal(unique(fcl2), 2)
  expect_true(all(fcl3 >= 2))
})

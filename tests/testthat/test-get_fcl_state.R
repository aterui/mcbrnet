test_that("get_fcl_state returns length(vector) = ncol(N)", {
  n_sp = 3
  n_patch = sample(3:10, 1)
  N <- matrix(rpois(n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
  N[,N[1,] == 0] <- 0
  expect_length(get_fcl_state(N = N), ncol(N))
})

test_that("get_fcl_state returns value from 0-3", {
  # matrix with FCL = 0, 1, 2, 2.5, and 3
  N <- matrix(c(0, 0, 0,
                1, 0, 0,
                1, 1, 0,
                1, 0, 1,
                1, 1, 1), nrow = 3, ncol = 5)

  out <- get_fcl_state(N = N)
  expect_identical(out, c(0.0, 1.0, 2.0, 2.5, 3.0))
})

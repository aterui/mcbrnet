
# setup -------------------------------------------------------------------

n_patch <- round(runif(1, 5, 15))

# sample matrix
m <- matrix(rpois(n_patch * 3, 2),
            nrow = 3,
            ncol = n_patch)
m[, sample(n_patch, 3)] <- 0

# preference
delta <- runif(n_patch)

# fcl
fcl <- fun_get_fcl(x = m, delta = delta)

# expected minimum
exp_fcl <- ifelse(m[2, ] > 0 | m[3, ] > 0, 2, 0) +
           ifelse(m[2, ] == 0 & m[3, ] == 0 & m[1, ] > 0, 1, 0)


# test --------------------------------------------------------------------

test_that("get fcl", {
  expect_equal(sum(fcl >= exp_fcl),
               n_patch)
})

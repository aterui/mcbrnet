
test_that("check adjtodist output", {

  # setup -------------------------------------------------------------------

  n_patch <- 10
  net <- brnet(n_patch = n_patch,
               p_branch = 0)

  x <- net$adjacency_matrix

  # test --------------------------------------------------------------------

  expect_equal(max(adjtodist(x)), n_patch - 1)

})

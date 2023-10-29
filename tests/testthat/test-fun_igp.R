
test_that("fun_igp", {

  # setup -------------------------------------------------------------------

  ## number of patch
  n_patch <- round(runif(1, 5, 15))

  ## sample matrices
  m_zero <- matrix(0,
                   nrow = 3,
                   ncol = n_patch)

  m_100 <- matrix(100,
                  nrow = 3,
                  ncol = n_patch)

  m_noc <- m_100
  m_noc[2, ] <- 0

  ## sample output
  y_zero <- fun_igp(m_zero)
  y_e_zero <- fun_igp(m_100, e = c(0, 0, 0))
  y_noc <- fun_igp(m_noc, h = c(0, 0, 0))
  y_100 <- fun_igp(m_100, h = c(0, 0, 0))

  ## reference
  zero_ref <- matrix(0,
                     nrow = 3,
                     ncol = n_patch)

  delta_ref <- rep(1, n_patch)

  # test --------------------------------------------------------------------

  expect_equal(y_zero$m_n_hat,
               zero_ref)

  expect_equal(y_noc$delta,
               delta_ref)

})

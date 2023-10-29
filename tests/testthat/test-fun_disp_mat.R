
test_that("fun_disp_mat", {

  # setup -------------------------------------------------------------------

  # matrix
  n_patch <- round(runif(1, 5, 15))
  m_distance <- dist(cbind(runif(n_patch, 0, 10),
                           runif(n_patch, 0, 10)),
                     diag = TRUE,
                     upper = TRUE)
  m_distance <- data.matrix(m_distance)

  m_dispersal <- exp(-m_distance)
  diag(m_dispersal) <- 0

  # sample
  m1 <- fun_disp_mat(n_patch = n_patch,
                     theta = 1,
                     landscape_size = 10)

  m2 <- fun_disp_mat(n_patch = n_patch,
                     theta = 1,
                     distance_matrix = m_distance)

  m3 <- fun_disp_mat(n_patch = n_patch,
                     theta = 1,
                     distance_matrix = m_distance,
                     dispersal_matrix = m_dispersal)

  # test --------------------------------------------------------------------

  # fun_disp_mat: ifelse structure
  expect_error(fun_disp_mat(n_patch = n_patch,
                            theta = 1))

  # fun_disp_mat: distance/dispersal matrix
  ## only distance matrix provided
  expect_true(is.matrix(m2$m_distance))

  ## distance & dispersal matrices provided
  expect_true(is.null(m3$m_distance))
  expect_true(is.matrix(m3$m_dispersal))

  ## m1 & m2 outputs equivalency
  expect_equal(m2$m_distance, m_distance)
  expect_equal(m3$m_dispersal, m_dispersal)

  ## fun_disp_mat: df_xy_coord
  expect_true(is.null(m2$df_xy_coord))
  expect_true(is.null(m3$df_xy_coord))

})

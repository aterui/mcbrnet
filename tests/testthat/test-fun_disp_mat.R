
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
                   distance_matrix = m_distance,
                   dispersal_matrix = m_dispersal)

m2 <- fun_disp_mat(n_patch = n_patch,
                   theta = 1,
                   landscape_size = 10)


# test --------------------------------------------------------------------

test_that("check ifelse structure", {
  expect_error(fun_disp_mat(n_patch = n_patch,
                            theta = 1,
                            dispersal_matrix = m_dispersal))
  expect_error(fun_disp_mat(n_patch = n_patch,
                            theta = 1))
})

test_that("distance/dispersal matrix", {
  expect_equal(m1$m_distance,
               m_distance)
  expect_equal(m1$m_dispersal,
               m_dispersal)
  expect_true(is.null(m1$df_xy_coord))
})

test_that("df_xy_coord", {
  expect_false(is.null(m2$df_xy_coord))
})


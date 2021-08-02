test_that("dispersal_n output dim", {

  n_patch = sample(3:10, 1)
  n_sp = 3
  v_theta = c(1, 1, 1)

  x_coord = runif(n_patch, 0, 5)
  y_coord = runif(n_patch, 0, 5)
  dist_mat = data.matrix(dist(cbind(x_coord,
                                    y_coord)))

  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
  diag(m_b_dispersal) <- 0
  # species C
  m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
  diag(m_c_dispersal) <- 0
  # species P
  m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
  diag(m_p_dispersal) <- 0

  N <- matrix(rpois(n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
  N[,N[1,] == 0] <- 0

  out <- dispersal_n(N = N, dist_mat = dist_mat,
                     m_b_dispersal = m_b_dispersal,
                     m_c_dispersal = m_c_dispersal,
                     m_p_dispersal = m_p_dispersal,
                     v_theta = v_theta,
                     v_p_dispersal = c(0.1, 0.1, 0.1))

  expect_length(out, length(N))
})

test_that("dispersal_n output not NA", {
  n_patch = sample(3:10, 1)
  n_sp = 3
  v_theta = c(1, 1, 1)

  x_coord = runif(n_patch, 0, 5)
  y_coord = runif(n_patch, 0, 5)
  dist_mat = data.matrix(dist(cbind(x_coord,
                                    y_coord)))

  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
  diag(m_b_dispersal) <- 0
  # species C
  m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
  diag(m_c_dispersal) <- 0
  # species P
  m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
  diag(m_p_dispersal) <- 0

  N <- matrix(rpois(n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
  N[,N[1,] == 0] <- 0

  out <- dispersal_n(N = N, dist_mat = dist_mat,
                     m_b_dispersal = m_b_dispersal,
                     m_c_dispersal = m_c_dispersal,
                     m_p_dispersal = m_p_dispersal,
                     v_theta = v_theta,
                     v_p_dispersal = c(0.1, 0.1, 0.1))

  out_na <- all(!is.na(out))
  out_gte_0 <- all(out >=0)
  expect_equal(out_na, TRUE)
  expect_equal(out_gte_0, TRUE)
})




# setup -------------------------------------------------------------------

n_species <- round(runif(1, 5, 15))
n_patch <- round(runif(1, 40, 60))
alpha <- runif(1, 0.5, 1)

m_n_hat <- matrix(rpois(n = n_species * n_patch,
                        lambda = 10),
                  nrow = n_species,
                  ncol = n_patch)

v_p_dispersal <- runif(n_species, 0, 1)

net <- brnet(n_patch = n_patch,
             p_branch = 0.5)

m_dispersal <- exp(-alpha * net$distance_matrix)
diag(m_dispersal) <- 0


m_n_prime <- fun_dispersal(x = m_n_hat,
                           v_p_dispersal = v_p_dispersal,
                           m_dispersal = m_dispersal)


# test --------------------------------------------------------------------

test_that("dispersal process", {
  expect_equal(rowSums(m_n_hat),
               rowSums(m_n_prime))
})

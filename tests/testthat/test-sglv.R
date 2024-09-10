
test_that("sglv", {

  # parameter ---------------------------------------------------------------

  ## set food web parameters
  ## s, number of species
  ## b, number of basal species
  ## l, expected number of links
  ## patch, number of patches
  ## theta, trophic specialization param
  s <- round(runif(1, min = 15, max = 20))
  b <- round(runif(1, min = 1, max = round(s / 2)))
  l <- round(runif(1, min = s - b, max = sum(b:(s - 1))))
  patch <- round(runif(1, min = 1, max = 10))
  theta <- runif(1, min = 0, max = 10)
  A <- ppm(s, b, l, theta)
  alpha <- to_alpha(A)

  ## find r values
  R <- findr(alpha, k0 = 10, lambda0 = 0.01)

  ## format
  r <- fun_to_m(R[, 1],
                n_species = s,
                n_patch = patch,
                param_attr = "species")$m_x

  ## for comparison
  x <- fun_to_m(R[, 2],
                n_species = s,
                n_patch = patch,
                param_attr = "species")$m_x


  net <- brnet(patch)

  # run simulation ----------------------------------------------------------

  ## R version
  set.seed(123)
  cout1 <- sglv(n_species = s,
                n_patch = patch,
                n_timestep = 200,
                r = r,
                dispersal = list(adj = net$adjacency_matrix,
                                 phi = runif(1, 0, 0.1),
                                 m = runif(1, 0, 0.001)),
                alpha = alpha,
                threshold = runif(1, 0, 0.001),
                cpp = FALSE)

  ## C++ version
  set.seed(123)
  cout2 <- sglv(n_species = s,
                n_patch = patch,
                n_timestep = 200,
                r = r,
                dispersal = list(adj = net$adjacency_matrix,
                                 phi = runif(1, 0, 0.1),
                                 m = runif(1, 0, 0.001)),
                alpha = alpha,
                threshold = runif(1, 0, 0.001),
                cpp = TRUE)

  x_prime <- c(x)
  x0_r <- colMeans(cout1[((nrow(cout1) - 1) / 2):nrow(cout1), -1])
  x0_cpp <- colMeans(cout2[((nrow(cout2) - 1) / 2):nrow(cout2), -1])
  rho1 <- cor(x_prime, x0_r)
  rho2 <- cor(x_prime, x0_cpp)

  expect_equal(x0_r, x0_cpp, tolerance = 1E-5)
  expect_gt(rho1, 0.99)
  expect_gt(rho2, 0.99)
})

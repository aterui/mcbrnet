
# setup -------------------------------------------------------------------

s <- round(runif(1, min = 15, max = 20))
b <- round(runif(1, min = 1, max = round(s / 2)))
l <- round(runif(1, min = s - b, max = sum(b:(s - 1))))
patch <- round(runif(1, min = 1, max = 10))
theta <- runif(1, min = 0, max = 10)
A <- ppm(s, b, l, theta)
alpha <- to_alpha(A)

R <- findr(alpha, k0 = 10, lambda0 = 0.01)
r <- fun_to_m(R[, 1],
              n_species = s,
              n_patch = patch,
              param_attr = "species")$m_x
x <- t(fun_to_m(R[, 2],
                n_species = s,
                n_patch = patch,
                param_attr = "species")$m_x)

cout <- sglv(n_species = s,
             n_patch = patch,
             n_timestep = 100,
             r = r,
             alpha = alpha,
             disturb = list(int = 0,
                            rate = 1,
                            s = 1),
             threshold = 0,
             cpp = FALSE)

x_prime <- round(c(x), 3)
x0 <- round(cout[nrow(cout), -1], 3)
rho <- cor(x_prime, x0)

# test --------------------------------------------------------------------

test_that("sglv", {
  expect_gt(rho, 0.99)
})


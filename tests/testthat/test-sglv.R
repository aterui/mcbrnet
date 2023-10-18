
# setup -------------------------------------------------------------------

s <- round(runif(1, min = 5, max = 20))
b <- round(runif(1, min = 1, max = s - 1))
l <- round(runif(1, min = s - b, max = sum(b:(s - 1))))
patch <- round(runif(1, min = 1, max = 10))
theta <- runif(1, min = 0, max = 10)
A <- ppm(s, b, l, theta)
alpha <- to_alpha(A)

# equilibrium density
x <- matrix(runif(s * patch), nrow = patch, ncol = s)

# intrinsic growth
r <- -x %*% alpha

cout <- sglv(n_species = s,
             n_patch = patch,
             n_timestep = 100,
             r = r,
             alpha = alpha,
             phi = 0,
             connectivity = matrix(0, patch, patch),
             disturb = list(int = 0,
                            rate = 1,
                            s = 1),
             threshold = 0)

x_prime <- round(c(x), 3)
x0 <- round(cout[nrow(cout), -1], 3)
rho <- cor(x_prime, x0)


# test --------------------------------------------------------------------

test_that("sglv", {
  expect_gt(rho, 0.99)
})


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

## sample output
y_zero <- fun_igp(m_zero)
y_e_zero <- fun_igp(m_100, e = c(0, 0, 0))
y_100 <- fun_igp(m_100, h = c(0, 0, 0))

## reference
zero_ref <- matrix(0,
                   nrow = 3,
                   ncol = n_patch)

z_ref <- matrix(0.5,
                nrow = 3,
                ncol = n_patch)
z_ref[2:3, ] <- 0.25

# test --------------------------------------------------------------------

test_that("N-zero matrix", {
  expect_equal(y_zero$z,
               zero_ref)
  expect_equal(y_zero$m_n_hat,
               zero_ref)
})

test_that("z values", {
  expect_equal(y_100$z,
               z_ref)
})


# setup -------------------------------------------------------------------

nsp <- round(runif(1, 5, 10))
np <- round(runif(1, 5, 10))

x1 <- runif(nsp)
m_sp <- fun_to_m(x = x1,
                 n_patch = np,
                 n_species = nsp,
                 param_attr = "species")

x2 <- runif(np)
m_p <- fun_to_m(x = x2,
                n_patch = np,
                n_species = nsp,
                param_attr = "patch")


# test --------------------------------------------------------------------

test_that("duplication works", {
  expect_equal(apply(m_sp$m_x, 1, unique), x1)
  expect_equal(apply(m_p$m_x, 2, unique), x2)
})

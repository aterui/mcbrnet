
# setup -------------------------------------------------------------------

n_patch <- round(runif(1, min = 5, max = 50))
p_branch <- 0
n_source <- round(runif(1, min = 1, max = n_patch))

x <- brnet(n_patch = n_patch,
           p_branch = p_branch)

y <- ptsource(x,
              n_source = n_source,
              p = 1,
              q = 0)

y <- y$df_patch

# test --------------------------------------------------------------------

test_that("accumulation", {
  expect_equal(max(dplyr::pull(y, impact)), n_source)
})

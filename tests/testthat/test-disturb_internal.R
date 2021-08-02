test_that("disturb_type = NULL, no disturbance", {
  N <- matrix(rpois(3*5, lambda = 100),
              nrow = 3, ncol = 5)
  out <- disturb_internal(N = N, disturb_type = NULL)
  expect_equal(out$N, N)
  expect_equal(length(out$patch_extinction), ncol(N))
})


test_that("extreme disturbance does not result in -N", {
n_patch = 5
N <- matrix(c(50, 50, 50, #patch 1, at bottom
              100, 100, 100, #next three patches are source
              100, 50, 10,
              50, 50, 50,
              100, 100, 100), # patch 5, confluence
            ncol = n_patch)

dist_mat <- matrix(
  c(0, 1, 2, 2, 1,
    1, 0, 3, 3, 2,
    2, 3, 0, 2, 1,
    2, 3, 2, 0, 1,
    1, 2, 1, 1, 0),
  ncol = n_patch,
  dimnames = list(
    c("patch1", "patch2", "patch3", "patch4", "patch5"),
    c("patch1", "patch2", "patch3", "patch4", "patch5")))
disturb_value <- boot::inv.logit(rnorm(n = n_patch,
                                       mean = boot::logit(0.999),
                                       sd = 0.0001))

out <- disturb_internal(N = N,
                        disturb_type = "regional",
                        disturb_value = disturb_value,
                        disturb_p = 1,
                        river_network_structure = FALSE)
expect_true(all(out$N > matrix(0, nrow = 3, ncol = n_patch)))
expect_equal(length(out$patch_extinction), ncol(N))}
)

test_that("regional disturbance reduces N correctly when env_val known",{
  n_patch = 5
  N <- matrix(100, # patch 5, confluence
              ncol = n_patch,
              nrow = 3)
  disturb_value = boot::inv.logit(c(-2, -1, 0, 1, 2))

out <- disturb_internal(N = N,
                        disturb_value = disturb_value,
                        disturb_type = "regional",
                        disturb_p = 1,
                        river_network_structure = TRUE)
# output is reduced by proper percentages (1 - env_scale)
expect_equal(round(colSums(out$N) / colSums(N), 3),
             round(1 - disturb_value, 3))
expect_equal(length(out$patch_extinction), ncol(N))
}
)

test_that("regional disturbance 2d; \n output < input \nvariation in output col \nncol output = ncol input", {

  n_patch = 10
  N <- matrix(100,
              ncol = n_patch,
              nrow = 3)
  disturb_value <- boot::inv.logit(rnorm(n_patch,
                         mean = boot::logit(0.99),
                         sd = 0.1))

  out <- disturb_internal(N = N,
                          disturb_type = "regional",
                          disturb_p = 1,
                          river_network_structure = FALSE,
                          disturb_value = disturb_value)
  expect_lte(sum(out$N), sum(N))
  expect_false(all(colSums(out$N)==colSums(out$N)[1]))
  expect_equal(length(out$patch_extinction), ncol(N))
          }
)

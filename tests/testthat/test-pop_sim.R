test_that("pop_sim out dim and type", {
  N <- matrix(rpois(3*5, lambda = 100),
              nrow = 3, ncol = 5)
  out <- pop_sim(N = N,
                 P_pref = 0.25, fixed_P_pref = TRUE,
                 alphabc = 4, betabc = 20, ebc = 2,
                 alphap = 4, betap = 20, ebp = 2, ecp = 2,
                 v_s0 = c(0.75, 0.75, 0.75),
                 k = 150, r_max = 2.5, b = 0.01)$N
  expect_equal(dim(out), dim(N))
  expect_type(out, "double")
  expect_equal(is.matrix(out), TRUE)
})

test_that("pop_sim ignore 'P_pref' when 'fixed_P_pref = FALSE",{
  N <- matrix(rpois(3*5, lambda = 100),
              nrow = 3, ncol = 5)
  out1 <- pop_sim(N = N,
                 P_pref = 0.25, fixed_P_pref = FALSE,
                 alphabc = 4, betabc = 20, ebc = 2,
                 alphap = 4, betap = 20, ebp = 2, ecp = 2,
                 v_s0 = c(0.75, 0.75, 0.75),
                 k = 150, r_max = 2.5, b = 0.01)$N
  out2 <- pop_sim(N = N,
                  P_pref = NULL, fixed_P_pref = FALSE,
                  alphabc = 4, betabc = 20, ebc = 2,
                  alphap = 4, betap = 20, ebp = 2, ecp = 2,
                  v_s0 = c(0.75, 0.75, 0.75),
                  k = 150, r_max = 2.5, b = 0.01)$N
  expect_equal(out1, out2)
})

test_that("pop_sim fails when fixed_P_pref = TRUE but no value for P_pref supplied", {
  expect_error(pop_sim(N = N,
          P_pref = NULL, fixed_P_pref = TRUE,
          alphabc = 4, betabc = 20, ebc = 2,
          alphap = 4, betap = 20, ebp = 2, ecp = 2,
          v_s0 = c(0.75, 0.75, 0.75),
          k = 150, r_max = 2.5, b = 0.01),
          regexp = "fixed_P_pref == TRUE, but no value for `P_pref` supplied")

})

test_that("k_n_upstream returns vectors of length n_upstream", {
  out1 <- k_n_upstream(r_max = 5, n_upstream = 5, n_patch = 1)
  out5 <- k_n_upstream(r_max = 5, n_upstream = c(1,2,3,4,5), n_patch = 5)

  expect_length(out1$b, 1)
  expect_length(out5$b, 5)
  expect_type(out1$k, "double")
  expect_type(out1$b, "double")
})

test_that("k_n_upstream k >= k_base + k_c", {
  k_base = 150
  k_c = 10

  expect_gte(k_n_upstream(k_min_exponent = 1.1, k_max_exponent = 1.1,
    r_max = 5, n_upstream = 1, n_patch = 1)$k,
            k_base + k_c)
})

test_that("k_n_upstream returns list length 2", {
  out <- k_n_upstream(r_max = 5, n_upstream = 1, n_patch = 1)
  expect_length(out, 2)
  expect_type(out, "list")
})

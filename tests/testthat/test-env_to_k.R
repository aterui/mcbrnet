test_that("env_to_k returns positive", {
  out <- env_to_k(c(-100, 0, 100), k_base = 0.001)
  expect_identical(all(out >= 0), TRUE)
})


test_that("env_to_k returns k_base when length(env) == 1", {
  k_base = 150
  out <- env_to_k(0, k_base = k_base)
  expect_identical(out, k_base)
})

test_that("env_to_k returns [0.5, 1.5] when variation is high", {
  k_base = 100
  env = rnorm(1000, sd = 100)
  out <- env_to_k(env = env, k_base = k_base)
  min_out <- min(out)
  max_out <- max(out)
  expect_gte(min_out, 50)
  expect_lte(min_out, 150)
})

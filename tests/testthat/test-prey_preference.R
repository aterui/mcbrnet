test_that("prey_preference returns vaue from 0-1", {
  out <- prey_preference(1, 1, 1, 1)
  expect_gte(out, 0)
  expect_lte(out, 1)
})

#----------------------------------------------

test_that("prey_preference stops when NA input", {
  expect_error(prey_preference(1, 1, 1, NA))
})

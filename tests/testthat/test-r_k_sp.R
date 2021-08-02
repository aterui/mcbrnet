test_that("r_k_sp; k and r_max output >=0", {
  expect_identical(all(r_k_sp()$k >= 0), TRUE)
  expect_identical(all(r_k_sp()$r_max >= 0), TRUE)
})

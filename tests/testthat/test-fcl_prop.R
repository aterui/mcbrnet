test_that("fcl_prop output names", {
  names.out <- names(
    fcl_prop(
      rep(c(0, 1, 2, 2.5, 3),
          each = 3)))
  expect_identical(names.out,
                   c("No_spp", "B_only", "B_C", "B_P", "B_C_P"))
})

test_that("fcl_prop output values", {
  out <- fcl_prop(rep(c(0, 1, 2, 2.5, 3), each = 3))
  expect_equal(out, c(0.2, 0.2, 0.2, 0.2, 0.2),
               ignore_attr = TRUE)
})

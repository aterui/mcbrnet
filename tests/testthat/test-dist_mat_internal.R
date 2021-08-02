test_that("dist_mat_internal input checks", {
  expect_message(
    dist_mat_internal(
      dist_mat = NULL,
      landscape_size = NULL,
      n_patch = 10),
    regexp = "No distance matrix supplied")

   expect_error(
    dist_mat_internal(
      dist_mat = c(1:10),
      landscape_size = NULL,
      n_patch = 10),
    regexp = "distance matrix should be provided as matrix")

   expect_error(
    dist_mat_internal(
      dist_mat = matrix(0, ncol = 9, nrow = 9),
      landscape_size = NULL,
      n_patch = 10),
    regexp = "invalid dimension:")

   expect_error(
     dist_mat_internal(
       dist_mat = matrix(1, ncol = 10, nrow = 10),
       landscape_size = NULL,
       n_patch = 10),
     regexp = "invalid distance matrix: diagonal")
})

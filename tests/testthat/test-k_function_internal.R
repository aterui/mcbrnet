test_that("k_function_internal 'NULL' returns k_base", {
  k_base = 150
  n_patch = 20
  out <- suppressMessages(k_function_internal(
    k_function = NULL,
                      k_base = k_base,
                      r_max = 2.5,
                      n_upstream = NULL,
                      n_patch = n_patch,
                      river_network_structure = TRUE,
                      environment_value = NULL))
  expect_equal(out$k, k_base,
               ignore_attr = TRUE)
  expect_named(out, c("k", "b"))
  expect_length(out$k, length(k_base))
})

# dist_mat_internal(dist_mat = NULL,
#                   landscape_size = 10,
#                   n_patch = 20)
#
# k_function_internal(k_function = NULL,
#                     k_base = 150,
#                     r_max = 2.5,
#                     n_upstream = NULL,
#                     n_patch = 20,
#                     river_network_structure = TRUE,
#                     environment_value = NULL)

test_that("k_function_internal patches upstream, n_patch needs to be supplied", {
  expect_error(k_function_internal(k_function = "patches-upstream",
                                     k_base = 150,
                                     r_max = 2.5,
                                     n_upstream = NULL,
                                     n_patch = 20,
                                     river_network_structure = TRUE,
                                     environment_value = NULL),
                 regexp = "number of patches upstream must be supplied to n_upstream argument")
  })

test_that("k_function_internal 'patches upstream' returns list and length = n_patch", {
  n_patch = 20
  out <- suppressMessages(
    k_function_internal(
      k_function = "patches-upstream",
      k_base = 150,
      r_max = 2.5,
      n_upstream = sample(1:20, 20),
      n_patch = n_patch,
      k_c = 10,
      k_min_exponent = 1.1,
      k_max_exponent = 1.35,
      river_network_structure = TRUE,
      environment_value = NULL))
  expect_type(out, "list")
  expect_named(out, c("k", "b"))
  expect_length(out$k, n_patch)
})


test_that("k_function_internal 'environment' returns list length n_patch", {
  n_patch = 20
  k_base = 150
  environment_value = rnorm(n_patch)
  out <- suppressMessages(
    k_function_internal(k_function = "environment",
                        k_base = k_base,
                        r_max = 2.5,
                        #n_upstream = sample(1:20, 20),
                        n_patch = n_patch,
                        #k_c = 10,
                        #k_min_exponent = 1.1,
                        #k_max_exponent = 1.35,
                        #river_network_structure = TRUE,
                        environment_value = environment_value))
  expect_type(out, "list")
  expect_named(out, c("k", "b"))
  expect_length(out$k, n_patch)
})

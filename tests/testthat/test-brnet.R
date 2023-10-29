
test_that("brnet", {

  # setup -------------------------------------------------------------------

  n_patch <- round(runif(1, 10, 50))
  mu_disturb <- runif(1, 0, 1)
  mu_env <- rnorm(1)

  net1 <- brnet(n_patch = n_patch,
                p_branch = 0.5,
                sd_env_lon = 0,
                sd_env_source = 0,
                sd_disturb_lon = 0,
                mean_env_source = mu_env,
                mean_disturb_source = mu_disturb,
                sd_disturb_source = .1)

  net2 <- brnet(n_patch = n_patch,
                p_branch = 0.5,
                sd_env_lon = 0,
                sd_env_source = 0,
                sd_disturb_lon = 0,
                mean_env_source = mu_env,
                mean_disturb_source = mu_disturb,
                sd_disturb_source = 0)

  linear <- brnet(n_patch = n_patch,
                  p_branch = 0,
                  mean_disturb_source = mu_disturb,
                  sd_disturb_lon = 0,
                  sd_disturb_source = 0.1)


  # test --------------------------------------------------------------------

  # basic output
  expect_error(brnet(n_patch = 100, p_branch = 1))
  expect_equal(net1$df_patch$n_patch_upstream[1], n_patch)
  expect_equal(max(linear$distance_matrix), n_patch - 1)

  # local environment
  expect_equal(unique(round(net1$df_patch$environment, 3)), round(mu_env, 3))

  # disturbance
  expect_equal(dplyr::n_distinct(net1$df_patch$disturbance),
               dplyr::n_distinct(net1$df_patch$branch_id))
  expect_equal(unique(round(net2$df_patch$disturbance, 3)), round(mu_disturb, 3))
  expect_equal(dplyr::n_distinct(linear$df_patch$disturbance), 1)
})

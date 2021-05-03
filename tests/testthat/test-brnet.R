
# setup -------------------------------------------------------------------

  context("brnet")

  n_patch <- 10
  net <- brnet(n_patch = n_patch,
               p_branch = 0.5,
               sd_env_lon = 0,
               sd_env_source = 0,
               mean_env_source = 10)

  linear <- brnet(n_patch = n_patch,
                  p_branch = 0)

# test --------------------------------------------------------------------

  test_that("check brnet output", {
    expect_error(brnet(n_patch = 100, p_branch = 1))

    expect_equal(net$df_patch$n_patch_upstream[1], n_patch)
    expect_equal(unique(round(net$df_patch$environment, 2)), 10)
    expect_equal(max(linear$distance_matrix), n_patch - 1)
  })

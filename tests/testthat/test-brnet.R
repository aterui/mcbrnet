
# setup -------------------------------------------------------------------

  context("brnet")

  n_patch <- 10
  net <- brnet(n_patch = n_patch,
               p_branch = 0.5,
               sd_env_lon = 0,
               sd_env_source = 0,
               mean_env_source = 0)

  linear <- brnet(n_patch = n_patch,
                  p_branch = 0)

# test --------------------------------------------------------------------

  test_that("check environmental/niche parameters", {
    expect_equal(net$df_patch$n_patch_upstream[1], n_patch)
    expect_equal(unique(net$df_patch$environment), 0)
    expect_equal(max(linear$distance_matrix), n_patch - 1)
  })


# setup -------------------------------------------------------------------

context("brnet")

n_patch <- 10
net1 <- brnet(n_patch = n_patch,
              p_branch = 0.5,
              sd_env_lon = 0,
              sd_env_source = 0,
              mean_env_source = 10,
              mean_disturb_source = 0.9,
              sd_disturb_source = .1,
              plot = FALSE)

net2 <- brnet(n_patch = n_patch,
              p_branch = 0.5,
              sd_env_lon = 0,
              sd_env_source = 0,
              mean_env_source = 10,
              mean_disturb_source = 0.9,
              sd_disturb_source = 0,
              plot = FALSE)

linear <- brnet(n_patch = n_patch,
                p_branch = 0,
                mean_disturb_source = 0.9,
                sd_disturb_source = 0.1,
                plot = FALSE)

# test --------------------------------------------------------------------

test_that("check brnet output", {
  expect_error(brnet(n_patch = 100, p_branch = 1))

  expect_equal(dplyr::n_distinct(net1$df_patch$disturbance),
               dplyr::n_distinct(net1$df_patch$branch_id))
  expect_equal(unique(round(net2$df_patch$disturbance, 1)), 0.9)
  expect_equal(dplyr::n_distinct(linear$df_patch$disturbance), 1)
  expect_equal(net1$df_patch$n_patch_upstream[1], n_patch)
  expect_equal(unique(round(net1$df_patch$environment, 2)), 10)
  expect_equal(max(linear$distance_matrix), n_patch - 1)
})

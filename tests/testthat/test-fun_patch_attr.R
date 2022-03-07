
# setup -------------------------------------------------------------------

x <- round(runif(1, 50, 150))
y <- runif(1, 0.2, 0.8)

net <- brnet(n_patch = x,
             p_branch = y,
             randomize_patch = FALSE,
             plot = F)

m_adj <- net$adjacency_matrix
branch <- dplyr::n_distinct(net$df_patch$branch_id)
m_dist <- net$distance_matrix
v_wa <- net$df_patch$n_patch_upstream

## patch environment
mu_env <- runif(1)
v_env0 <- fun_patch_attr(x = m_adj,
                         n_branch = branch,
                         mean_source = mu_env,
                         sd_source = 0,
                         sd_lon = 0,
                         m_distance = m_dist,
                         rho = 1,
                         v_wa = v_wa)

v_env1 <- fun_patch_attr(x = m_adj,
                         n_branch = branch,
                         mean_source = mu_env,
                         sd_source = 1,
                         sd_lon = 0,
                         m_distance = m_dist,
                         rho = 1,
                         v_wa = v_wa)


## patch disturbance
mu_disturb <- runif(1, 0, 1)
v_disturb0 <- fun_patch_attr(x = m_adj,
                             n_branch = branch,
                             mean_source = mu_disturb,
                             sd_source = 0,
                             sd_lon = 0,
                             m_distance = m_dist,
                             rho = 1,
                             v_wa = v_wa)

v_disturb1 <- fun_patch_attr(x = m_adj,
                             n_branch = branch,
                             mean_source = mu_disturb,
                             sd_source = 1,
                             sd_lon = 0,
                             m_distance = m_dist,
                             rho = 1,
                             v_wa = v_wa)


# test --------------------------------------------------------------------

test_that("patch environment no variation", {
  expect_equal(unique(round(v_env0, 4)), round(mu_env, 4))
})

test_that("patch disturbance no variation", {
  expect_equal(unique(round(v_disturb0, 4)), round(mu_disturb, 4))
})

test_that("patch environment source variation", {
  expect_equal(dplyr::n_distinct(v_env1), branch)
})

test_that("patch disturbance source variation", {
  expect_equal(dplyr::n_distinct(v_disturb1), branch)
})

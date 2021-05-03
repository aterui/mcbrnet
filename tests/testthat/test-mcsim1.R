
# setup -------------------------------------------------------------------

  context("mcsim")

  n_patch <- 50
  n_species <- 10
  mean_env <- rnorm(n_patch)
  niche_optim <- rnorm(n_species)
  sd_niche_width <- 1:n_species * 0.5
  K <- 1:n_patch

  output <- mcsim(n_species = n_species,
                  n_patch = n_patch,
                  mean_env = mean_env,
                  niche_optim = niche_optim,
                  sd_niche_width = sd_niche_width,
                  carrying_capacity = K)


# test --------------------------------------------------------------------

  test_that("check environmental/niche parameters", {
    expect_error(mcsim(n_species = n_species, niche_optim = c(0, 1)))
    expect_error(mcsim(n_species = n_species, sd_niche_width = c(0.1, 1)))

    expect_equal(unique(output$df_dynamics$mean_env), mean_env)
    expect_equal(unique(output$df_dynamics$niche_optim), niche_optim)
    expect_equal(unique(output$df_dynamics$carrying_capacity), K)

    expect_equal(output$df_species$niche_optim, niche_optim)
    expect_equal(output$df_species$sd_niche_width, sd_niche_width)
    expect_equal(output$df_patch$mean_env, mean_env)
    expect_equal(output$df_patch$carrying_capacity, K)
  })


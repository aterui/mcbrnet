
test_that("mean niche", {

  # setup -------------------------------------------------------------------

  n_patch <- round(runif(1, 40, 60))
  n_species <- round(runif(1, 5, 15))
  mean_env <- rnorm(n_patch)
  niche_optim <- rnorm(n_species)
  sd_niche_width <- 1:n_species * 0.5
  cc <- 10 * 1:n_patch
  distance_mat <- data.matrix(dist(cbind(runif(n_patch, 0, 10),
                                         runif(n_patch, 0, 10)),
                                   diag = TRUE,
                                   upper = TRUE))

  output1 <- mcsim(n_species = n_species,
                   n_patch = n_patch,
                   mean_env = mean_env,
                   niche_optim = niche_optim,
                   sd_niche_width = sd_niche_width,
                   carrying_capacity = cc)

  output2 <- mcsim(n_species = n_species,
                   n_patch = n_patch,
                   mean_env = mean_env,
                   niche_optim = niche_optim,
                   sd_niche_width = sd_niche_width,
                   carrying_capacity = cc,
                   distance_matrix = distance_mat,
                   interaction_type = "random",
                   min_alpha = 0,
                   max_alpha = 0.5)

  # test --------------------------------------------------------------------

  ## mean niche
  expect_error(mcsim(n_species = n_species,
                     niche_optim = rep(0, n_species + 1)))

  expect_equal(unique(output1$df_dynamics$niche_optim),
               niche_optim)

  expect_equal(output1$df_species$niche_optim,
               niche_optim)

  expect_equal(unique(output2$df_dynamics$niche_optim),
               niche_optim)

  expect_equal(output2$df_species$niche_optim,
               niche_optim)

  ## sd niche
  expect_error(mcsim(n_species = n_species, sd_niche_width =  rep(0, n_species + 1)))

  expect_equal(output1$df_species$sd_niche_width,
               sd_niche_width)

  expect_equal(output2$df_species$sd_niche_width,
               sd_niche_width)

  ## environment
  expect_equal(unique(output1$df_dynamics$mean_env),
               mean_env)

  expect_equal(output1$df_patch$mean_env,
               mean_env)

  expect_equal(unique(output2$df_dynamics$mean_env),
               mean_env)

  expect_equal(output2$df_patch$mean_env,
               mean_env)

  ## carrying capacity
  expect_error(mcsim(n_patch = n_patch, carrying_capacity =  rep(10, n_patch + 1)))

  expect_equal(unique(output1$df_dynamics$carrying_capacity),
               cc)

  expect_equal(output1$df_patch$carrying_capacity,
               cc)

  expect_equal(unique(output2$df_dynamics$carrying_capacity),
               cc)

  expect_equal(output2$df_patch$carrying_capacity,
               cc)
})

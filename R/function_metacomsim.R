
mcsim <- function(n_warmup = 200,
                  n_burnin = 200,
                  n_timestep = 1000,
                  propagule_interval = NULL,
                  n_species = 5,
                  n_patch = 5,
                  carrying_capacity = 100,
                  interaction_type = "constant",
                  alpha = 0,
                  min_alpha = NULL,
                  max_alpha = NULL,
                  r0 = 4,
                  sd_niche_width = 1,
                  optim_min = -1,
                  optim_max = 1,
                  m_distance = NULL,
                  landscape_size = 10,
                  mu_env = 0,
                  sd_env = 1,
                  phi = 1,
                  spatial_env_cor = FALSE,
                  p_dispersal = 0.1,
                  theta = 1,
                  plot = FALSE
                  ) {

# Required library --------------------------------------------------------
  library(mvtnorm)
  library(ggplot2)
  library(dplyr)

# parameter setup ---------------------------------------------------------

  # carrying capacity
  if (length(carrying_capacity) == 1) {
    print("Only one carrying capacity is given: the model will assume carrying capacities are the same at all habitat patches")
    m_k <- matrix(carrying_capacity, nrow = n_species, ncol = n_patch)
  } else {
    if(length(carrying_capacity) != n_patch) stop("carrying_capacity must have the length of one or n_patch")
    m_k <- matrix(rep(carrying_capacity, each = n_species), nrow = n_species, ncol = n_patch)
  }

  # species niche
  v_mu <- runif(n_species, optim_min, optim_max)
  m_mu <- matrix(rep(v_mu, n_patch), nrow = n_species, ncol = n_patch) # optimal value

  # environmental variation among patches
  if (length(mu_env) == 1) {
    print("Only one environmental value is given: the model will assume environmental conditions are the same at all habitat patches")
    v_mu_z <- rep(mu_env, n_patch)
  } else {
    if(n_patch != length(mu_env)) stop("mu_env must have the length of n_patch")
    v_mu_z <- mu_env
  }

  # interaction matrix
  if (interaction_type == "constant") {
    if (alpha < 0 | length(alpha) != 1) stop("invalid value of alpha - the value must be a positive scalar")
    alpha <- alpha
  }

  if (interaction_type == "random") {
    if (min_alpha < 0 | max_alpha < 0 | min_alpha > max_alpha | is.null(min_alpha) | is.null(max_alpha)) stop("invalid values of min_alpha and/or min_alpha - the values must be positive values")
    alpha <- runif(n_species * n_species, min_alpha, max_alpha)
  }

  m_interaction <- matrix(alpha, nrow = n_species, ncol = n_species)
  diag(m_interaction) <- 1

  # dispersal matrix
  if (is.null(m_distance)) {
    print("Distance matrix m_distance is not given: create a square landscape with landscape_size (default: 10) on a side")
    v_x_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    v_y_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    m_distance <- data.matrix(dist(cbind(v_x_coord, v_y_coord), diag = TRUE, upper = TRUE))
    m_dispersal <- data.matrix(exp(-theta*m_distance))
    diag(m_dispersal) <- 0
  } else {
    if (is.matrix(m_distance) == 0) stop("Distance matrix should be provided as matrix")
    if (nrow(m_distance) == n_patch) stop("Invalid dimension of distance matrix m_distance: m_distance must have a dimension of n_patch*n_patch")
    m_dispersal <- exp(-theta * m_distance)
    diag(m_dispersal) <- 0
  }


  # dynamics ----------------------------------------------------------------

  n_sim <- n_warmup + n_burnin + n_timestep
  n_discard <- n_warmup + n_burnin

  # environment
  var_env <- sd_env^2
  if (spatial_env_cor == TRUE) {
    m_sigma <- var_env * exp(-phi * m_distance)
  }else{
    m_sigma <- var_env * diag(x = 1, nrow = n_patch, ncol = n_patch)
  }
  m_z <- mvtnorm::rmvnorm(n_sim, mean = v_mu_z, sigma = m_sigma)

  # community
  colname <- c("N", "species", "niche_optim", "r_xt", "patch", "mean_env", "env", "timestep")
  m_dynamics <- matrix(NA, nrow = n_species * n_patch * n_timestep, ncol = length(colname))
  colnames(m_dynamics) <- colname
  st_row <- seq(1, nrow(m_dynamics), by = n_species * n_patch)

  if (is.null(propagule_interval)) propagule_interval <- max(c(1, floor(n_warmup / 10)))
  propagule <- seq(propagule_interval, max(c(1, n_warmup)), by = propagule_interval)

  m_n <- matrix(rpois(n_species * n_patch, 0.5), nrow = n_species, ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (n in 1:n_sim) {
    if (n_warmup > 0) {
      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n_species * n_patch, 0.5), nrow = n_species, ncol = n_patch)
      }
    }

    m_z_xt <- matrix(rep(m_z[n, ], each = n_species), nrow = n_species, ncol = n_patch)
    m_r_xt <- r0 * exp(-((m_mu - m_z_xt) / (sqrt(2) * sd_niche_width))^2)
    m_n_hat <- (m_n * m_r_xt)/(1 + ((r0 - 1) / m_k) * (m_interaction %*% m_n) )

    m_e_hat <- m_n_hat * p_dispersal
    v_e_sum <- rowSums(m_e_hat)
    m_i_raw <- m_e_hat %*% m_dispersal
    v_i_sum <- rowSums(m_i_raw)
    v_i_sum[v_i_sum == 0] <- 1
    m_i_prob <- m_i_raw / v_i_sum
    m_i_hat <- m_i_prob * v_e_sum
    m_n_hat <- m_n_hat + m_i_hat - m_e_hat

    m_n <- matrix(rpois(n = n_species * n_patch, lambda = m_n_hat), nrow = n_species, ncol = n_patch)

    if(n > n_discard){
      row_id <- st_row[n - n_discard]:(st_row[n - n_discard] + n_species * n_patch - 1)
      m_dynamics[row_id, ] <- cbind(c(m_n),
                                    rep(1:n_species, n_patch),
                                    c(m_mu),
                                    c(m_r_xt),
                                    rep(1:n_patch, each = n_species),
                                    rep(v_mu_z, each = n_species),
                                    c(m_z_xt),
                                    I(n - n_discard))
    }
    setTxtProgressBar(pb, n)
  }
  close(pb)
    if (plot == TRUE) {
      sample_patch <- sample(1:n_patch, size = min(c(n_patch, 5)), replace = FALSE)
      sample_species <- sample(1:n_species, size = min(c(n_species, 5)), replace = FALSE)

      g <- data.frame(m_dynamics) %>%
              dplyr::filter(patch %in% sample_patch, species %in% sample_species) %>%
              ggplot2::ggplot() +
              ggplot2::facet_grid(rows = vars(species), cols = vars(patch)) +
              ggplot2::geom_line(mapping = aes(x = timestep, y = m_n, color = abs(niche_optim - env))) +
              ggplot2::scale_color_viridis_c(alpha = 0.5) +
              ggplot2::labs(color = "Environmental \ndeviation")
      print(g)
    }

    return(list(dynamics = data.frame(m_dynamics),
                distance_matrix = m_distance))
}


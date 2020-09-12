#' Simulate metacommunity dynamics
#'
#' @param n_species numeric scalar. Number of species in a metacommunity.
#' @param n_patch numeric scalar. Number of patches in a metacommunity.
#' @param n_warmup numeric scalar. Number of time-steps for warm-up.
#' @param n_burnin numeric scalar. Number of time-steps for burn-in.
#' @param n_timestep numeric scalar. Number of time-steps to be saved.
#' @param propagule_interval numeric scalar. Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param carrying_capacity numeric scalar or vector (length should be equal to the number of patches). Carrying capacities of individual patches.
#' @param interaction_type character scalar. \code{"constant"} or \code{"random"}. \code{"constant"} assumes the single interaction strength of alpha for all pairs of species. \code{"random"} draws random numbers from a uniform distribution with \code{min_alpha} and \code{max_alpha}.
#' @param alpha species interaction strength.
#' @param min_alpha numeric scalar. Minimum value of a uniform distribution that generates alpha.
#' @param max_alpha numeric scalar. Maximum value of a uniform distribution that generates alpha.
#' @param r0 numeric scalar or vector (length should be equal to the number of species). Maximum population growth rate of the Beverton-Holt model.
#' @param sd_niche_width numeric scalar. Niche width of species. Higher values indicate greater niche width.
#' @param optim_min numeric scalar. Minimum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species.
#' @param optim_max numeric scalar. Maximum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species.
#' @param distance_matrix numeric matrix. Distance matrix indicating distance between habitat patches. If \code{NULL}, a square landscape with randomly distributed patches will be generated. Default \code{NULL}.
#' @param landscape_size numeric scalar. Length of a landscape on a side. Active only when \code{dispersal_matrix = NULL}.
#' @param mu_env numeric scalar or vector (with length equal to the number of patches). Mean environmental values of patches.
#' @param sd_env numeric scalar. Temporal SD of environmental variation at each patch.
#' @param phi numeric scalar. Parameter describing distance decay of spatial autocorrelation in temporal environmental variation (\eqn{\rho = exp(-\phi d)}, where \code{d} is the distance between patches).
#' @param spatial_env_cor logical. Indicates whether spatial autocorrelation in temporal environmental variation is considered or not. Default \code{FALSE}.
#' @param p_dispersal numeric scalar or vector (length should be equal to the number of species). Probability of dispersal (success probability of a binomial distribution).
#' @param theta numeric scalar. Dispersal parameter describing dispersal capability of species as \eqn{exp(-\theta d)}, where \code{d} is the distance between patches.
#' @param plot logical. If \code{TRUE}, results are plotted.
#'
#' @return \code{dynamics} temporal metacommunity dynamics.
#' @return \code{distance_matrix} distance matrix.
#' @return \code{interaction_matrix} species interaction matrix.
#'
#' @importFrom dplyr %>% filter
#' @importFrom ggplot2 ggplot vars labeller geom_line aes scale_color_viridis_c labs facet_grid label_both
#' @importFrom stats dist rpois
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' mcsim(n_warmup = 200, n_burnin = 200, n_timestep = 1000)
#'
#' @export
#'
mcsim <- function(n_species = 5,
                  n_patch = 5,
                  n_warmup = 200,
                  n_burnin = 200,
                  n_timestep = 1000,
                  propagule_interval = NULL,
                  carrying_capacity = 100,
                  interaction_type = "constant",
                  alpha = 0,
                  min_alpha = NULL,
                  max_alpha = NULL,
                  r0 = 4,
                  sd_niche_width = 1,
                  optim_min = -1,
                  optim_max = 1,
                  distance_matrix = NULL,
                  landscape_size = 10,
                  mu_env = 0,
                  sd_env = 0.1,
                  phi = 1,
                  spatial_env_cor = FALSE,
                  p_dispersal = 0.1,
                  theta = 1,
                  plot = FALSE
) {

  # parameter setup ---------------------------------------------------------

  # carrying capacity
  if (length(carrying_capacity) == 1) {
    print("Only one carrying capacity is given: the model will assume carrying capacities are the same at all habitat patches")
    m_k <- matrix(carrying_capacity, nrow = n_species, ncol = n_patch)
  } else {
    if (length(carrying_capacity) != n_patch) stop("carrying_capacity must have length of one or n_patch")
    m_k <- matrix(rep(x = carrying_capacity, each = n_species), nrow = n_species, ncol = n_patch)
  }

  # species niche
  v_mu <- runif(n = n_species, min = optim_min, max = optim_max)
  m_mu <- matrix(rep(x = v_mu, times = n_patch), nrow = n_species, ncol = n_patch)

  if (any(r0 < 1)) stop("r0 must be greater than or equal to one")
  if (length(r0) == 1) {
    v_r0 <- rep(x = r0, times = n_species)
    m_r0 <- matrix(rep(x = r0, times = n_species * n_patch), nrow = n_species, ncol = n_patch)
  } else {
    if (length(r0) != n_species) stop("r0 must have length of one or n_species")
    v_r0 <- r0
    m_r0 <- matrix(rep(x = v_r0, times = n_patch), nrow = n_species, ncol = n_patch)
  }

  # environmental variation among patches
  if (length(mu_env) == 1) {
    print("Only one environmental value is given: the model will assume environmental conditions are the same at all habitat patches")
    v_mu_z <- rep(x = mu_env, times = n_patch)
  } else {
    if (length(mu_env) != n_patch) stop("mu_env must have length of n_patch")
    v_mu_z <- mu_env
  }

  # interaction matrix
  if (interaction_type == "constant") {
    if (alpha < 0 | length(alpha) != 1) stop("invalid value of alpha - the value must be a positive scalar")
    alpha <- alpha
  }

  if (interaction_type == "random") {
    if (min_alpha < 0 | max_alpha < 0) stop("invalid values of min_alpha and/or max_alpha - values must be positive values")
    if (is.null(min_alpha) | is.null(max_alpha)) stop("provide min_alpha and max_alpha")
    if (min_alpha > max_alpha) stop("max_alpha must exceed min_alpha")
    alpha <- runif(n = n_species * n_species, min = min_alpha, max = max_alpha)
  }

  m_interaction <- matrix(alpha, nrow = n_species, ncol = n_species)
  diag(m_interaction) <- 1

  # dispersal matrix
  if (is.null(distance_matrix)) {
    print("Distance matrix is not given: generate a square landscape with landscape_size (default: 10) on a side")
    v_x_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    v_y_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    m_distance <- data.matrix(dist(cbind(v_x_coord, v_y_coord), diag = TRUE, upper = TRUE))
    m_dispersal <- data.matrix(exp(-theta * m_distance))
    diag(m_dispersal) <- 0
  } else {
    m_distance <- distance_matrix
    if (is.matrix(m_distance) == 0) stop("Distance matrix should be provided as matrix")
    if (nrow(m_distance) != n_patch) stop("Invalid dimension of distance matrix m_distance: distance_matrix must have a dimension of n_patch * n_patch")
    m_dispersal <- exp(-theta * m_distance)
    diag(m_dispersal) <- 0
  }

  if (length(p_dispersal) == 1) {
    print("Only one dispersal probability is given: the model will assume dispersal probability is the same for all species")
    v_p_dispersal <- rep(x = p_dispersal, times = n_species)
  } else {
    if (length(p_dispersal) != n_species) stop("p_dispersal must have length of one or n_species")
    v_p_dispersal <- p_dispersal
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
  colname <- c("abundance", "species", "niche_optim", "r_xt", "patch", "mean_env", "env", "timestep")
  m_dynamics <- matrix(NA, nrow = n_species * n_patch * n_timestep, ncol = length(colname))
  colnames(m_dynamics) <- colname
  st_row <- seq(from = 1, to = nrow(m_dynamics), by = n_species * n_patch)

  if (is.null(propagule_interval)) propagule_interval <- max(c(1, ceiling(n_warmup / 10)))
  propagule <- seq(from = propagule_interval, to = max(c(1, n_warmup)), by = propagule_interval)

  m_n <- matrix(rpois(n_species * n_patch, 0.5), nrow = n_species, ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (n in 1:n_sim) {
    if (n_warmup > 0) {
      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n_species * n_patch, 0.5), nrow = n_species, ncol = n_patch)
      }
    }

    m_z_xt <- matrix(rep(x = m_z[n, ], each = n_species), nrow = n_species, ncol = n_patch)
    m_r_xt <- m_r0 * exp(- ((m_mu - m_z_xt) / (sqrt(2) * sd_niche_width))^2)
    m_n_hat <- (m_n * m_r_xt) / (1 + ((m_r0 - 1) / m_k) * (m_interaction %*% m_n))

    m_e_hat <- m_n_hat * v_p_dispersal
    v_e_sum <- rowSums(m_e_hat)
    m_i_raw <- m_e_hat %*% m_dispersal
    v_i_sum <- rowSums(m_i_raw)
    v_i_sum[v_i_sum == 0] <- 1
    m_i_prob <- m_i_raw / v_i_sum
    m_i_hat <- m_i_prob * v_e_sum
    m_n_hat <- m_n_hat + m_i_hat - m_e_hat

    m_n <- matrix(rpois(n = n_species * n_patch, lambda = m_n_hat), nrow = n_species, ncol = n_patch)

    if (n > n_discard) {
      row_id <- st_row[n - n_discard]:(st_row[n - n_discard] + n_species * n_patch - 1)
      m_dynamics[row_id, ] <- cbind(c(m_n),
                                    rep(x = 1:n_species, times = n_patch),
                                    c(m_mu),
                                    c(m_r_xt),
                                    rep(x = 1:n_patch, each = n_species),
                                    rep(x = v_mu_z, each = n_species),
                                    c(m_z_xt),
                                    I(n - n_discard))
    }
    setTxtProgressBar(pb, n)
  }
  close(pb)


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {
    sample_patch <- sample(1:n_patch, size = min(c(n_patch, 5)), replace = FALSE)
    sample_species <- sample(1:n_species, size = min(c(n_species, 5)), replace = FALSE)

    g <- dplyr::as_tibble(m_dynamics) %>%
           dplyr::filter(.data$patch %in% sample_patch, .data$species %in% sample_species) %>%
           ggplot() +
           facet_grid(rows = vars(.data$species), cols = vars(.data$patch),
                      labeller = labeller(.rows = label_both, .cols = label_both)) +
           geom_line(mapping = aes(x = .data$timestep, y = .data$abundance, color = abs(.data$niche_optim - .data$env))) +
           scale_color_viridis_c(alpha = 0.8) +
           labs(color = "Environmental \ndeviation")
    print(g)
  }


  # return ------------------------------------------------------------------

  colnames(m_distance) <- rownames(m_distance) <- sapply(X = seq_len(n_patch), function(x) paste0("patch", x))
  colnames(m_interaction) <- rownames(m_interaction) <- sapply(X = seq_len(n_species), function(x) paste0("sp", x))

  species_df <- dplyr::as_tibble(m_dynamics) %>%
                  dplyr::group_by(species) %>%
                  dplyr::summarise(mean_abundance = mean(abundance)) %>%
                  dplyr::right_join(dplyr::tibble(species = 1:n_species, r0 = v_r0, niche_optim = v_mu, p_dispersal = v_p_dispersal), by = "species")

  return(list(dynamics_df = dplyr::as_tibble(m_dynamics),
              species_df = species_df,
              distance_matrix = m_distance,
              interaction_matrix = m_interaction))
}

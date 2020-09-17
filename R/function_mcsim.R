#' Simulate metacommunity dynamics
#'
#' @param n_species numeric value. Number of species in a metacommunity.
#' @param n_patch numeric value. Number of patches in a metacommunity.
#' @param n_warmup numeric value. Number of time-steps for warm-up. Default \code{200}.
#' @param n_burnin numeric value. Number of time-steps for burn-in. Default \code{200}.
#' @param n_timestep numeric value. Number of time-steps to be saved. Default \code{1000}.
#' @param propagule_interval numeric value. Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param carrying_capacity numeric value (length should be one or equal to \code{n_patch}). Carrying capacities of individual patches. Default \code{100}.
#' @param interaction_type character string. \code{"constant"} or \code{"random"}. \code{"constant"} assumes the unique interaction strength of \code{alpha} for all pairs of species. \code{"random"} draws random numbers from a uniform distribution with \code{min_alpha} and \code{max_alpha}.
#' @param alpha numeric value. Species interaction strength. Enabled if \code{interaction_type = "constant"}. Default \code{0}.
#' @param min_alpha numeric value. Minimum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param max_alpha numeric value. Maximum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param r0 numeric value (length should be one or equal to \code{n_species}). Maximum reproductive number of the Beverton-Holt model.
#' @param niche_optim numeric value (length should be one or equal to \code{n_species}). Niche optimum of species (environmental value that maximizes the reproductive number). Default \code{NULL}.
#' @param sd_niche_width numeric value (length should be one or equal to \code{n_species}). Niche width of species. Higher values indicate greater niche width.
#' @param optim_min numeric value. Minimum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species. Enabled if \code{niche_optim = NULL}.
#' @param optim_max numeric value. Maximum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species. Enabled if \code{niche_optim = NULL}.
#' @param xy_coord data frame. Each row should correspond to an individual patch, with x and y coordinates (columns). Defualt \code{NULL}.
#' @param distance_matrix numeric value. Distance matrix indicating distance between habitat patches. If \code{NULL}, a square landscape with randomly distributed patches will be generated. Default \code{NULL}.
#' @param landscape_size numeric value. Length of a landscape on a side. Enabled if \code{dispersal_matrix = NULL}.
#' @param mean_env numeric value (length should be one or equal to \code{n_patch}). Mean environmental values of patches.
#' @param sd_env numeric value. Standard deviation of temporal environmental variation at each patch.
#' @param spatial_env_cor logical. If \code{TRUE}, spatial autocorrelation in temporal environmental fluctuation is considered. Default \code{FALSE}.
#' @param phi numeric value. A parameter describing the distance decay of spatial autocorrelation in temporal environmental fluctuation. Enabled if \code{spatial_env_cor = TRUE}.
#' @param p_dispersal numeric value (length should be one or equal to \code{n_species}). Probability of dispersal.
#' @param theta numeric value. Dispersal parameter describing dispersal capability of species.
#' @param plot logical. If \code{TRUE}, five sample patches and species of \code{df_dynamics} are plotted.
#'
#' @return \code{df_dynamics} a data frame containing simulated metacommunity dynamics.
#' @return \code{df_species} a data frame containing species attributes.
#' @return \code{df_patch} a data frame containing patch attributes.
#' @return \code{df_diversity} a data frame containing diversity metrics.
#'
#' @importFrom dplyr %>% filter
#' @importFrom ggplot2 ggplot vars labeller geom_line aes scale_color_viridis_c labs facet_grid label_both
#' @importFrom stats dist rpois
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data
#'
#' @section Reference: see \href{https://github.com/aterui/mcbrnet}{github page} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' mcsim(n_patch = 5, n_species = 5)
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
                  niche_optim = NULL,
                  sd_niche_width = 1,
                  optim_min = -1,
                  optim_max = 1,
                  xy_coord = NULL,
                  distance_matrix = NULL,
                  landscape_size = 10,
                  mean_env = 0,
                  sd_env = 0.1,
                  spatial_env_cor = FALSE,
                  phi = 1,
                  p_dispersal = 0.1,
                  theta = 1,
                  plot = FALSE
) {

  # parameter setup ---------------------------------------------------------

  # carrying capacity
  if (length(carrying_capacity) == 1) {
    message("single value of carrying_capacity is given: assume carrying capacities are the same at all habitat patches")
    m_k <- matrix(carrying_capacity, nrow = n_species, ncol = n_patch)
  } else {
    if (length(carrying_capacity) != n_patch) stop("carrying_capacity must have length of one or n_patch")
    m_k <- matrix(rep(x = carrying_capacity, each = n_species), nrow = n_species, ncol = n_patch)
  }

  # species niche

  ## niche optimum
  if (is.null(niche_optim)) {
    message("niche_optim is not given: generate species niche optimum randomly")
    v_mu <- runif(n = n_species, min = optim_min, max = optim_max)
    m_mu <- matrix(rep(x = v_mu, times = n_patch), nrow = n_species, ncol = n_patch)
  } else {
    if (length(niche_optim) == 1) {
      message("single value of niche_optim is given: assume niche optimum are the same for all species")
      v_mu <- rep(x = niche_optim, times = n_species)
      m_mu <- matrix(niche_optim, nrow = n_species, ncol = n_patch)
    } else {
      if (length(niche_optim) != n_species) stop("niche_optim must have length of one or n_species")
      v_mu <- niche_optim
      m_mu <- matrix(rep(x = niche_optim, times = n_patch), nrow = n_species, ncol = n_patch)
    }
  }

  ## niche width
  if (length(sd_niche_width) == 1) {
    message("single value of sd_niche_width is given: assume niche width are the same for all species")
    m_sd_niche_width <- matrix(sd_niche_width, nrow = n_species, ncol = n_patch)
  } else {
    if (length(sd_niche_width) != n_species) stop("sd_niche_width must have length of one or n_species")
    m_sd_niche_width <- matrix(rep(x = sd_niche_width, times = n_patch), nrow = n_species, ncol = n_patch)
  }

  ## maximum reproductive number
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
  if (length(mean_env) == 1) {
    message("single value of mean_env is given: assume environmental conditions are the same at all habitat patches")
    v_mu_z <- rep(x = mean_env, times = n_patch)
  } else {
    if (length(mean_env) != n_patch) stop("mean_env must have length of one or n_patch")
    v_mu_z <- mean_env
  }

  # interaction matrix
  if (interaction_type == "constant") {
    if (alpha < 0 | length(alpha) != 1) stop("invalid value of alpha - the value must be a positive scalar")
    m_interaction <- matrix(alpha, nrow = n_species, ncol = n_species)
    diag(m_interaction) <- 1
  } else {
    if (interaction_type != "random") stop("invalid interaction_type")
    if (is.null(min_alpha) | is.null(max_alpha)) stop("provide min_alpha and max_alpha")
    if (min_alpha < 0 | max_alpha < 0) stop("invalid values of min_alpha and/or max_alpha - values must be positive values")
    if (min_alpha > max_alpha) stop("max_alpha must exceed min_alpha")
    alpha <- runif(n = n_species * n_species, min = min_alpha, max = max_alpha)
    m_interaction <- matrix(alpha, nrow = n_species, ncol = n_species)
    diag(m_interaction) <- 1
  }

  # dispersal matrix
  if (is.null(xy_coord) & is.null(distance_matrix)) {
    message("neither xy_coord nor distance_matrix is given: generate a square landscape with landscape_size (default: 10) on a side")
    v_x_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    v_y_coord <- runif(n = n_patch, min = 0, max = landscape_size)
    df_xy_coord <- dplyr::tibble(x_coord = v_x_coord, y_coord = v_y_coord)
    m_distance <- data.matrix(dist(df_xy_coord, diag = TRUE, upper = TRUE))
    m_dispersal <- data.matrix(exp(-theta * m_distance))
    diag(m_dispersal) <- 0
  } else {
    if (!is.null(xy_coord) & is.null(distance_matrix)) {
      if (nrow(xy_coord) != n_patch) stop("row numbers must match n_patch")
      if (ncol(xy_coord) != 2) stop("the number of columns must be two, describing x- and y-cooridnates")
      colnames(xy_coord) <- c("x_coord", "y_coord")
      df_xy_coord <- dplyr::as_tibble(xy_coord)
      m_distance <- data.matrix(dist(df_xy_coord, diag = TRUE, upper = TRUE))
      m_dispersal <- data.matrix(exp(-theta * m_distance))
      diag(m_dispersal) <- 0
    } else {
      if (!is.null(xy_coord)) message("both xy_coord and distance_matrix are given: argument xy_coord is ignored")
      if (is.matrix(distance_matrix) == 0) stop("distance matrix should be provided as matrix")
      if (nrow(distance_matrix) != n_patch) stop("invalid dimension: distance_matrix must have a dimension of n_patch * n_patch")
      df_xy_coord <- NULL
      m_distance <- distance_matrix
      m_dispersal <- exp(-theta * m_distance)
      diag(m_dispersal) <- 0
    }
  }

  if (length(p_dispersal) == 1) {
    message("single value of dispersal probability is given: assume dispersal probability is the same for all species")
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
  colname <- c("timestep", "patch", "mean_env", "env", "carrying_capacity",
               "species", "niche_optim", "r_xt", "abundance")
  m_dynamics <- matrix(NA, nrow = n_species * n_patch * n_timestep, ncol = length(colname))
  colnames(m_dynamics) <- colname
  st_row <- seq(from = 1, to = nrow(m_dynamics), by = n_species * n_patch)

  if (is.null(propagule_interval)) propagule_interval <- ceiling(n_warmup / 10)
  propagule <- seq(from = propagule_interval, to = max(c(1, n_warmup)), by = propagule_interval)

  m_n <- matrix(rpois(n = n_species * n_patch, lambda = 0.5), nrow = n_species, ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for (n in 1:n_sim) {
    if (n_warmup > 0) {
      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n = n_species * n_patch, lambda = 0.5), nrow = n_species, ncol = n_patch)
      }
    }

    m_z_xt <- matrix(rep(x = m_z[n, ], each = n_species), nrow = n_species, ncol = n_patch)
    m_r_xt <- m_r0 * exp(- ((m_mu - m_z_xt) / (sqrt(2) * m_sd_niche_width))^2)
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
      m_dynamics[row_id, ] <- cbind(I(n - n_discard),
                                    rep(x = 1:n_patch, each = n_species),
                                    rep(x = v_mu_z, each = n_species),
                                    c(m_z_xt),
                                    c(m_k),
                                    rep(x = 1:n_species, times = n_patch),
                                    c(m_mu),
                                    c(m_r_xt),
                                    c(m_n))
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
           labs(color = "Deviation \nfrom niche optimum")
    print(g)
  }


  # return ------------------------------------------------------------------

  colnames(m_distance) <- rownames(m_distance) <- sapply(X = seq_len(n_patch), function(x) paste0("patch", x))
  colnames(m_interaction) <- rownames(m_interaction) <- sapply(X = seq_len(n_species), function(x) paste0("sp", x))

  df_species <- dplyr::as_tibble(m_dynamics) %>%
                  dplyr::group_by(.data$species) %>%
                  dplyr::summarise(mean_abundance = mean(.data$abundance)) %>%
                  dplyr::right_join(dplyr::tibble(species = 1:n_species,
                                                  r0 = v_r0,
                                                  niche_optim = v_mu,
                                                  niche_width = sd_niche_width,
                                                  p_dispersal = v_p_dispersal),
                                    by = "species")

  df_dyn <- dplyr::as_tibble(m_dynamics)

  df_patch <- df_dyn %>%
                dplyr::group_by(.data$patch) %>%
                dplyr::summarise(alpha_div = sum(.data$abundance > 0) / n_timestep) %>%
                dplyr::left_join(dplyr::tibble(patch = seq_len(n_patch),
                                               mean_env = mean_env,
                                               carrying_capacity = carrying_capacity,
                                               connectivity = rowSums(m_dispersal)),
                                 by = "patch")

  alpha_div <- sum(df_dyn$abundance > 0) / (n_timestep * n_patch)
  gamma_div <- mean(tapply(X = df_dyn$species[df_dyn$abundance > 0],
                           INDEX = df_dyn$timestep[df_dyn$abundance > 0],
                           FUN = dplyr::n_distinct))
  beta_div <- gamma_div - alpha_div

  return(list(df_dynamics = df_dyn,
              df_species = df_species,
              df_patch = df_patch,
              df_diversity = dplyr::tibble(alpha_div, beta_div, gamma_div),
              df_xy_coord = df_xy_coord,
              distance_matrix = m_distance,
              interaction_matrix = m_interaction))
}

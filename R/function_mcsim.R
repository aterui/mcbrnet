#' Simulate metacommunity dynamics
#'
#' @param n_species Number of species in a metacommunity.
#' @param n_patch Number of patches in a metacommunity.
#' @param n_warmup Number of time-steps for warm-up. Default \code{200}.
#' @param n_burnin Number of time-steps for burn-in. Default \code{200}.
#' @param n_timestep Number of time-steps to be saved. Default \code{1000}.
#' @param propagule_interval Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param propagule_seed Propagule mean density (intensity parameter in a Poisson distribution). Default \code{0.5}.
#' @param carrying_capacity Carrying capacities of individual patches. Length must be one or equal to \code{n_patch}. Default \code{100}.
#' @param xi Hassell exponent. Undercompensation (xi < 1), no compensation (xi = 1; reduced to Beverton-Holt), and overcompensation (xi > 1).
#' @param interaction_type \code{"constant"} or \code{"random"}. \code{"constant"} assumes the unique interaction strength of \code{alpha} for all pairs of species. \code{"random"} draws random numbers from a uniform distribution with \code{min_alpha} and \code{max_alpha}.
#' @param alpha Species interaction strength. Enabled if \code{interaction_type = "constant"}. Default \code{0}.
#' @param min_alpha Minimum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param max_alpha Maximum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param r0 Maximum reproductive number of the Beverton-Holt model. Length must be one or equal to \code{n_species}.
#' @param niche_optim Niche optimum of species (environmental value that maximizes the reproductive number). Length must be one or equal to \code{n_species}. Default \code{NULL}.
#' @param min_optim Minimum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species. Enabled if \code{niche_optim = NULL}.
#' @param max_optim Maximum value of a uniform distribution that generates optimal environmental values of simulated species. Values are randomly assigned to species. Enabled if \code{niche_optim = NULL}.
#' @param sd_niche_width Niche width of species. Higher values indicate greater niche width. Length must be one or equal to \code{n_species}.
#' @param min_niche_width Minimum value of a uniform distribution that generates niche width values of simulated species. Values are randomly assigned to species. Enabled if \code{sd_niche_width = NULL}.
#' @param max_niche_width Maximum value of a uniform distribution that generates niche width values of simulated species. Values are randomly assigned to species. Enabled if \code{sd_niche_width = NULL}.
#' @param niche_cost Determine the cost of wide niche (smaller values imply greater costs of wider niche). Default \code{1}.
#' @param xy_coord Site coordinates. Must be provided as a data frame in which each row corresponds to an individual site with x and y coordinates (columns). Defualt \code{NULL}.
#' @param distance_matrix Distance matrix indicating distance between habitat patches. If provided, the distance matrix will be used to generate dispersal matrix and to calculate distance decay of environmental correlations. Default \code{NULL}.
#' @param dispersal_matrix Dispersal matrix to be used to simulate dispersal process. Override distance_matrix. Default \code{NULL}.
#' @param p_disturb Disturbance probability.
#' @param i_disturb Disturbance intensity expressed as proportional mortality (0 to 1). Length must be one or equal to \code{n_patch}.
#' @param landscape_size Length of a landscape on a side. Enabled if \code{dispersal_matrix = NULL}.
#' @param mean_env Mean environmental values of patches. Length must be one or equal to \code{n_patch}.
#' @param sd_env Standard deviation of temporal environmental variation at each patch.
#' @param spatial_env_cor If \code{TRUE}, spatial autocorrelation in temporal environmental fluctuation is considered. Default \code{FALSE}.
#' @param phi Parameter describing the distance decay of spatial autocorrelation in temporal environmental fluctuation. Enabled if \code{spatial_env_cor = TRUE}.
#' @param p_dispersal Dispersal probability. Length must be one or equal to \code{n_species}.
#' @param theta Rate parameter of exponential dispersal kernel.
#' @param zeta Species sensitivity to environmental pollutants.
#' @param q Concentration of environmental pollutants.
#' @param plot If \code{TRUE}, five sample patches and species of \code{df_dynamics} are plotted.
#'
#' @importFrom dplyr %>% filter
#' @importFrom ggplot2 ggplot vars labeller geom_line aes scale_color_viridis_c labs facet_grid label_both
#' @importFrom stats dist rbinom rgeom rnorm rpois runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/index.html}{package webpage} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' \dontrun{
#' # not run
#' mcsim(n_patch = 5, n_species = 5)
#' }
#'
#' @export
#'
mcsim <- function(n_species = 5,
                  n_patch = 5,
                  n_warmup = 200,
                  n_burnin = 200,
                  n_timestep = 1000,
                  propagule_interval = NULL,
                  propagule_seed = 0.5,
                  carrying_capacity = 100,
                  xi = 1,
                  interaction_type = "constant",
                  alpha = 0,
                  min_alpha = NULL,
                  max_alpha = NULL,
                  r0 = 4,
                  niche_optim = NULL,
                  min_optim = -1,
                  max_optim = 1,
                  sd_niche_width = NULL,
                  min_niche_width = 0.1,
                  max_niche_width = 1,
                  niche_cost = 1,
                  xy_coord = NULL,
                  distance_matrix = NULL,
                  dispersal_matrix = NULL,
                  p_disturb = 0,
                  i_disturb = 0,
                  landscape_size = 10,
                  mean_env = 0,
                  sd_env = 0.1,
                  spatial_env_cor = FALSE,
                  phi = 1,
                  p_dispersal = 0.1,
                  theta = 1,
                  zeta = 0,
                  q = 0,
                  plot = FALSE
) {

  # define function ---------------------------------------------------------

  fun_r <- function(r0, z, mu, sd, nu) {
    r0 * exp(- ((sd^2) / (2 * nu^2))) * exp(- ((mu - z) / (sqrt(2) * sd))^2)
  }

  # parameter setup ---------------------------------------------------------

  # carrying capacity ####
  ## internal function; see "fun_to_m.R"
  list_k <- fun_to_m(x = carrying_capacity,
                     n_species = n_species,
                     n_patch = n_patch,
                     param_attr = "patch")

  m_k <- list_k$m_x

  ## niche optimum ####
  ## internal function; see "fun_to_m.R"
  if (is.null(niche_optim)) {

    message("niche_optim is not given: generate species niche optimum randomly")

    v_mu <- runif(n = n_species,
                  min = min_optim,
                  max = max_optim)

    m_mu <- matrix(rep(x = v_mu,
                       times = n_patch),
                   nrow = n_species,
                   ncol = n_patch)

  } else {

    list_mu <- fun_to_m(x = niche_optim,
                        n_species = n_species,
                        n_patch = n_patch,
                        param_attr = "species")

    v_mu <- list_mu$v_x
    m_mu <- list_mu$m_x

  }

  ## niche width ####
  ## internal function; see "fun_to_m.R"
  if (is.null(sd_niche_width)) {

    message("sd_niche_width is not given: generate species niche width randomly")

    v_sd_niche_width <- runif(n = n_species,
                              min = min_niche_width,
                              max = max_niche_width)

    m_sd_niche_width <- matrix(rep(x = v_sd_niche_width,
                                   times = n_patch),
                               nrow = n_species,
                               ncol = n_patch)

  } else {

    list_sd_niche_width <- fun_to_m(x = sd_niche_width,
                                    n_species = n_species,
                                    n_patch = n_patch,
                                    param_attr = "species")

    v_sd_niche_width <- list_sd_niche_width$v_x
    m_sd_niche_width <- list_sd_niche_width$m_x

  }

  ## maximum reproductive number ####
  ## internal function; see "fun_to_m.R"
  if (any(r0 < 1)) stop("r0 must be greater than or equal to one")

  list_r0 <- fun_to_m(x = r0,
                      n_species = n_species,
                      n_patch = n_patch,
                      param_attr = "species")

  v_r0 <- list_r0$v_x
  m_r0 <- list_r0$m_x

  ## environmental variation among patches ####
  ## internal function; see "fun_to_m.R"
  list_mu_z <- fun_to_m(x = mean_env,
                        n_species = n_species,
                        n_patch = n_patch,
                        param_attr = "patch")

  v_mu_z <- list_mu_z$v_x

  ## interaction matrix ####
  ## internal function; see "fun_int_mat.R"
  m_interaction <- fun_int_mat(n_species = n_species,
                               alpha = alpha,
                               min_alpha = min_alpha,
                               max_alpha = max_alpha,
                               interaction_type = interaction_type)

  ## dispersal matrix ####
  ## internal function; see "fun_disp_mat.R"
  list_dispersal <- fun_disp_mat(n_patch = n_patch,
                                 landscape_size = landscape_size,
                                 theta = theta,
                                 xy_coord = xy_coord,
                                 distance_matrix = distance_matrix,
                                 dispersal_matrix = dispersal_matrix)

  m_distance <- list_dispersal$m_distance
  m_dispersal <- list_dispersal$m_dispersal
  df_xy_coord <- list_dispersal$df_xy_coord

  ## internal function; see "fun_to_v.R"
  v_p_dispersal <- fun_to_v(x = p_dispersal,
                            n = n_species)

  ## disturbance ####
  if (p_disturb > 1 | p_disturb < 0) stop("p_disturb must be 0 to 1")
  if (any(i_disturb > 1) | any(i_disturb < 0)) stop("i_disturb must be 0 to 1")

  v_i_disturb <- fun_to_v(x = i_disturb,
                          n = n_patch)

  ## local pollutants ####
  if (zeta < 0) stop("zeta must be > 0")

  v_zeta <- fun_to_v(x = zeta,
                     n = n_species)

  v_q <- fun_to_v(x = q,
                  n = n_patch)

  # dynamics ----------------------------------------------------------------

  ## object setup ####

  ## number of replicates
  n_sim <- n_warmup + n_burnin + n_timestep
  n_discard <- n_warmup + n_burnin

  ## environment
  var_env <- sd_env^2

  if (spatial_env_cor == TRUE) {

    if (is.null(m_distance)) stop("Provide distance matrix to model spatial environmental autocorrelation")
    m_sigma <- var_env * exp(-phi * m_distance)

  }else{

    m_sigma <- var_env * diag(x = 1, nrow = n_patch, ncol = n_patch)

  }

  m_z <- mvtnorm::rmvnorm(n = n_sim,
                          mean = v_mu_z,
                          sigma = m_sigma)

  ## community object
  colname <- c("timestep",
               "patch_id",
               "mean_env",
               "env",
               "carrying_capacity",
               "species",
               "niche_optim",
               "r_xt",
               "abundance")

  m_dynamics <- matrix(NA,
                       nrow = n_species * n_patch * n_timestep,
                       ncol = length(colname))

  colnames(m_dynamics) <- colname

  st_row <- seq(from = 1,
                to = nrow(m_dynamics),
                by = n_species * n_patch)

  ## propagule
  if (is.null(propagule_interval)) {

    propagule_interval <- ceiling(n_warmup / 10)

  }

  propagule <- seq(from = propagule_interval,
                   to = max(c(1, n_warmup)),
                   by = propagule_interval)

  ## disturbance incidence
  psi <- rbinom(n = n_sim,
                size = 1,
                prob = p_disturb)
  psi[1:n_warmup] <- 0

  ## initial values
  m_n <- matrix(rpois(n = n_species * n_patch,
                      lambda = 0.5),
                nrow = n_species,
                ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

  ## temporal process ####
  for (n in seq_len(n_sim)) {

    if (n_warmup > 0) {

      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n = n_species * n_patch,
                                  lambda = propagule_seed),
                            nrow = n_species,
                            ncol = n_patch)
      }

    }

    ## time-specific local environment
    m_z_xt <- matrix(rep(x = m_z[n, ],
                         each = n_species),
                     nrow = n_species,
                     ncol = n_patch)

    ## time-specific intrinsic growth rate
    m_r_xt <- fun_r(r0 = m_r0,
                    mu = m_mu,
                    z = m_z_xt,
                    sd = m_sd_niche_width,
                    nu = niche_cost)

    m_r_xt <- m_r_xt * exp(-outer(v_zeta, v_q))

    ## internal community dynamics with competition
    m_n_hat <- (m_n * m_r_xt) / (1 + (m_r0 - 1) * ((m_interaction %*% m_n) / m_k)^xi)

    ## disturbance
    m_n_bar <- t(t(m_n_hat) * (1 - psi[n] * v_i_disturb))

    ## dispersal, internal function: see "fun_dispersal.R"
    m_n_prime <- fun_dispersal(x = m_n_bar,
                               v_p_dispersal = v_p_dispersal,
                               m_dispersal = m_dispersal)

    ## demographic stochasticity
    m_n <- matrix(rpois(n = n_species * n_patch,
                        lambda = m_n_prime),
                  nrow = n_species,
                  ncol = n_patch)

    if (n > n_discard) {

      row_id <- seq(from = st_row[n - n_discard],
                    to = st_row[n - n_discard] + n_species * n_patch - 1,
                    by = 1)

      m_dynamics[row_id, ] <- cbind(
        # timestep
        I(n - n_discard),
        # patch_id
        rep(x = seq_len(n_patch), each = n_species),
        # mean_env
        rep(x = v_mu_z, each = n_species),
        # env
        c(m_z_xt),
        # carrying_capacity
        c(m_k),
        # species
        rep(x = seq_len(n_species), times = n_patch),
        # niche_optim
        c(m_mu),
        # r_xt
        c(m_r_xt),
        # abundance
        c(m_n))

    }

    setTxtProgressBar(pb, n)

  }

  close(pb)


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {

    sample_patch <- sample(seq_len(n_patch),
                           size = min(c(n_patch, 5)),
                           replace = FALSE)

    sample_species <- sample(seq_len(n_species),
                             size = min(c(n_species, 5)),
                             replace = FALSE)

    g <- dplyr::as_tibble(m_dynamics) %>%
      dplyr::filter(.data$patch_id %in% sample_patch,
                    .data$species %in% sample_species) %>%
      ggplot() +
      facet_grid(rows = vars(.data$species),
                 cols = vars(.data$patch_id),
                 labeller = labeller(.rows = label_both,
                                     .cols = label_both)) +
      geom_line(mapping = aes(x = .data$timestep,
                              y = .data$abundance,
                              color = abs(.data$niche_optim - .data$env))) +
      scale_color_viridis_c(alpha = 0.8) +
      labs(color = "Deviation \nfrom niche optimum")

    print(g)

  }


  # return ------------------------------------------------------------------

  # dynamics
  df_dyn <- dplyr::as_tibble(m_dynamics)

  # species attributes
  df_species <- df_dyn %>%
    dplyr::group_by(.data$species) %>%
    dplyr::summarise(mean_abundance = mean(.data$abundance)) %>%
    dplyr::right_join(dplyr::tibble(species = seq_len(n_species),
                                    r0 = v_r0,
                                    niche_optim = v_mu,
                                    sd_niche_width = v_sd_niche_width,
                                    p_dispersal = v_p_dispersal),
                      by = "species")

  # patch attributes
  df_patch <- df_dyn %>%
    dplyr::group_by(.data$patch_id) %>%
    dplyr::summarise(alpha_div = sum(.data$abundance > 0) / n_timestep) %>%
    dplyr::left_join(dplyr::tibble(patch_id = seq_len(n_patch),
                                   mean_env = mean_env,
                                   carrying_capacity = carrying_capacity,
                                   disturbance = v_i_disturb),
                     by = "patch_id")

  # diversity metrics
  alpha_div <- sum(df_dyn$abundance > 0) / (n_timestep * n_patch)

  gamma_div <- df_dyn %>%
    dplyr::group_by(.data$timestep) %>%
    dplyr::summarise(gamma_t = dplyr::n_distinct(.data$species[.data$abundance > 0])) %>%
    dplyr::pull(.data$gamma_t) %>%
    mean()

  if (gamma_div == 0 | alpha_div == 0) {
    beta_div <- NA
  } else {
    beta_div <- gamma_div / alpha_div
  }

  # return
  return(list(df_dynamics = df_dyn,
              df_species = df_species,
              df_patch = df_patch,
              df_diversity = dplyr::tibble(alpha_div, beta_div, gamma_div),
              df_xy_coord = df_xy_coord,
              distance_matrix = m_distance,
              interaction_matrix = m_interaction))
}

#' Simulate meta-food web dynamics with intraguild predation
#'
#' @param n_patch Number of patches in a metacommunity.
#' @param n_warmup Number of time-steps for warm-up. Default \code{200}.
#' @param n_burnin Number of time-steps for burn-in. Default \code{200}.
#' @param n_timestep Number of time-steps to be saved. Default \code{1000}.
#' @param r_b Maximum reproductive rate of basal species
#' @param conv_eff Energetic conversion efficiency; must be given by the order of basal to ig-prey, basal to ig-predator, and ig-prey to ig-predator
#' @param attack_rate Attack rate. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param handling_time Handling time. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param s0 Background survival rate. Must be given by the order of basal, ig-prey, and ig-predator
#' @param propagule_interval Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param carrying_capacity Carrying capacities of individual patches. Length must be one or equal to \code{n_patch}. Default \code{100}.
#' @param xy_coord Data frame for site coordinates. Each row should correspond to an individual patch, with x and y coordinates (columns). Defualt \code{NULL}.
#' @param distance_matrix Distance matrix indicating distance between habitat patches. If provided, the distance matrix will be used to generate dispersal matrix and to calculate distance decay of environmental correlations. Default \code{NULL}.
#' @param dispersal_matrix Dispersal matrix to be used to simulate dispersal process. Override distance_matrix. Default \code{NULL}.
#' @param p_disturb Disturbance probability.
#' @param m_disturb Disturbance-induced proportional mortality.
#' @param landscape_size Length of a landscape on a side. Enabled if \code{dispersal_matrix = NULL}.
#' @param p_dispersal Probability of dispersal. Length must be one or equal to \code{n_species}.
#' @param theta Dispersal parameter describing dispersal capability of species.
#' @param plot If \code{TRUE}, five sample patches and species of \code{df_dynamics} are plotted.
#'
#' @importFrom dplyr %>% filter
#' @importFrom ggplot2 ggplot vars labeller geom_line aes scale_color_viridis_c labs facet_grid label_both
#' @importFrom stats dist rpois
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

igpsim <- function(n_patch = 5,
                   n_warmup = 200,
                   n_burnin = 200,
                   n_timestep = 1000,
                   r_b = 5,
                   conv_eff = rep(5, 3),
                   attack_rate = rep(0.5, 3),
                   handling_time = rep(1, 3),
                   s0 = rep(0.8, 3),
                   propagule_interval = NULL,
                   carrying_capacity = 100,
                   xy_coord = NULL,
                   distance_matrix = NULL,
                   dispersal_matrix = NULL,
                   p_disturb = 0,
                   m_disturb = 0,
                   landscape_size = 10,
                   p_dispersal = 0.1,
                   theta = 1,
                   plot = FALSE
) {

  # parameter setup ---------------------------------------------------------

  ## number of species ####
  n_species <- 3

  ## functional response ####
  v_a <- fun_to_v(x = attack_rate,
                  n = n_species)

  v_e <- fun_to_v(x = conv_eff,
                  n = n_species)

  v_h <- fun_to_v(x = handling_time,
                  n = n_species)

  ## carrying capacity ####
  v_k <- fun_to_v(x = carrying_capacity,
                  n = n_patch)

  ## background survival ####
  v_s0 <- fun_to_v(x = s0,
                   n = n_species)

  ## dispersal matrix ####
  list_dispersal <- fun_disp_mat(n_patch = n_patch,
                                 landscape_size = landscape_size,
                                 theta = theta,
                                 xy_coord = xy_coord,
                                 distance_matrix = distance_matrix,
                                 dispersal_matrix = dispersal_matrix)

  m_distance <- list_dispersal$m_distance
  m_dispersal <- list_dispersal$m_dispersal
  df_xy_coord <- list_dispersal$df_xy_coord

  v_p_dispersal <- fun_to_v(x = p_dispersal,
                            n = n_species)

  ## disturbance ####
  if (p_disturb > 1 | p_disturb < 0) stop("p_disturb must be 0 to 1")
  if (any(m_disturb > 1) | any(m_disturb < 0)) stop("m_disturb must be 0 to 1")

  v_disturb <- fun_to_v(x = m_disturb,
                        n = n_patch)

  # dynamics ----------------------------------------------------------------

  ## object setup ####

  ## number of replicates
  n_sim <- n_warmup + n_burnin + n_timestep
  n_discard <- n_warmup + n_burnin

  ## community objects
  colname <- c("timestep",
               "patch_id",
               "carrying_capacity",
               "species",
               "abundance",
               "fcl")

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
                                  lambda = 0.5),
                            nrow = n_species,
                            ncol = n_patch)
      }

    }

    ## internal community dynamics with intraguild predation
    list_n_hat <- fun_igp(x = m_n,
                          r_b = r_b,
                          k = v_k,
                          e = v_e,
                          a = v_a,
                          h = v_h,
                          s0 = v_s0)

    ## disturbance
    m_n_hat <- t(t(list_n_hat$m_n_hat) * (1 - psi[n] * v_disturb))

    ## dispersal, internal function: see "fun_dispersal.R"
    m_n_prime <- fun_dispersal(x = m_n_hat,
                               v_p_dispersal = v_p_dispersal,
                               m_dispersal = m_dispersal)

    ## demographic stochasticity
    m_n <- matrix(rpois(n = n_species * n_patch,
                        lambda = m_n_prime),
                  nrow = n_species,
                  ncol = n_patch)

    ## food chain length
    v_fcl <- fun_get_fcl(x = m_n,
                         delta = list_n_hat$delta)

    if (n > n_discard) {

      row_id <- seq(from = st_row[n - n_discard],
                    to = st_row[n - n_discard] + n_species * n_patch - 1,
                    by = 1)

      m_dynamics[row_id, ] <- cbind(# timestep
                                    I(n - n_discard),
                                    # patch id
                                    rep(x = seq_len(n_patch), each = n_species),
                                    # carrying capacity
                                    c(v_k),
                                    # species id
                                    rep(x = seq_len(n_species), times = n_patch),
                                    # abundance
                                    c(m_n),
                                    # food chain length
                                    rep(x = v_fcl, each = n_species))

    }

    setTxtProgressBar(pb, n)

  }

  close(pb)


  # visualization -----------------------------------------------------------

  # dynamics
  df_dyn <- dplyr::as_tibble(m_dynamics) %>%
    dplyr::mutate(species = dplyr::case_when(.data$species == 1 ~ "basal",
                                             .data$species == 2 ~ "ig-prey",
                                             .data$species == 3 ~ "ig-predator"),
                  species = factor(.data$species,
                                   levels = c("basal",
                                              "ig-prey",
                                              "ig-predator")))

  if (plot == TRUE) {

    sample_patch <- sample(seq_len(n_patch),
                           size = min(c(n_patch, 5)),
                           replace = FALSE)

    sample_species <- sample(seq_len(n_species),
                             size = min(c(n_species, 5)),
                             replace = FALSE)

    g <- df_dyn %>%
      dplyr::filter(.data$patch_id %in% sample_patch) %>%
      ggplot() +
      facet_grid(rows = vars(.data$species),
                 cols = vars(.data$patch_id),
                 labeller = labeller(.rows = label_both,
                                     .cols = label_both)) +
      geom_line(mapping = aes(x = .data$timestep,
                              y = .data$abundance,
                              color = .data$species))

    print(g)

  }


  # return ------------------------------------------------------------------

  # species attributes
  df_species <- dplyr::as_tibble(m_dynamics) %>%
    dplyr::group_by(.data$species) %>%
    dplyr::summarise(mean_abundance = mean(.data$abundance)) %>%
    dplyr::right_join(dplyr::tibble(species = seq_len(n_species),
                                    s0 = s0,
                                    p_dispersal = v_p_dispersal),
                      by = "species") %>%
    dplyr::mutate(species = dplyr::case_when(.data$species == 1 ~ "basal",
                                             .data$species == 2 ~ "ig-prey",
                                             .data$species == 3 ~ "ig-predator"),
                  species = factor(.data$species,
                                   levels = c("basal",
                                              "ig-prey",
                                              "ig-predator")))

  # patch attributes
  df_patch <- df_dyn %>%
    dplyr::group_by(.data$patch_id) %>%
    dplyr::summarise(fcl = mean(.data$fcl)) %>%
    dplyr::left_join(dplyr::tibble(patch_id = seq_len(n_patch),
                                   carrying_capacity = carrying_capacity,
                                   disturbance = v_disturb),
                     by = "patch_id")

  # interaction
  df_int <- dplyr::tibble(interaction = c("ig-prey on basal",
                                          "ig-predator on basal",
                                          "ig-predator on ig-prey"),
                          conv_eff = v_e,
                          attack_rate = v_a,
                          handling_time = v_h)

  # return
  return(list(df_dynamics = df_dyn,
              df_species = df_species,
              df_patch = df_patch,
              df_int = df_int,
              df_xy_coord = df_xy_coord,
              distance_matrix = m_distance)
  )
}

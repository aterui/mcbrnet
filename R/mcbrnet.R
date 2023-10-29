#' Generate a random branching network
#'
#' @param n_patch Number of patches in a network.
#' @param p_branch Branching probability (success probability of a geometric distribution).
#' @param mean_env_source Mean value of environmental condition at upstream terminals.
#' @param sd_env_source SD of environmental condition at upstream terminals.
#' @param rho Strength of spatial autocorrelation in environmental condition.
#' @param sd_env_lon SD of longitudinal environmental noise.
#' @param mean_disturb_source Mean disturbance strength at headwaters. The value is assumed to represent the proportional mortality (0 - 1.0) at the patch level.
#' @param sd_disturb_source SD of disturbance strength at headwaters. The SD is defined in a logit scale with a normal distribution.
#' @param sd_disturb_lon SD of longitudinal noise of disturbance strength. The SD is defined in a logit scale with a normal distribution.
#' @param randomize_patch Whether randomize patches or not. If \code{FALSE}, the function may generate a biased network with ordered patches. Default \code{TRUE}.
#' @param n_patch_free Whether imposing a constraint on \code{n_patch}. If \code{TRUE}, the number of patches a random variable following a negative binomial distribution.
#'
#' @importFrom dplyr %>%
#' @importFrom grDevices grey
#' @importFrom graphics par text
#' @importFrom stats complete.cases rbinom rgeom rnorm runif
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' \dontrun{
#' # not run
#' brnet(n_patch = 10, p_branch = 0.5)
#' }
#'
#' @export

brnet <- function(n_patch = 50,
                  p_branch = 0.5,
                  mean_env_source = 0,
                  sd_env_source = 1,
                  rho = 1,
                  sd_env_lon = 0.1,
                  mean_disturb_source = 0.9,
                  sd_disturb_source = 1,
                  sd_disturb_lon = 0.1,
                  randomize_patch = TRUE,
                  n_patch_free = FALSE) {


  # define variables --------------------------------------------------------

  ## internal function: see "fun_get_n_branch.R"
  n_branch <- fun_get_n_branch(n_patch = n_patch,
                               p_branch = p_branch)


  # adjacency matrix --------------------------------------------------------

  list_adj <- fun_m_adj(n_patch = n_patch,
                        p_branch = p_branch,
                        n_branch = n_branch,
                        n_patch_free = n_patch_free)

  v_n_patch_branch <- list_adj$v_n_patch_branch
  m_adj <- list_adj$m_adj


  # distance matrix ---------------------------------------------------------

  ## internal function: see "function_adjtodist.R"
  m_distance <- m_adj %>%
    igraph::graph_from_adjacency_matrix("undirected") %>%
    igraph::distances()


  # upstream watershed area ------------------------------------------------

  ## internal function: see "fun_wa.R"
  m_wa <- fun_wa(x = m_adj)
  v_wa <- rowSums(m_wa)


  # environmental condition -------------------------------------------------

  ## internal function: see "fun_patch_attr.R"
  v_env <- fun_patch_attr(x = m_adj,
                          n_branch = n_branch,
                          mean_source = mean_env_source,
                          sd_source = sd_env_source,
                          sd_lon = sd_env_lon,
                          m_distance = m_distance,
                          rho = rho,
                          v_wa = v_wa,
                          logit = FALSE)


  # disturbance -------------------------------------------------------------

  if (!(mean_disturb_source <= 1 && mean_disturb_source >= 0)) {

    stop("mean_disturb_source must be between 0 and 1")

  }

  ## internal function: see "fun_patch_attr.R"
  v_disturb <- fun_patch_attr(x = m_adj,
                              n_branch = n_branch,
                              mean_source = mean_disturb_source,
                              sd_source = sd_disturb_source,
                              sd_lon = sd_disturb_lon,
                              m_distance = m_distance,
                              rho = 1,
                              v_wa = v_wa,
                              logit = TRUE)


  # randomize nodes ---------------------------------------------------------

  branch <- unlist(lapply(seq_len(n_branch),
                          function(x) rep(x, each = v_n_patch_branch[x])))

  patch <- seq_len(n_patch)

  if (randomize_patch == TRUE) {

    if (n_branch > 1) {

      df_id <- dplyr::tibble(branch = as.character(c(1, resample(2:n_branch)))) %>%
        dplyr::left_join(data.frame(patch,
                                    branch = as.character(branch)),
                         by = "branch")
      v_wa <- v_wa[df_id$patch]
      v_env <- v_env[df_id$patch]
      v_disturb <- v_disturb[df_id$patch]
      m_adj <- m_adj[df_id$patch, df_id$patch]
      m_distance <- m_distance[df_id$patch, df_id$patch]

    } else {

      df_id <- dplyr::tibble(branch = 1,
                             patch = seq_len(n_patch))

    }
  } else {

    df_id <- dplyr::tibble(branch = branch,
                           patch = patch)

  }

  df_patch <- dplyr::tibble(patch_id = seq_len(n_patch),
                            branch_id = as.numeric(df_id$branch),
                            environment = c(v_env),
                            disturbance = c(boot::inv.logit(v_disturb)),
                            n_patch_upstream = c(v_wa))


  # return ------------------------------------------------------------------

  list_brnet <- list(adjacency_matrix = m_adj,
                     distance_matrix = m_distance,
                     df_patch = df_patch)

  class(list_brnet) <- "brnet"

  return(list_brnet)
}

#' Create point source disturbance
#'
#' @param x `brnet()` object
#' @param n_source Number of point sources
#' @param p Reduction factor for downstream propagation
#' @param q Reduction factor for upstream propagation
#' @param pattern Spatial pattern of point sources. "random", "cluster", "upstream", "downstream"
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @importFrom dplyr %>% filter
#' @importFrom rlang .data
#'
#' @export

ptsource <- function(x,
                     n_source,
                     p,
                     q,
                     pattern = "random") {


  # parameter check ---------------------------------------------------------

  if (p < 0 || p > 1 || q < 0 || q > 1) stop("p & q must be 0 - 1")
  if (nrow(x$distance_matrix) < n_source) stop("n_source must be < n_patch")
  if (!inherits(x, what = "brnet")) stop("x must be a 'brnet' object")


  # prepare basic objects ---------------------------------------------------

  ## basic objects
  v_patch_order <- x$distance_matrix[1, ]
  n_patch <- nrow(x$adjacency_matrix)

  v_d_to_root <- sort(x$distance_matrix[1, ])
  m_adj <- x$adjacency_matrix[order(v_patch_order), order(v_patch_order)]
  m_dist <- x$distance_matrix[order(v_patch_order), order(v_patch_order)]
  v_wa <- x$df_patch$n_patch_upstream[order(v_patch_order)]

  # extract upstream adjacency
  m_adj_up <- m_adj
  m_adj_up[lower.tri(m_adj_up)] <- 0

  # upstream tributary proportion
  m_wa_prop <- t(apply(X = m_adj_up,
                       MARGIN = 1,
                       FUN = function(k) {
                         (k * v_wa) / ifelse(sum(k) == 0,
                                             yes = 1,
                                             no = sum(k * v_wa))
                       }))


  # impact ------------------------------------------------------------------

  if (pattern == "random") {
    pts <- resample(seq_len(n_patch), size = n_source)
  }

  if (pattern == "cluster") {
    s <- resample(seq_len(n_patch), size = 1)
    distance <- m_dist[s, ]
    distance[distance == 0] <- 2
    prob <- 1 / distance
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  if (pattern == "downstream") {
    prob <- v_wa
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  if (pattern == "upstream") {
    prob <- 1 / v_wa
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  v_d <- v_u <- v_d_attr <- v_u_attr <- rep(0, n_patch)
  v_d[pts] <- v_u[pts] <- 1
  v_d_attr[pts] <- v_u_attr[pts] <- 1

  for (i in seq_len(max(v_d_to_root))) {
    ## upstream propagation
    v_u <- t(m_wa_prop) %*% (q * v_u)
    v_u_attr <- v_u + v_u_attr

    ## downstream propagation
    v_d <- m_wa_prop %*% (p * v_d)
    v_d_attr <- v_d + v_d_attr
  }

  v_u_attr[pts] <- v_u_attr[pts] - 1
  v_attr <- v_u_attr + v_d_attr


  # export ------------------------------------------------------------------

  v_source <- rep(0, n_patch)
  v_source[pts] <- 1

  df0 <- dplyr::tibble(patch_id = x$df_patch$patch_id[order(v_patch_order)],
                       impact = c(v_attr),
                       point_source = v_source) %>%
    dplyr::arrange(.data$patch_id) %>%
    dplyr::right_join(y = x$df_patch,
                      by = "patch_id")

  x$df_patch <- df0

  return(x)

}

#' Dispersal probability matrix after accounting for network fragmentation
#'
#' @param x 'brnet' object or adjacency matrix
#' @param rate rate parameter of exponential dispersal kernel
#' @param pattern fragmentation pattern
#' @param p passability of fragmented edges (probability)
#' @param n_barrier number of barriers
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export

frgm <- function(x,
                 rate,
                 pattern = "random",
                 p = NULL,
                 n_barrier) {


  # adjacency graph & distance matrix ---------------------------------------

  if (inherits(x, what = "brnet")) {
    m_adj <- x$adjacency_matrix
  } else {
    m_adj <- x
  }

  if (dplyr::n_distinct(dim(m_adj)) != 1) stop("invalid dimension in the adjacency matrix")

  g0 <- m_adj %>%
    igraph::graph.adjacency("undirected")

  df_g0 <- igraph::as_data_frame(g0)

  if (inherits(x, what = "brnet")) {
    m_dist <- x$distance_matrix
  } else {
    m_dist <- igraph::distances(g0)
  }


  # basic numbers -----------------------------------------------------------

  n_patch <- unique(dim(m_adj))
  n_edge <- n_patch - 1

  m_disp <- exp(-rate * m_dist)
  diag(m_disp) <- 0


  # passability -------------------------------------------------------------

  if (any(p < 0 | p > 1)) stop("p must be 0 - 1.")

  if (n_barrier > n_edge) stop("n_barrier exceeds the number of edges in the graph")

  if (pattern == "random") {

    barrier <- resample(seq_len(n_edge), size = n_barrier)

  } else {

    if (pattern == "downstream" || pattern == "upstream") {

      if (!inherits(x, what = "brnet")) stop("x must be class 'brnet'")

      v_wa <- df_g0 %>%
        dplyr::as_tibble() %>%
        dplyr::left_join(x$df_patch,
                         by = c("from" = "patch_id")) %>%
        dplyr::left_join(x$df_patch,
                         by = c("to" = "patch_id")) %>%
        dplyr::rowwise() %>%
        dplyr::summarise(n_patch_upstream = min(.data$n_patch_upstream.x,
                                                .data$n_patch_upstream.y)) %>%
        dplyr::pull()

      z <- ifelse(pattern == "downstream", 1, 0)
      prob <- z * v_wa + (1 - z) * (1 / v_wa)

      barrier <- resample(1:n_edge,
                          size = n_barrier,
                          prob = prob)

    } else {

      stop("pattern must be either 'random', 'upstream', or 'downstream'")

    }
  }

  # cumulative fragmentation effect -----------------------------------------

  if (!(length(p) == 1 || length(p) == n_barrier))
    stop(paste("invalid length in p; length must be one or match n_barrier,",
               n_barrier))

  v_p <- rep(1, n_edge)
  v_p[barrier] <- p

  ## negative log passability; igraph does not allow negative values for E()$weight
  igraph::E(g0)$weight <- -log(v_p)
  m_frgm <- exp(-igraph::distances(g0))


  # append frgm matrix ------------------------------------------------------

  if (inherits(x, what = "brnet")) {

    x$df_edge <- dplyr::tibble(patch_x = seq_len(n_patch)[1:(n_patch - 1)],
                               patch_y = seq_len(n_patch)[2:n_patch],
                               passability = v_p)
    x$frgm_matrix <- m_frgm
    x$dispersal_matrix_frgm <- m_disp * m_frgm

    return(x)

  } else {

    if (is.null(m_disp)) stop("dispersal matrix must be provided")
    m_disp_frgm <- m_disp * m_frgm

    return(list(df_edge = dplyr::tibble(patch_x = seq_len(n_patch)[1:(n_patch - 1)],
                                        patch_y = seq_len(n_patch)[2:n_patch],
                                        passability = v_p),
                frgm_matrix = m_frgm,
                dispersal_matrix_frgm = m_disp_frgm))

  }

}

#' Convert an adjacency matrix to a distance matrix
#'
#' @param x Adjacency matrix to be converted
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

adjtodist <- function(x) {

  if (!is.matrix(x)) {
    stop("Adjacency matrix must be provided as a matrix")
  }

  m_adj <- x

  m_distance <- matrix(0,
                       nrow = nrow(x),
                       ncol = ncol(x))

  m_identity <- diag(x = 1,
                     nrow = nrow(x),
                     ncol = ncol(x))

  for (i in seq_len(nrow(x))) {

    m_identity <- m_identity %*% m_adj
    m_distance[m_identity != 0 & m_distance == 0] <- i

    if (length(which(m_distance == 0)) == 0) break

  }

  diag(m_distance) <- 0

  return(m_distance)

}

#' Visualize a branching network
#'
#' @param x `brnet()` object
#' @param patch_color Type of patch (vertex) label (either \code{"env"}, \code{"disturbance"}, \code{"other"} or any color). Default \code{"env"}.
#' @param edge_color Edge color
#' @param edge_weight Type of edge weight (\code{"passability"}). Default \code{NULL}.
#' @param value_col Patch values. Must be specified if \code{patch_color = "other"}.
#' @param color_label Color legend title
#' @param patch_label Type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{"none"}, no label will be shown on patches in the plot. Default \code{"none"}.
#' @param patch_size Patch (vertex) size in the plot.
#' @param ... Arguments passed to \code{ggraph::geom_node_label()}.
#'
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

ggbrnet <- function(x,
                    patch_color = "env",
                    edge_color = "black",
                    edge_weight = "none",
                    value_col = NULL,
                    color_label = NULL,
                    patch_label = "none",
                    patch_size = 3,
                    ...) {

  if (!inherits(x, what = "brnet")) stop("x must be a 'brnet' object")

  # patch attributes --------------------------------------------------------

  adj <- igraph::graph.adjacency(x$adjacency_matrix,
                                 mode = "undirected")

  patch_attr <- x$df_patch

  ## patch color
  if (patch_color == "env") {

    igraph::V(adj)$patch_value <- patch_attr$environment
    if (is.null(color_label)) color_label <- "Environment"

  }

  if (patch_color == "disturb") {

    igraph::V(adj)$patch_value <- patch_attr$disturbance
    if (is.null(color_label)) color_label <- "Disturbance"

  }

  if (patch_color == "other") {

    if (is.null(value_col)) stop("Provide 'value_col' if patch_color = 'other'")
    igraph::V(adj)$patch_value <- unlist(patch_attr[, value_col])
    if (is.null(color_label)) color_label <- "Value"

  }

  ## edge color
  if (edge_weight == "passability") {

    edge_attr <- x$df_edge
    igraph::E(adj)$p <- edge_attr$passability

  }

  ## patch label
  if (patch_label == "patch") {

    igraph::V(adj)$patch_id <- patch_attr$patch_id

  } else {

    if (patch_label == "branch") {

      igraph::V(adj)$patch_id <- patch_attr$branch_id

    } else {

      if (patch_label == "n_upstream") {

        igraph::V(adj)$patch_id <- patch_attr$n_patch_upstream

      } else {

        if (patch_label != "none") message("patch_label should be 'patch', 'branch', or 'n_upstream' to display")
        igraph::V(adj)$patch_id <- ""

      }

    }

  }


  # plot --------------------------------------------------------------------

  g <- ggraph::ggraph(adj,
                      layout = igraph::layout_as_tree(adj,
                                                      root = 1,
                                                      flip.y = FALSE))

  if (!(patch_color %in% c("env", "disturb", "other"))) {
    ## single patch color
    g <- g +
      ggraph::geom_edge_link(color = grey(0.5)) +
      ggraph::geom_node_point(size = patch_size,
                              color = patch_color) +
      ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                              fill = NA,
                              size = 3,
                              label.size = 0,
                              ...) +
      ggplot2::theme_void()

  } else {

    if (!(edge_weight %in% c("passability"))) {
      ## variable patch color with single edge color
      g <- g +
        ggraph::geom_edge_link(color = edge_color) +
        ggraph::geom_node_point(ggplot2::aes(color = .data$patch_value),
                                size = patch_size) +
        ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                                fill = NA,
                                size = 3,
                                label.size = 0,
                                ...) +
        MetBrewer::scale_color_met_c("Hiroshige",
                                     direction = -1) +
        ggplot2::labs(color = color_label) +
        ggplot2::theme_void()

    } else {
      ## variable patch & edge color
      g <- g +
        ggraph::geom_edge_link(ggplot2::aes(alpha = .data$p,
                                            color = edge_color)) +
        ggraph::geom_node_point(ggplot2::aes(color = .data$patch_value),
                                size = patch_size) +
        ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                                fill = NA,
                                size = 3,
                                label.size = 0,
                                ...) +
        MetBrewer::scale_color_met_c("Hiroshige",
                                     direction = -1) +
        ggplot2::labs(color = color_label) +
        ggplot2::theme_void()

    }

  }

  return(g)

}

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
#' @param impact Concentration of environmental pollutants.
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
                  impact = 0,
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
  if (p_disturb > 1 || p_disturb < 0) stop("p_disturb must be 0 to 1")
  if (any(i_disturb > 1) || any(i_disturb < 0)) stop("i_disturb must be 0 to 1")

  v_i_disturb <- fun_to_v(x = i_disturb,
                          n = n_patch)

  ## local pollutants ####
  if (zeta < 0) stop("zeta must be > 0")

  v_zeta <- fun_to_v(x = zeta,
                     n = n_species)

  v_imp <- fun_to_v(x = impact,
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

  } else {

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

    m_r_xt <- m_r_xt * exp(-outer(v_zeta, v_imp))

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

  if (gamma_div == 0 || alpha_div == 0) {
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
#' @param s Strength of switching between basal and ig-prey (confined to 0 - 1). Switching to basal species is more likely to occur with higher values.
#' @param propagule_interval Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param propagule_seed Propagule mean density (intensity parameter in a Poisson distribution). Should be given as a scalar or vector. If given as a vector, the elements should appear in order of basal, intra-guild prey, and intra-guild predator.
#' @param carrying_capacity Carrying capacities of individual patches. Length must be one or equal to \code{n_patch}. Default \code{100}.
#' @param xy_coord Data frame for site coordinates. Each row should correspond to an individual patch, with x and y coordinates (columns). Defualt \code{NULL}.
#' @param distance_matrix Distance matrix indicating distance between habitat patches. If provided, the distance matrix will be used to generate dispersal matrix and to calculate distance decay of environmental correlations. Default \code{NULL}.
#' @param dispersal_matrix Dispersal matrix to be used to simulate dispersal process. Override distance_matrix. Default \code{NULL}.
#' @param p_disturb Disturbance probability.
#' @param i_disturb Disturbance-induced proportional mortality.
#' @param phi_disturb Temporal precision of disturbance-induced proportional mortality. Set \code{Inf} to assume no temporal variability.
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

igpsim <- function(n_patch = 5,
                   n_warmup = 200,
                   n_burnin = 200,
                   n_timestep = 1000,
                   r_b = 10,
                   conv_eff = 0.9,
                   attack_rate = 0.05,
                   handling_time = 0.5,
                   s = 0,
                   propagule_interval = NULL,
                   propagule_seed = c(100, 10, 1),
                   carrying_capacity = 100,
                   xy_coord = NULL,
                   distance_matrix = NULL,
                   dispersal_matrix = NULL,
                   p_disturb = 0,
                   i_disturb = 0,
                   phi_disturb = Inf,
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

  ## propagule ####
  v_seed <- fun_to_v(x = propagule_seed,
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
  if (p_disturb > 1 || p_disturb < 0) stop("p_disturb must be 0 to 1")
  if (any(i_disturb > 1) || any(i_disturb < 0)) stop("i_disturb must be 0 to 1")
  if (phi_disturb <= 0) stop("phi_disturb must be positive")

  v_disturb <- fun_to_v(x = i_disturb,
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
               "disturbance",
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

  ## spatio-temporal variation in disturbance
  shape1 <- phi_disturb * v_disturb
  shape2 <- phi_disturb * (1 - v_disturb)

  if (phi_disturb == Inf) {
    t_disturb <- matrix(v_disturb,
                        nrow = n_patch,
                        ncol = n_sim)
  } else {
    t_disturb <- matrix(stats::rbeta(n_sim * n_patch,
                                     shape1,
                                     shape2),
                        nrow = n_patch,
                        ncol = n_sim)
  }

  ## initial values
  m_n <- matrix(rpois(n = n_species * n_patch,
                      lambda = v_seed),
                nrow = n_species,
                ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

  ## temporal process ####
  for (n in seq_len(n_sim)) {

    if (n_warmup > 0) {

      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n = n_species * n_patch,
                                  lambda = v_seed),
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
                          s = s)

    ## disturbance
    m_n_hat <- t(t(list_n_hat$m_n_hat) * (1 - psi[n] * t_disturb[, n]))

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

      m_dynamics[row_id, ] <- cbind(
        # timestep
        I(n - n_discard),
        # patch id
        rep(x = seq_len(n_patch), each = n_species),
        # carrying capacity
        c(v_k),
        # disturbance
        rep(x = psi[n] * t_disturb[, n], each = n_species),
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

#' Spatial generalized Lotka-Volterra model
#'
#' @param n_species Integer. Number of species
#' @param n_patch Integer. Number of habitat patches
#' @param n_timestep Integer. Number of time step for a simulation run
#' @param interval Numeric. Time interval used for a numerical ODE solver
#' @param r Numeric. Intrinsic growth rates of modeled species
#' @param alpha n_species x n_species interaction matrix
#' @param dispersal List. This list must contain the following parameters as named elements: \code{adj} n_patch x n_patch adjacency matrix; \code{phi} dispersal rate; \code{m} dispersal mortality rate.
#' @param disturb List. This list must contain the following parameters as named elements: \code{int} disturbance intensity; \code{rate} disturbance frequency; \code{s} scale parameter that controls the duration of disturbance.
#' @param n0 List. Initial densities for modeled species that are randomly generated by \code{runif(n_species, min = min, max = max)}. This list must contain \code{min} and \code{max} as named elements.
#' @param threshold Numeric. Extinction threshold. Species will be removed from a simulation if species density goes below this value
#' @param cpp Logical. If \code{TRUE}, Rcpp function is used for the numerical solver
#' @param ... Additional arguments passed to \code{\link{ode}}
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

sglv <- function(n_species,
                 n_patch,
                 n_timestep = 100,
                 interval = 0.01,
                 r,
                 alpha,
                 dispersal = list(adj = matrix(0,
                                               nrow = n_patch,
                                               ncol = n_patch),
                                  phi = 0,
                                  m = 0),
                 disturb = list(int = 1,
                                rate = (1 / n_timestep) * 10,
                                s = interval * 10),
                 n0 = list(min = 0,
                           max = 1),
                 threshold = 1E-4,
                 cpp = TRUE,
                 ...) {

  # verify inputs -----------------------------------------------------------

  # r input
  if (!inherits(r, "matrix")) {
    if (length(r) != 1 && length(r) != n_species)
      stop("r must have length one or n_species")
  }

  # alpha input
  if (!all(dim(alpha) == n_species))
    stop("alpha's dimensions must be n_species x n_species")

  # dispersal input
  if (!all(sort(names(dispersal)) == sort(c("adj", "phi", "m"))))
    stop("Elements in the dispersal list must be 'adj', 'phi', and 'm'")

  with(dispersal, {
    if (!all(dim(adj) == n_patch))
      stop("adj's dimensions must be n_patch x n_patch")

    if (phi < 0 || m < 0)
      stop("'phi' and 'm' must be >= 0")
  })

  # disturbance input
  if (!all(sort(names(disturb)) == sort(c("int", "rate", "s"))))
    stop("Elements in the dispersal list must be 'int', 'rate', and 's'")

  with(disturb, {
    if (length(int) != 1 && length(int) != n_patch)
      stop("disturbance intensity (int) must have length one or n_patch")

    if (int < 0 || rate < 0 || s < 0)
      stop("'int,' 'rate,' and 's' must be >= 0")
  })

  # n0 input
  if (!all(sort(names(n0)) == sort(c("min", "max"))))
    stop("Elements in the n0 list must be 'min' and 'max'")

  with(n0, {
    if (min < 0 || max < 0 || min > max)
      stop("invalid inputs in n0")
  })

  # disturbance setup -------------------------------------------------------

  # n possible disturbance points
  # + 100 to ensure cumsum(x) > n_timestep
  # remove points > n_timestep
  psi <- with(disturb, cumsum(rexp(n = ceiling((1 / rate) + 100), rate)))
  psi <- psi[psi <= n_timestep]

  # pseudo time steps
  time <- seq(0, n_timestep, by = disturb$s)
  t_psi <- sapply(psi, FUN = function(x) which.min(abs(time - x)))

  # create dummy time-series for signals
  signal <- data.frame(time = time,
                       psi = rep(0, length(time)))

  signal$psi[t_psi] <- 1

  # linear interpolation function
  input <- stats::approxfun(signal, rule = 2)

  # n_patch x n_species disturbance intensity matrix
  ## E: fun_to_m() returns n_species x n_patch matrix
  ## v_e: E into vector (patch1; sp1, sp2, ...) x ... (patch n; sp1, sp2, ...)
  E <- with(disturb, fun_to_m(int,
                              n_species = n_species,
                              n_patch = n_patch,
                              param_attr = "patch"))

  v_e <- c(E$m_x)

  # intrinsic growth --------------------------------------------------------

  ## v_r: vector (patch1; sp1, sp2, ...) ... (patch n; sp1, sp2, ...)
  ## length(v_r) = n_species x n_patch
  if (inherits(r, "matrix")) {

    if (!all(dim(r) == c(n_species, n_patch)))
      stop("Dimension mismatch in r; r must be a n_species x n_patch matrix")

    v_r <- c(r)

  } else {

    R <- fun_to_m(r,
                  n_species = n_species,
                  n_patch = n_patch,
                  param_attr = "species")

    v_r <- c(R$m_x)

  }


  # interaction -------------------------------------------------------------

  ## large interaction matrix
  ## dim: (n_species + n_patch) x (n_species + n_patch)
  A <- kronecker(diag(n_patch), alpha)

  # dispersal ---------------------------------------------------------------

  # pick id of non-zero columns
  C1 <- C0 <- with(dispersal, adj)
  id_nz <- which(colSums(C0) != 0)

  # scale adjacency matrix
  if (length(id_nz) > 0) {

    if (length(id_nz) == 1) {
      x <- sum(C0[, id_nz])
    } else {
      x <- colSums(C0[, id_nz])
    }

    C1[, id_nz] <- t(t(C0[, id_nz]) / x)
  }

  # apply dispersal rate and mortality
  ## (phi - m) * diag(n_species): immigration, disp. mortality m
  ## -phi * diag(n_species): emigration
  m_off <- with(dispersal, (phi - m) * diag(n_species))
  m_diag <- with(dispersal, -phi * diag(n_species))

  ## expand to large connectivity matrix
  ## dim: (n_species + n_patch) x (n_species + n_patch)
  C_off <- kronecker(C1, m_off)
  C_diag <- kronecker(diag(n_patch), m_diag)
  C <- C_diag + C_off

  # ode ---------------------------------------------------------------------

  # n: n_patch x n_species abundance vector
  # N: n_patch x n_species abundance matrix
  # R: n_patch x n_species growth rate matrix
  # psi: indicator parameter scalar
  # E: n_patch x n_species disturbance intensity matrix
  # A: n_species x n_species interaction matrix
  # phi: migration rate
  # C: n_patch x n_patch connectivity matrix

  derivr <- function(t, n, parms) {
    with(parms, {
      psi <- input(t)
      dn <- n * (r - psi * e + t(A) %*% n) + C %*% n
      list(dn)
    })
  }

  # define absorbing condition
  ## root function
  rootfun <- function(t, n, parms) {
    return(n - threshold)
  }

  ## extinction: triggered when "n - threshold = 0"
  eventfun <- function(t, n, parms) {
    n <- ifelse(n <= threshold, 0, n)
    return(n)
  }

  # parameter list
  parms <- list(r = v_r,
                e = v_e,
                A = A,
                C = C)

  # initial values for a state variable
  n_init <- with(n0, stats::runif(n_species * n_patch,
                                  min = min,
                                  max = max))

  # time-series
  times <- seq(0, n_timestep, by = interval)

  # run ode solver
  cout <- deSolve::ode(y = n_init,
                       times = times,
                       func = ifelse(cpp, deriv, derivr),
                       parms = parms,
                       events = list(func = eventfun, root = TRUE),
                       rootfun = rootfun,
                       ...)

  return(cout)
}

#' Find suitable intrinsic growth rates
#'
#' @param alpha Numeric. n_species x n_species interaction matrix
#' @param k0 Numeric. Total carrying capacity for basal species combined
#' @param lambda0 Numeric. Initial lambda value for the exponential decay of equilibrium densities with trophic position
#' @param interval Numeric. Increment of lambda value.
#' @param sigma Numeric. Degree of noise added to equilibrium densities
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

findr <- function(alpha,
                  k0,
                  lambda0 = 1E-3,
                  interval = 0.01,
                  sigma = 0.01) {

  # half interaction matrix
  alpha0 <- alpha
  alpha0[lower.tri(alpha0, diag = TRUE)] <- 0
  alpha0[alpha0 != 0] <- 1

  # basal species id
  id_basal <- which(colSums(alpha0) == 0)
  n_basal <- length(id_basal)
  n_c <- ncol(alpha) - n_basal

  # trophic position
  v_k <- runif(length(id_basal))
  f_k <- k0 * (v_k / sum(v_k))
  tp <- attr(alpha, "tp")

  ## k0 varies by basal species
  ## for basal species, use f_k
  ## for consumers, use mean k0 (`mean(f_k)`)
  w_k0 <- rep(mean(f_k), ncol(alpha))
  w_k0[1:n_basal] <- f_k

  # initialize lambda and r
  lambda <- lambda0
  eps <- exp(rnorm(n_c, mean = 0, sd = sigma))

  x <- w_k0 * exp(-lambda * (tp - 1))
  x[-id_basal] <- x[-id_basal] * eps
  r <- drop(- x %*% alpha)

  # loop until all consumer's r < 0
  while(any(r[-id_basal] >= 0)) {
    lambda <- lambda + interval
    x <- w_k0 * exp(-lambda * (tp - 1))
    x[-id_basal] <- x[-id_basal] * eps
    r <- drop(- x %*% alpha)
  }

  cout <- cbind(r, x, tp)
  rownames(cout) <- 1:ncol(alpha)
  colnames(cout) <- c("r", "equilibrium", "tp")
  attr(cout, "lambda") <- lambda

  return(cout)
}

#' Preferential prey model
#'
#' @param n_species Integer. Number of species
#' @param n_basal Integer. Number of basal species
#' @param l Interger. Expected number of links in the upper triangle
#' @param theta Numeric. Scale parameter of an exponential distribution. Smaller values indicate greater trophic specialization.
#' @param cannibal Logical. If \code{TRUE}, cannibalism allowed
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

ppm <- function(n_species,
                n_basal,
                l,
                theta,
                cannibal = FALSE) {
  # n_species: number of species in the pool
  # n_basal: number of basal species
  # l: expected number of links
  # theta: scale parameter

  # verify inputs
  if (n_basal < 1) stop("At least one basal species needed to construct a food web")
  if (n_basal >= n_species) stop("n_basal must be smaller than n_species")

  # number of consumers
  n_c <- n_species - n_basal

  # tp: trophic position
  # assign trophic position for basal species
  # assign -1 for consumers as initial values (to be updated)
  tp <- rep(-1, n_species)
  tp[seq_len(n_basal)] <- 1

  # beta parameter for beta distribution
  # determined so that E(L) = L
  if (l < n_c)
    stop("l must be at least equal to the number of consumers (n_species - n_basal)")

  if (cannibal) {
    max_l <- sum((n_basal + 1):n_species)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 1)
      stop("n_species + n_basal must be equal to or greater than 1")

    b <- ((n_species + n_basal - 1) * n_c) / (2 * (l - n_c)) - 1
  } else {
    max_l <- sum((n_basal + 1):n_species - 1)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 3)
      stop("n_species + n_basal must be equal to or greater than 3")

    b <- ((n_species + n_basal - 3) * n_c) / (2 * (l - n_c)) - 1
  }

  # vector for link proportions for all consumers
  v_xi <- stats::rbeta(n_c,
                       shape1 = 1,
                       shape2 = b)

  # vector for initial prey choice for all consumers
  v_i0 <- c(rep(-1, n_basal),
            sapply(X = (n_basal + 1):n_species,
                   FUN = function(j) resample(seq_len(j - 1), size = 1)))

  # realized number of prey nodes for all consumers
  # kappa follows a beta-binomial distribution
  # size parameter is the possible maximum number of prey nodes
  if (cannibal) {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = n_basal:(n_species - 1),
                               prob = v_xi))
  } else {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = (n_basal - 1):(n_species - 2),
                               prob = v_xi))
  }

  # alpha0: S x S interaction binary matrix
  # initialized with all zero
  # update the initial consumer's first prey
  alpha0 <- matrix(0, n_species, n_species)
  alpha0[v_i0[n_basal + 1], n_basal + 1] <- 1

  for (j in (n_basal + 1):n_species) {
    alpha0 <- extra_prey(alpha0 = alpha0,
                         j = j,
                         i0 = v_i0[j],
                         tp = tp,
                         kappa = v_kappa[j],
                         theta = theta,
                         cannibal = cannibal)

    v_n_prey <- colSums(alpha0)
    v_n_prey[v_n_prey == 0] <- 1

    tp_new <- (tp %*% alpha0) / v_n_prey + 1
    tp[j] <- tp_new[j]
  }

  attr(alpha0, "tp") <- tp

  return(alpha0)
}

#' Extra prey function
#'
#' @inheritParams ppm
#' @param alpha0 Interaction matrix
#' @param j Integer. Consumer's index
#' @param i0 Integer. Index for the first prey
#' @param tp Numeric. Initial trophic positions
#' @param kappa Integer. Number of extra prey items
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

extra_prey <- function(alpha0,
                       j,
                       i0,
                       tp,
                       theta,
                       kappa,
                       cannibal = FALSE) {
  # alpha0:      adjacency matrix defining trophic interactions
  # j:      consumer's index (j - 1 is the number of possible prey)
  # i0:     index of the first prey chosen (1 =< i0 < j)
  # tp:     vector of trophic positions for possible prey nodes (length = S)
  # kappa:  number of extra prey nodes in addition to the first prey i0
  # theta:  scale parameter for an exponential decay of prey preference

  # verify inputs
  if (length(tp) != ncol(alpha0)) stop(paste("'tp' length seems incorrect;",
                                             "the length must be",
                                             ncol(alpha0)))

  if (!(i0 < j && 1 <= i0)) stop("i0 must be a non-zero intger smaller than j")

  if (kappa >= j) stop("kappa must be an integer smaller than j")

  # probability of consumer j picks prey i
  lambda <- 1 / theta
  tp[j] <- tp[i0] + 1
  p_ij <- exp(- lambda * abs(tp[i0] - tp))

  # pick kappa prey species out of j - 1 nodes (without cannibalism) or j nodes
  if (cannibal) {
    i_index <- 1:j
    p_ij <- p_ij[1:j]
  } else {
    i_index <- 1:(j - 1)
    p_ij <- p_ij[1:(j - 1)]
  }

  # exclude i0 from resample() because i0 is the first pick
  # pick 'kappa' samples from index[-i0]
  # weighted by p_ij
  if (length(i_index) > 1) {
    i_pick <- resample(i_index[-i0],
                       size = kappa,
                       prob = p_ij[-i0])

    i <- c(i0, i_pick)
  } else {
    # when j = 2 with no cannibalism
    # no choice but i0 available as prey
    i <- i0
  }

  alpha0[i, j] <- 1

  return(alpha0)
}

#' Apply conversion efficiency and attack rate
#'
#' @inheritParams extra_prey
#' @param attack List for attack rates. Specify minimum and maximum values for a uniform distribution.
#' @param convert List for conversion efficiency. Specify minimum and maximum values for a uniform distribution.
#' @param mortal List for mortality (or intraspecific competition). Specify minimum and maximum values for a uniform distribution.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

to_alpha <- function(alpha0,
                     attack = list(min = 0,
                                   max = 1),
                     convert = list(min = 0,
                                    max = 1),
                     mortal = list(min = 0,
                                   max = 1)) {

  # identify basal species
  u_alpha0 <- alpha0
  basal <- which(colSums(alpha0) == 0)

  # scale by # of resources
  su_alpha0 <- t(t(alpha0[, -basal]) / colSums(alpha0)[-basal])

  # generate random parameters
  ## matrix for attack rates
  a <- with(attack, matrix(stats::runif(ncol(alpha0)^2, min = min, max = max),
                           nrow = nrow(alpha0),
                           ncol = ncol(alpha0)))

  ## matrix for conversion efficiency
  b <- with(convert, matrix(stats::runif(ncol(alpha0)^2, min = min, max = max),
                            nrow = nrow(alpha0),
                            ncol = ncol(alpha0)))

  ## vector for intraspecific competition
  m <- with(mortal, stats::runif(ncol(alpha0), min = min, max = max))

  # interaction matrix
  # diag(alpha) represent net effects of cannibalism
  u_alpha0[, -basal] <- a[, -basal] * su_alpha0
  alpha <- b * u_alpha0 - t(u_alpha0)

  # intraspecific interaction
  # added to net effects of cannibalism as "diag(alpha) - m"
  diag(alpha) <- diag(alpha) - m

  return(alpha)
}

#' Generate a food web based on the preferential prey model
#'
#' @inheritParams ppm
#' @inheritParams to_alpha
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodweb <- function(n_species,
                    n_basal,
                    l,
                    theta,
                    cannibal = FALSE,
                    attack = list(min = 0,
                                  max = 1),
                    convert = list(min = 0,
                                   max = 1),
                    mortal = list(min = 0,
                                  max = 1)) {

  alpha0 <- ppm(n_species = n_species,
                n_basal = n_basal,
                theta = theta,
                l = l,
                cannibal = cannibal)

  alpha <- to_alpha(alpha0 = alpha0,
                    attack = attack,
                    convert = convert,
                    mortal = mortal)

  return(alpha)
}

#' Maximum trophic position
#'
#' @param n_species Integer. Number of species
#' @param n Numeric. Vector of species density
#' @param alpha Numeric. Interaction matrix
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

max_tp <- function(n_species, n, alpha) {

  alpha[lower.tri(alpha, diag = TRUE)] <- 0

  # initialize trophic positions for basal species
  tp <- rep(-1, n_species)

  id_basal <- which(colSums(alpha) == 0)
  n_basal <- length(id_basal)
  tp[id_basal] <- 1

  if (all(n[id_basal] == 0)) {

    # if basal species are all absent
    return(0)

  } else {

    # total gain from prey (prey converted to predator density)
    g <- drop(n %*% alpha)

    # fractional contribution of each prey
    ## frac excludes basal species
    frac <- do.call(rbind,
                    lapply(1:n_species,
                           function(i) {
                             ## num: contribution of each prey
                             ## den: total prey converted into consumer
                             ## frac0: factional contribution of each prey
                             num <- alpha[i, (n_basal + 1):n_species] * n[i]
                             den <- g[(n_basal + 1):n_species]
                             frac0 <- num / den

                             ## 0 / 0 returns NaN
                             ## replace with 0 if num == 0 && den == 0
                             k <- intersect(which(num == 0), which(den == 0))
                             frac0[k] <- 0

                             return(frac0)
                           }))

    ## add columns for basal species (all zero)
    frac <- cbind(matrix(0, nrow = nrow(frac), ncol = n_basal), frac)

    # update trophic positions recursively
    for (i in (n_basal + 1):n_species) {
      value <- tp %*% frac + 1
      tp[i] <- ifelse(g[i] == 0, 0, value[i])
    }

    max_tp <- max(tp[n > 0])

    return(max_tp)
  }

}

#' Spatiotemporal average of maximum trophic positions using \code{\link{sglv}} output
#'
#' @param n Numeric. Matrix output from \code{\link{sglv}}
#' @param n_species Integer. Number of species
#' @param n_patch Integer. Number of habitat patches
#' @param alpha Numeric. Interaction matrix
#' @param start Integer. Initial time step to be included in the calculation
#' @param end Integer. Last time step to be included in the calculation
#' @param full Logical. If \code{TRUE}, a full spatiotemporal matrix of food chain length returned as attributes.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodchain <- function(n,
                      n_species,
                      n_patch,
                      alpha,
                      start = max(n[, 1]) - 99,
                      end = max(n[, 1]),
                      full = FALSE) {

  # verify input ------------------------------------------------------------
  if (!inherits(n, "deSolve"))
    stop("n must be class 'deSolve'")

  if (ncol(n) != n_species * n_patch + 1)
    stop(paste("Dimension mismatch; n_species x n_patch must be",
               ncol(n) - 1))

  if (!all(dim(alpha) == n_species))
    stop(paste("Dimension mismatch; alpha's dimentions must be",
               n_species,
               "x",
               n_species))

  # subset data -------------------------------------------------------------
  times <- seq(start, end, by = 1)
  index_t <- which(n[, "time"] %in% times)
  x <- n[index_t, -which(colnames(n) == "time")]

  # tidy format -------------------------------------------------------------
  ## time vector
  v_t <- rep(times, times = n_species * n_patch)

  ## species vector
  v_sp <- rep(seq_len(n_species), each = length(times) * n_patch)

  ## patch vector
  v_patch <- rep(rep(seq_len(n_patch), each = length(times)), times = n_species)

  ## combine: g specifies a unique combination of time and patch
  y <- data.table::data.table(t = v_t,
                              species = v_sp,
                              patch = v_patch,
                              g = paste0("t",
                                         sprintf("%05d", v_t),
                                         "p",
                                         sprintf("%05d", v_patch)),
                              n = c(x))

  # food chain length -------------------------------------------------------
  ## vector - time (t) and patch-specific (j) food chains
  v_fcl <- with(y, tapply(y, INDEX = g,
                          FUN = function(d) {
                            with(d, max_tp(n_species = n_species,
                                           n = n,
                                           alpha = alpha)
                            )
                          }))

  ## matrix - time (t) and patch-specific (j) food chains
  ## row: time, col: patch
  m_fcl <- matrix(v_fcl,
                  nrow = length(times),
                  ncol = n_patch,
                  byrow = TRUE)

  colnames(m_fcl) <- paste0("patch", seq_len(n_patch))
  rownames(m_fcl) <- times

  fcl <- mean(v_fcl)

  if (full) {
    attr(fcl, "fcl_matrix") <- m_fcl
  }

  return(fcl)
}

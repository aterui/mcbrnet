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
#' @param plot Whether a plot should be shown or not. If \code{FALSE}, a plot of the generated network will not be shown. Default \code{TRUE}.
#' @param patch_color Type of patch (vertex) label (either \code{"env"}, \code{"disturbance"} or any color code). Default \code{"env"}.
#' @param patch_label Type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{"none"}, no label will be shown on patches in the plot. Default \code{"none"}.
#' @param patch_size Patch (vertex) size in the plot.
#' @param n_patch_free Whether imposing a constraint on \code{n_patch}. If \code{TRUE}, the number of patches a random variable following a negative binomial distribution.
#' @param ... Arguments passed to \code{ggraph::geom_node_label()}.
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
                  plot = FALSE,
                  patch_color = "env",
                  patch_label = "none",
                  patch_size = 3,
                  n_patch_free = FALSE,
                  ...) {


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


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {

    g <- ggbrnet(x = list(adjacency_matrix = m_adj,
                          distance_matrix = m_distance,
                          df_patch = df_patch),
                 patch_color = patch_color,
                 patch_label = patch_label,
                 patch_size = patch_size,
                 ...)

    print(g)

  }

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

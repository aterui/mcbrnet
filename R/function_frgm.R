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
#'

frgm <- function(x,
                 rate,
                 pattern = "random",
                 p = NULL,
                 n_barrier) {


  # adjacency graph & distance matrix ---------------------------------------

  if(inherits(x, what = "brnet")) {
    m_adj <- x$adjacency_matrix
  } else {
    m_adj <- x
  }

  if(dplyr::n_distinct(dim(m_adj)) != 1) stop("invalid dimension in the adjacency matrix")

  g0 <- m_adj %>%
    igraph::graph.adjacency("undirected")

  df_g0 <- igraph::as_data_frame(g0)

  if(inherits(x, what = "brnet")) {
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

    if (pattern == "downstream"|pattern == "upstream") {

      if (!inherits(x, what = "brnet")) stop("x must be class 'brnet'")

      v_wa <- df_g0 %>%
        dplyr::as_tibble() %>%
        dplyr::left_join(x$df_patch,
                         by = c("from" = "patch_id")) %>%
        dplyr::left_join(x$df_patch,
                         by = c("to" = "patch_id")) %>%
        dplyr::rowwise() %>%
        dplyr::summarise(n_patch_upstream = min(n_patch_upstream.x,
                                                n_patch_upstream.y)) %>%
        dplyr::pull()

      z <- ifelse(pattern == "downstream", 1, 0)
      prob <- z * v_wa + (1 - z) * (1 / v_wa)

      barrier <- sample(n_edge,
                        size = n_barrier,
                        prob = prob)

    } else {

      stop("pattern must be either 'random', 'upstream', or 'downstream'")

    }
  }

  # cumulative fragmentation effect -----------------------------------------

  if (!(length(p) == 1 | length(p) == n_barrier)) stop(paste("invalid length in p;
                                                       length must be one or
                                                       match n_barrier,",
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

    return(list(df_edge <- dplyr::tibble(patch_x = seq_len(n_patch)[1:(n_patch - 1)],
                                         patch_y = seq_len(n_patch)[2:n_patch],
                                         passability = v_p),
                frgm_matrix = m_frgm,
                dispersal_matrix_frgm = m_disp_frgm))

  }

}

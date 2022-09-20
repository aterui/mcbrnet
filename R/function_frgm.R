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
    igraph::graph_from_adjacency_matrix("undirected")

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
  }

  if (pattern == "cluster") {
    s <- resample(seq_len(n_edge), size = 1)
    barrier <- order(m_dist[s, ])[seq_len(n_barrier)]
  }

  if (pattern == "downstream"|pattern == "upstream") {

    if (!inherits(x, what = "brnet")) stop("x must be class 'brnet'")

    v_wa <- x$df_patch$n_patch_upstream
    s <- ifelse(pattern == "downstream",
                order(-v_wa),
                order(v_wa))

    barrier <- s[seq_len(n_barrier)]

  }


  # cumulative fragmentation effect -----------------------------------------

  if (!(length(p) == 1 | length(p) == n_barrier)) stop(paste("invalid length in p;
                                                       length must be one or
                                                       match n_barrier,",
                                                       n_barrier))

  v_p <- rep(1, n_edge)
  v_p[barrier] <- p

  igraph::E(g0)$weight <- -log(v_p)
  m_frgm <- exp(-igraph::distances(g0))


  # append frgm matrix ------------------------------------------------------

  if (inherits(x, what = "brnet")) {

    x$frgm_matrix <- m_frgm
    x$dispersal_matrix_frgm <- m_disp * m_frgm

    return(x)

  } else {

    if (is.null(m_disp)) stop("dispersal matrix must be provided")
    m_disp_frgm <- m_disp * m_frgm

    return(list(frgm_matrix = m_frgm,
                dispersal_matrix_frgm = m_disp_frgm))

  }

}

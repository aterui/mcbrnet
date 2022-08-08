#' Return a dispersal probability matrix after accounting for network fragmentation
#'
#' @param x
#' @param dispersal_matrix adjacency matrix
#' @param pattern
#' @param p
#' @param n_barrier
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
                 dispersal_matrix,
                 patten = "random",
                 p = NULL,
                 n_barrier) {

  ## adjacency graph
  m_adj <- ifelse(class(x) == "brnet", x$adjacency_matrix, x)
  g0 <- m_adj %>%
    graph_from_adjacency_matrix("undirected")

  ## distance matrix
  m_dist <- ifelse(class(x) == "brnet",
                   x$distance_matrix,
                   igraph::distances(g0))

  ## basic numbers
  n_patch <- unique(dim(m_adj))
  n_edge <- n_patch - 1
  m_disp <- dispersal_matrix

  if(dplyr::n_distinct(dim(m_adj)) != 1) stop("invalid dimension in the adjacency matrix")

  if (!is.null(p)) {

    ## if p specified
    if(any(p < 0 | p > 1)) stop("p must be 0 - 1.")
    if(n_edge != length(p)) stop("invalid length in p;
                                  resistance p must match
                                  the dimension of
                                  the adjacency matrix")

    v_p <- p
    E(g0)$weight <- -log(v_p)
    m_frag <- exp(-distances(g0))

  } else {

    ## if p not specified
    if (n_barrier > n_edge) stop("n_barrier exceeds the number of edges in the graph")

    if (pattern == "random") {
      barrier <- resample(seq_len(n_edge), size = n_barrier)
    }

    if (pattern == "cluster") {
      s <- resample(seq_len(n_edge), size = 1)
      barrier <- order(m_dist[s, ])[seq_len(n_barrier)]
    }

    if (pattern == "downstream"|"upstream") {

      if (class(x) != "brnet") stop("x must be class 'brnet'")

      v_wa <- x$df_patch$n_patch_upstream
      s <- ifelse(pattern == "downstream",
                  order(-v_wa),
                  order(v_wa))

      barrier <- s[seq_len(n_barrier)]

    }

  }

}

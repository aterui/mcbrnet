#' Return a dispersal probability matrix after accounting for network fragmentation
#'
#' @param adjacency_matrix adjacency matrix
#' @param dispersal_matrix adjacency matrix
#' @param omega adjacency matrix
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

frgm <- function(adjacency_matrix,
                 dispersal_matrix
                 omega) {

  m_adj <- adjacency_matrix
  m_disp <- dispersal_matrix

  if(dplyr::n_distinct(dim(m_adj)) != 1) stop("invalid dimension in the adjacency matrix")
  if(unique(dim(m_adj)) != length(omega)) stop('invalid length in "prob"; must match the dimension of the adjacency matrix')

  n_patch <- unique(dim(m_adj))
  v_seq <- 1:n_patch

  m_p <- sapply(1:n_patch, function(i) {
    y <- m_adj %>%
      graph_from_adjacency_matrix("undirected") %>%
      shortest_paths(from = i) %>%
      .$vpath %>%
      lapply(FUN = function(x) cumprod(omega[x[-1]])[length(x) - 1]) %>%
      unlist()

    v_seq[v_seq != i] <- y

    return(v_seq)
  }) %>%
    t()

  diag(m_p) <- 0

  m_disp_prime <- m_p * m_disp

  return(m_disp_prime)
}

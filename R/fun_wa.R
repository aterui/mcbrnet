#' Internal function: watershed area
#'
#' @param x Adjacency matrix to be converted
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_wa <- function(x,
                   n_patch) {

  m_adj_up <- x
  m_adj_up[lower.tri(m_adj_up)] <- 0

  m_wa <- matrix(0,
                 ncol = n_patch,
                 nrow = n_patch)

  m_identity <- diag(1,
                     nrow = n_patch,
                     ncol = n_patch)

  for (i in seq_len(n_patch)) {

    m_identity <- m_identity %*% m_adj_up
    m_wa[m_identity != 0 & m_wa == 0] <- 1

    if (length(which(m_wa[upper.tri(m_wa)] == 0)) == 0) break

  }

  diag(m_wa) <- 1

  return(m_wa)
}


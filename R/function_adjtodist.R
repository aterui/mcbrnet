#' Convert an adjacency matrix to a distance matrix
#'
#' @param x adjacency matrix to be converted
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

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

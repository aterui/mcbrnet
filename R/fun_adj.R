#' Internal function: return adjacency in a linear habitat
#'
#' @param n Number of nodes
#' @param start_id Start node ID
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_adj <- function(n, start_id = 1) {

  if (n == 1) m_y <- cbind(1, NA)

  if (n == 2) m_y <- cbind(c(1, 2), c(2, 1))

  if (n > 2) {

    y1 <- c(1, rep(2:(n - 1), each = 2), n)
    y2 <- c(2, sapply(2:(n - 1), function(i) c(i - 1, i + 1)), n - 1)
    m_y <- cbind(y1, y2)

  }

  colnames(m_y) <- c("patch_1", "patch_2")

  m_y <- m_y + (start_id - 1)

  return(m_y)
}

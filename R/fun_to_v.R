#' Internal function: to vector
#'
#' @param x scalar or vector
#' @param n replication (n_species or n_patch)
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_to_v <- function(x, n) {

  if (!is.vector(x)) stop("x must be a vector")

  if (length(x) == 1) {

    message("single value is given for an species or patch attribute: assume identical across species or patches")
    v_k <- rep(x = x,
               times = n)

  } else {

    if (length(x) != n) stop("x must have a length of one or n_species/n_patch")
    v_k <- x

  }

  return(v_k)

}

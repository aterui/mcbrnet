#' Internal function: to vector
#'
#' @param x scalar or vector
#' @param n_species number of species
#' @param n_patch number of patches
#' @param param_attr indicate species or patch attribute
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_to_m <- function(x,
                     n_species,
                     n_patch,
                     param_attr) {

  # check input -------------------------------------------------------------

  if (!is.vector(x)) stop("x must be a vector")


  # conversion to matrix ----------------------------------------------------

  ## patch-wise
  if (param_attr == "patch") {

    if (length(x) == 1) {

      message("single value is given for an species or patch attribute: assume identical across species or patches")

      v_x <- rep(x = x, times = n_patch)
      m_x <- matrix(x,
                    nrow = n_species,
                    ncol = n_patch)

    } else {

      if (length(x) != n_patch) stop("x must have a length of one or n_patch")

      v_x <- x
      m_x <- matrix(rep(x = x,
                        each = n_species),
                    nrow = n_species,
                    ncol = n_patch)

    } # ifelse

  } # if (param_attr == "patch")

  ## species-wise
  if (param_attr == "species") {

    if (length(x) == 1) {

      message("single value is given for an species or patch attribute: assume identical across species or patches")

      v_x <- rep(x = x, times = n_species)
      m_x <- matrix(x,
                    nrow = n_species,
                    ncol = n_patch)

    } else {

      if (length(x) != n_species) stop("x must have a length of one or n_species")

      v_x <- x
      m_x <- matrix(rep(x = x,
                        times = n_patch),
                    nrow = n_species,
                    ncol = n_patch)

    } # ifelse

  } # if (param_attr == "species")


  # export ------------------------------------------------------------------

  return(list(v_x = v_x,
              m_x = m_x))

}

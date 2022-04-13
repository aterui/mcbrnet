#' Internal function: get n_branch
#'
#' @param n_patch Number of patches
#' @param p_branch Branching probability
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_get_n_branch <- function(n_patch,
                             p_branch) {

  if (p_branch > 0 & p_branch < 1) {

    repeat {

      n_branch <- rbinom(n = 1, size = n_patch, prob = p_branch)
      if (n_branch %% 2 == 1) break

    }

  } else {

    if (p_branch == 0) n_branch <- 1

    if (p_branch == 1) {

      if (n_patch %% 2 == 0) stop("n_patch must be an odd number when p_branch = 1")

      n_branch <- n_patch

    }

  }

  return(n_branch)

}

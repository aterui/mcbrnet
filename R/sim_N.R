#' Simulate abundance matrix, N
#'
#' @description This is a helper function to simulate an abundance matrix N. This is not used in the main package, and is included to trobleshoot problems and for illustrative purposes.
#'
#' @param n_patch Integer, number of patches (columns in N)
#' @param fcl_5_states Logical. If `TRUE`, `n_patch` must be >= 5, forces the first 5 columns of `N` to represent all possible community compositions in order: No species, B only, B+P, B+C, B+P+C.
#'
#' @return `N` abundance matrix with `nrow = 3` and `ncol = n_patch`.
#' @export
#'
#' @examples
#'
#' # simulate N with all 5 community composition possibilities.
#' sim_N(n_patch = 5, fcl_5_states = TRUE)
#'
sim_N <- function(n_patch = 5, fcl_5_states = FALSE){
  n_sp = 3
  if(fcl_5_states == FALSE){
  N = matrix(rpois(n_patch * n_sp,
                   lambda = c(100, 50, 10)),
             nrow = n_sp,
             ncol = n_patch)
  }
  if(fcl_5_states == TRUE){
    if(n_patch < 5) stop("n_patch must be >= 5 for all FCL states to be present")
    N = matrix(rpois(n_patch * n_sp,
                     lambda = c(100, 50, 10)),
               nrow = n_sp,
               ncol = n_patch)
    N[,1] <- 0
    N[2:3,2] <- 0
    N[2,3] <- 0
    N[3,4] <- 0
  }
  N
}

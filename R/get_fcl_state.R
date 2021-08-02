#' Calculates FCL "state" variable
#'
#' @description Internal function used by `igp_sim()` to describe community composition in each patch for each time step
#'
#' @param N Matrix of population abundances, where `nrow = 3`, and `ncol = n_patch`.
#'
#' @details This function calculates the food chain length state of each patch at each time step in the IGP simulation. This is used internally.
#'
#' @return numeric vector of length = `n_patch`. FCL state 0 = no species, 1 = B only, 2 = B + C, 2.5 = B + P, 3 = B + C + P
#' @export
#'
#' @examples
#' # simulate abundance matrix
#' # nrow = n_sp = 3
#' # ncol = n_patch = 10
#' N <- matrix(rpois(3 * 10, lambda = 1), nrow = 3, ncol = 10)
#'
#' # make sure abundance of C and P = 0 when abundance of B = 0
#' N[,N[1,] == 0] <- 0
#'
#' # calculate food chain length state
#' get_fcl(N)
#'
get_fcl_state <- function(N){
  if(!is.matrix(N)){
    stop("get_fcl() requires N to be a matrix")
  }
  N[,N[1,] == 0] <- 0
  patch_fcl = N
  patch_fcl[patch_fcl>0] = 1
  patch_fcl = patch_fcl*c(1, 2, 3)
  fcl = colSums(patch_fcl)
  fcl[fcl == 3] <- 2
  fcl[fcl == 4] <- 2.5
  fcl[fcl == 6] <- 3
  fcl
}

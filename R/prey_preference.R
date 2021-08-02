#' Calculate preference between two prey resources.
#'
#' @description Internal function used by `pop_sim()` when `igp_sim(P_pref = NULL)`
#'
#' @param e1 Conversion efficiency for prey resource 1. Most often this `ebp`, or the efficiency of turning B biomass into new P biomass.
#' @param e2 Conversion efficiency for prey resource 2. Most often this `ecp`, or the efficiency of turning C biomass into new P biomass.
#' @param N1 The population abundance of prey resource 1. Usually this is prey B
#' @param N2 The population abundance of prey resource 2. Usually this is prey C
#'
#'@details The preference of prey resource 1 over prey resource 2. Values < 0.5 indicate the resource 2 is preferred, value > 0.5 indicate resource 1 is preferred.
#'
#'Preference is calculated as `e1*N1 / (e1*N1 + e2*N2)`
#'
#' @return value between 0 and 1 describing search preference for resource B by predator P.
#' @export
#'
#' @examples
#' prey_preference(e1 = 2, e2 = 1, N1 = 100, N2 = 50)
#'
prey_preference <- function(e1, e2, N1, N2){
  if(any(is.na(c(e1, e2, N1, N2)))){
    stop("all inputs to `prey_preference()` need to be defined,\n one or more input values is 'NA' or 'NaN'" )
  }
  if(any(c(e1, e2, N2, N1) < 0)){
    stop("input to 'prey_preference()' should be positive values,\n at least one is negative")
  }
  if(!identical(length(e1), length(e2), length(N1), length(N2))){
    stop("length of inputs to 'prey_preference()' are not all equal")
  }
  e1*N1 / (e1*N1 + e2*N2)
}

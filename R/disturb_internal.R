#' Determine if a disturbance occurred and reduce population abundances in a time step.
#'
#' @param N Abundance matrix with `nrow = n_sp = 3` and `ncol = n_patch`.
#' @param disturb_value Vector of numbers from 0 to 1 describing the magnitude of disturbance `i`.
#' This is optimized to work with branching river networks generated from the the `mcbrnet` package by taking the output of `df_patch$disturbance` from `brnet()` as an input. See `?brnet` for more details.
#' For 2D habitats, this number is sampled randomly by converting the `disturb_mag_mean` and `disturb_mag_sd` arguments using the inverse logit function. Currently, the value for each patch is randomly sampled, i.e., there is no spatial correlation. This has not been optimized for 2D habitats
#' @param disturb_type Character vector of one of `c("regional", or NULL)`. `"regional"` applies a disturbance to all patches within the network.
#' If `disturb_type = NULL`, no disturbances will occur in the simulation.
#' Future releases are intended to include an additional disturbance type for branching river networks: `"point-source"`. Which is intended to apply a disturbance randomly to individual patches, and the magnitude of disturbance will decay downstream. This is not currently implemented.
#' @param disturb_p probability of disturbance occurring. If `disturb_type = "regional"`, disturbance for the entire meta-community is determined in each time step via `rbinom(n = 1, size = 1, prob = disturb_p)`.
#' @param river_network_structure Logical indicating if habitat architecture is branching (`TRUE`) or a 2D square (`FALSE`)


#' @return List with two elements. `N` = abundance matrix after accounting for disturbances. Values don't necessarily have to be integers. Numeric-double values will be converted to integer with `rpois()` in the simulation function. `patch_extinction` is an integer vector of `length = n_patch` indicating if a disturbance did happen (1) or did not happen (0).
#'
#' @export
#' @importFrom stats rbinom
#'
#' @examples
#' disturb_internal(N, disturb_type, disturb_p, disturb_value, river_network_structure)
#'
disturb_internal <- function(N,
                             disturb_type,
                             disturb_p,
                             disturb_value,
                             river_network_structure){
  n_patch <-  ncol(N)
  # no disturbances
  if(is.null(disturb_type)){
    N <- N
    patch_extinction <-  rep(0, n_patch)
  } else {
    if(disturb_type == "regional"){
      if(river_network_structure == FALSE){
        patch_extinction <- rbinom(n = 1, size = 1, prob = disturb_p)
      if(patch_extinction == 1){
        # Number of individuals after disturbance
        N <- N*(1 - disturb_value)[col(N)]
        patch_extinction <- rep(patch_extinction, n_patch)
      }
        }
      if(river_network_structure == TRUE){
        if(is.null(disturb_value)){
          stop(
            "disturb_type = regional but no `disturb_value` supplied")
        }
        patch_extinction <- rbinom(n = 1, size = 1, prob = disturb_p)
        if(patch_extinction == 1){
          N <- N*(1 - disturb_value)[col(N)]
          patch_extinction <- rep(patch_extinction, n_patch)
        }
      }
    }
  }
  N[N<0] <- 0
  return(list(N = N, patch_extinction = patch_extinction))
}


#' Calculate K based on environmental value
#'
#' @description function used internally by `k_function_internal` when `igp_sim(k_function = "environment")`. Inherits all arguments from `igp_sim()`.
#'
#' @param env Vector of environmental values. Inherited from `igp_sim()` function.
#' @param k_base Mean carrying capacity K.
#'
#' @details This function adds variation to the carrying capacity based on the environment value. When the variation in the `env` variable is high, the carrying capacity varies between 0.5 and 1.5 * `k_base`, with most of the density being concentrated at the extremes. When variation in `env` is very low, the observed k is approximately = `k_base`.
#'
#'  This function is designed to take the `df_patch$environment` output from the `mcbrnet::brnet()` as the `env` argument. See the `brnet()` function to control the variation of `environment` in head waters as well as longitudinal variation of environmental values between patches.
#'
#' @return Vector of carrying capacities
#'
#' @export
#'
#' @examples
#' env_to_k(env = rnorm(10), k_base = 100)
#'
env_to_k <- function(env,
                     k_base){
  if(any(is.na(c(env, k_base)))){
    stop("all inputs to `env_to_k()` need to be defined, \n one or more input values is 'NA' or 'NaN'" )
  }
  if(length(env) ==1){
    k = k_base
    k
  } else{
    # transform env to inverse logit scale
    # add 0.5 to put values between [0.5, 1.5]
  env_scale = (exp(env)/(1+exp(env))) + 0.5
  k = round(env_scale * k_base)
  # return vector of observed k values per patch
  k}
}

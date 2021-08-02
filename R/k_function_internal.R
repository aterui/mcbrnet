#' Internal function used to calculate carrying capacity within `igp_sim()`
#'
#' @description Function used internally by `igp_sim()`. Inherits all arguments from `igp_sim()` function call, or others calculated internally.
#'
#' @param k_function Character string, one of c("patches upstream", "environment", or NULL).
#' @param k_base Numeric input scaling the carrying capacity. See details for more information.
#' @param r_max Maximum reproductive rate of basal species, B
#' @param n_upstream Vector of length = n_patch describing the number (n+1) of upstream patches. Terminal upstream patches have a value of 1.
#' @param n_patch The number of total patches
#' @param k_c constant (default = 10) used in `k_n_upstream()` when `k_function =  "patches-upstream"`
#' @param k_min_exponent min exponent (default = 1.1) used in `k_n_upstream()` when `k_function =  "patches upstream"`
#' @param k_max_exponent min exponent (default = 1.35) used in `k_n_upstream()` when `k_function =  "patches upstream"`
#' @param river_network_structure Logical indicating if this is a branching river network `TRUE` or not `FALSE`.
#' @param environment_value Vector of environmental values. Designed to take `brnet()` function output.
#'
#' @return `k` vector of length = `n_patch` which describes carrying capacity for each patch
#' `b` parameter describing density dependent growth of basal species B.
#' @export
#'
#' @examples
#' k_function_internal(k_function, k_base, r_max, n_upstream, n_patch, k_c, k_min_exponent, k_max_exponent, river_network_structure, environment_value)
#'
k_function_internal <- function(k_function,
                                k_base,
                                r_max,
                                n_upstream,
                                n_patch,
                                k_c,
                                k_min_exponent,
                                k_max_exponent,
                                river_network_structure,
                                environment_value){
  if(is.null(k_function)){
    if(is.null(k_base)){
      stop("k_function is `NUll` and k_base argument is empty; \nsupply value for k_base")
    }
    if(length(k_base) == 1){
      message("only one value of k_base supplied, assuming it is the same in all patches")
      return(list(k = k_base,
                  b = calc_b(r_max = r_max, k = k_base)))
    }} else{  # for branching river networks
      if(k_function == "patches-upstream"){
        if(river_network_structure == FALSE){
          message("***CAUTION*** \nk_function = `patches upstream` but distance matrix for river network not supplied;\nrelationship between k ~ patches uncertain")
          n_upstream = rexp(n_patch) *n_patch / 6
        }
        if(is.null(n_upstream))
          stop("number of patches upstream must be supplied to n_upstream argument")
        if(length(n_upstream)!= n_patch)
          stop("length of n_upstream needs to equal n_patch")
        b_k_list = k_n_upstream(k_base = k_base,
                                k_c = k_c,
                                k_min_exponent = k_min_exponent,
                                k_max_exponent = k_max_exponent,
                                r_max = r_max,
                                n_upstream = n_upstream,
                                n_patch = n_patch)
        return(b_k_list)
      }
      if(k_function == "environment"){
        if(is.null(environment_value))
          stop("need to supply vector of environmental values")
        if(length(environment_value)!= n_patch)
          stop("length of `environment_value` should equal n_patch")
        k = env_to_k(env = environment_value, k_base = k_base)
        b = calc_b(r_max = r_max, k = k)
        message("environmental value scaled from 50% to 150% of k_base")
        return(list( k=k, b = b))
      }
    }
}

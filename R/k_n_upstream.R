
#' Carrying Capacity for basal species B, based on number of upstream nodes
#'
#' @description internal function used by `k_function_internal` when `igp_sim(k_function = "patches-upstream")`. Inherits all arguments from `igp_sim()`.
#'
#' @param k_base integer
#' @param k_c numeric
#' @param k_min_exponent minimum value for exponent describing productivity ~ watershed area relationship. default is 1.10 (Koening et al. 2019) numeric
#' @param k_max_exponent minimum value for exponent describing productivity ~ watershed area relationship. default is 1.35 (Koening et al. 2019) numeric
#' @param r_max maximum reproductive rate, numeric
#' @param n_upstream number of patches upstream + 1. i.e., the most upstream terminal patch has `n_upstream = 1`, integer
#' @param n_patch number of patches to calculate, should be same length as n_upstreamm, integer
#'
#'@importFrom stats runif
#'
#'@details The minimum carrying capacity is equal to `k_base` + `k_c`. The upstream terminal node The exponent for each patch is sampled from a uniform distribution from `k_min_exponent` to `k_max_exponent`.
#'k = k_base + k_c * n_upstream^k_exp
#'K is automatically rounded to the nearest integer.
#'b = (r_max - 1) /k
#'
#' @return list with 2 vector elements, each the length of `n_patch`. `k` = carrying capacity for that patch. `b` = parameter controlling realized strength of population growth at a given population size. b is used internally in `igp_sim()`
#' @export
#'
#' @examples k_n_upstream(k_base = 150, k_c = 10, k_min_exponent = 1.10,
#'  k_max_exponent = 1.35, r_max = 4,
#'   n_upstream = c(5, 10), n_patch = 2)
#'
#' @references Koenig, L.E., Helton, A.M., Savoy, P., Bertuzzo, E., Heffernan, J.B., Hall, R.O., Jr. and Bernhardt, E.S. (2019), Emergent productivity regimes of river networks. Limnol Oceanogr, 4: 173-181. https://doi.org/10.1002/lol2.10115
#'
k_n_upstream <- function(k_base = 150,
                         k_c = 10,
                         k_min_exponent = 1.10,
                         k_max_exponent = 1.35,
                         r_max,
                         n_upstream,
                         n_patch){

  # carrying capacity ####
  # power law with exponent between 1.1-1.35 (Koening et al. 2019)
  # k_base = "minimum" carrying capacity
  # k_c = constant multiplier for power law

  # k_exp = exponent for carrying capacity function

  k_exp = runif(n = n_patch,
                min = k_min_exponent,
                max = k_max_exponent)
  k <- round(
    k_base + k_c * n_upstream^k_exp)
  b = calc_b(r_max = r_max, k = k)
  list(k = k,
       b = b)
}

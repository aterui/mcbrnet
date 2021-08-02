#' Calculate density dependent growth parameter
#'
#' @description Internal function used by `k_function_internal` to calculate parameter `b`, which determines asymptotic level of basal species, B, recruitment
#'
#' @details b = (r_max - 1) / k
#'
#' @param r_max maximum per capita recruitment rate for basal species B
#' @param k carrying capacity
#'
#' @return parameter b
#' @export
#'
#' @examples
#' calc_b(r_max = 2, k = 100)
#'
calc_b <- function(r_max, k){
  (r_max - 1) / k
}

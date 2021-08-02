#' Maximum reproductive rate and carrying capacity for all three species
#'
#' @description This is a helper function to estimate reproductive rates and carrying capacity of species C and P.
#'
#' @param k_b Carrying capacity of basal species B
#' @param r_b Maximum reproductive rate of basal species, B
#' @param ebc conversion efficiency of turning B biomass into new C biomass
#' @param ebp conversion efficiency of turning B biomass into new P biomass
#' @param ecp conversion efficiency of turning C biomass into new p biomass
#' @param alphabc parameter controlling predation and reproductive rates of consumer species, C
#' @param alphap parameter controlling predation and reproductive rates of consumer species, P
#' @param betabc parameter controlling predation and reproductive rates of consumer species, C
#' @param betap parameter controlling predation and reproductive rates of consumer species, P
#' @param P_pref The preference that species P has for B over C. Default is `NULL`, and will be calculated based on other parameters (see below). This value can be fixed by the user i.e. `P_pref = 0.25` indicates that the predator spends ~25% of it's time searching for B, and `1 - P_pref = 75%`  of its time searching for C.
#'
#' @details This function calculates the maximum reproductive rate and carrying capacity for each species in the IGP simulation communities. This can be used to estimate what effect different parameters will have. Carrying capacities are calculated assuming prey resources are at their respective carrying capacities, hence the return values estimate upper ceiling of population abundances for a given set of parameters. This function does also not account for indirect effects, i.e. C predation upon B lowering resource abundance for P.
#'
#' Carrying capacity and maximum reproductive rate of the basal species is set by the user.
#'
#' Maximum reproductive rates for C are calculated as `r_c = ebc * alphac`. Carrying capacity is calculated as `k_c = k_b*(r_c - 1) / betac`
#'
#' Maximum reproductive rates for the predator are calculated separately for both of its prey, and depends on its preference of prey resource B over C. If `P_pref = NULL`, the default, the preference is calculated according to the `prey_preference()` function.
#'
#'
#'
#' @return a list with 3 named elements: `max_r` which has the maximum reproductive rates for B (`r_b`) and C (`r_c`), and the partitioned (`r_bp`, and `r_cp`) and total rates for P (`r_p_total`); `k` which has the carrying capacity for B (`k_b`) and C (`k_c`), and the partitioned and total rates for P; `input` which has the input values used for calculations.
#' @export
#'
#' @seealso `prey_preference`
#'
#' @examples
#' r_k_sp(k_b = 500, r_b = 2.5, ebc = 2,ebp = 2, ecp = 2, alphac = 4, alphap = 4, betac = 20, betap = 20, P_pref = NULL)
#'
r_k_sp <- function(k_b = 150,
                   r_b = 2.5,
                   ebc = 2,
                   ebp = 2,
                   ecp = 2,
                   alphabc = 4,
                   alphap = 4,
                   betabc = 20,
                   betap = 20,
                   P_pref = NULL){

  r_c = ebc * alphabc
  k_c = k_b*(r_c - 1) / betabc

  if(is.null(P_pref)){
    P_pref = prey_preference(e1 = ebp, e2 = ecp, N1 = k_b, N2 = k_c)
  }

  alphabp = P_pref * alphap
  alphacp = (1 - P_pref) * alphap

  r_bp = ebp*alphabp
  r_cp = ecp*alphacp
  r_p_total = r_bp + r_cp

  k_bp = k_b*(r_bp - 1) / betap
  if (k_bp < 0){
    k_bp = 0
  }

  k_cp = k_c*(r_cp - 1) / betap
  if (k_cp < 0){
    k_cp = 0
  }

  k_p_total = k_cp + k_bp

  max_r = signif(
    c(r_b = r_b, r_c = r_c,
      r_bp = r_bp, r_cp = r_cp, r_p_total = r_p_total),
    digits = 2)
  k = round(
    c(k_b = k_b, k_c = k_c,
        k_bp = k_bp, k_cp = k_cp,
        k_p_total = k_p_total))
  input = signif(
    c(k_b = k_b,
      r_b = r_b,
      ebc = ebc,
      ebp = ebp,
      ecp = ecp,
      alphabc =alphabc,
      alphap =alphap,
      betabc = betabc,
      betap = betap,
      P_pref = P_pref),
    digits = 2)
  list(max_r = max_r, k = k, input = input)
}

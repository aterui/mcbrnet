#' population dynamic function for each time step
#'
#' @description Internal function used by `igp_sim()` to simulate population dynamics for each time step. Inherits all arguments from `igp_sim()` call or other variables calculated internally.
#'
#' @param N Abundance matrix at time t
#' @param P_pref value for the predator preference of B over C, numeric between 0-1 or NULL
#' @param fixed_P_pref Logical indicating if P_pref is a fixed value `TRUE`, or if it varies based on prey abundances `FALSE`. This is determined internally based on the input to `P_pref` in `igp_sim()`
#' @param alphabc Parameter controlling the predation of resource B by consumer C
#' @param betabc Parameter controlling the predation of resource B by consumer C
#' @param ebc conversion efficiency of turning consumed B into new C
#' @param alphap Parameter controlling the predation of both resources by consumer P
#' @param betap Parameter controlling the predation of both resources by consumer P
#' @param ebp conversion efficiency of turning consumed B into new P
#' @param ecp conversion efficiency of turning consumed C into new P
#' @param v_s0 vector of base survival probabilities. created internally based on s0 argument in `igp_sim()`
#' @param r_max maximum per capita reproduction rate of basal species B
#' @param b parameter controlling asymptotic level of recruitment, calculated as (r_max - 1) / k
#' @param k carrying capacity for a patch.
#'
#' @details This function estimates the continuous effects of predation and survival probability through one time step (i.e. 1 season), and discrete reproduction at the end of the time step.
#'
#' @return `N` Abundance matrix N, at time t+1.
#' `lambda_b` Proportion of B in P's diet. Used to caluclate food chain length.
#' @export
#'
#' @examples
pop_sim <- function(N,
                    P_pref,
                    fixed_P_pref,
                    alphabc,
                    betabc,
                    ebc,
                    alphap,
                    betap,
                    ebp,
                    ecp,
                    v_s0,
                    b,
                    k,
                    r_max){
  # pop_sim() ####
  # predation ####
  if(fixed_P_pref == TRUE){
    if(is.null(P_pref)){
      stop("fixed_P_pref == TRUE, but no value for `P_pref` supplied")
    }
    P_pref = P_pref
  }
  if(fixed_P_pref == FALSE){
    P_pref = prey_preference(e1 = ebp, e2 = ecp, N1 = N[1,], N2 = N[2,])
  }

  # number of prey consumed
  #wij = number of prey i eaten by predator j
  wbc = alphabc * N[1,] * N[2,] / (betabc * N[2,] + N[1,])
  wbp = P_pref * (alphap * N[1,] * N[3,]/(betap * N[3,] + N[1,]))
  wcp = (1 - P_pref) * (alphap * N[2,] * N[3,]/(betap * N[3,] + N[2,]))

  # lambda: proportion of B in P diet
  w_dot_p = wbp + wcp
  lambda_b = wbp / w_dot_p

  # when P is exctinct, lamba_b = NaN
  # remove these from output
  lambda_b[is.nan(lambda_b)] <- 0

  # survival ####
  B_prime = v_s0[1] * (N[1,] - wbc - wbp)
  C_prime = v_s0[2] * (N[2,] - wcp)
  P_prime = v_s0[3] * N[3,]

  # reproduction ####
  B_t1 = (r_max / (1 + b * B_prime)) * B_prime
  C_t1 = (ebc * alphabc * N[1,] /
            (betabc * N[2,] + N[1,])) * # B converted to C
    C_prime # number of C
  P_t1 = ((P_pref * (ebp * alphap * N[1,] /
                       (betap * N[3,] + N[1,]))) + # B converted to P
            ((1 - P_pref) *
               (ecp * alphap * N[2,] /
                  (betap * N[3,] + N[2,])))) * # C converted to P
    P_prime # number of P

  N = rbind(B_t1, C_t1, P_t1)

  N[is.nan(N)] <- 0
  N[is.na(N)] <- 0
  N[N <0 ] <- 0
  return(list(N = N, lambda_b = lambda_b))
}

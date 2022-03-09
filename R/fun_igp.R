#' Internal function: igp function
#'
#' @param x community matrix for basal, intraguild prey (ig-prey), and intraguild predator (ig-predator)
#' @param r_b maximum reproductive rate of basal species
#' @param e energetic conversion efficiency; must be given by the order of basal to ig-prey, basal to ig-predator, and ig-prey to ig-predator
#' @param k carrying capacity
#' @param a attack rate; must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param h handling time; ; must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param s0 background survival rate; must be given by the order of basal, ig-prey, and ig-predator
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_igp <- function(x,
                    r_b = 5,
                    k = 100,
                    e = c(5, 5, 5),
                    a = c(0.5, 0.5, 0.5),
                    h = c(1, 1, 1),
                    s0 = rep(0.8, 3)
) {

  # check inputs ------------------------------------------------------------

  if (!is.matrix(x)) stop("error in x; x must be matrix")
  if (dim(x)[1] != 3) stop("error in x; dim(x)[1] must be 3")
  if (any(a > 1)) warning("a > 1; may yield negative density")
  if (a[1] + a[2] > 1) warning("a[1] + a[2] > 1; may yield negative density")
  if (any(a < 0) |
      any(e < 0) |
      any(x < 0) |
      any(h < 0) |
      any(s0 < 0) |
      any(s0 > 1)
      ) stop("negative values detected in parameters or community matrix")

  # define parameters -------------------------------------------------------

  ## basal species, intraguild prey, intraguild predator
  v_n_b <- x[1, ]
  v_n_c <- x[2, ]
  v_n_p <- x[3, ]

  ## growth rate / conversion efficiency
  e_bc <- e[1]
  e_bp <- e[2]
  e_cp <- e[3]

  ## attack rate
  a_bc <- a[1]
  a_bp <- a[2]
  a_cp <- a[3]

  ## handling time
  h_bc <- h[1]
  h_bp <- h[2]
  h_cp <- h[3]

  ## background survival
  s0_b <- s0[1]
  s0_c <- s0[2]
  s0_p <- s0[3]


  # trophic relationship ----------------------------------------------------

  # functional response ####
  ## v_w_xx: total number of prey eaten by predators
  ### denominator
  den_bc <- (v_n_c + a_bc * h_bc * v_n_b)
  den_bp <- (v_n_p + a_bp * h_bp * v_n_b)
  den_cp <- (v_n_p + a_cp * h_cp * v_n_c)

  v_w_bc <- ifelse(den_bc > 0, (a_bc * v_n_b * v_n_c) / den_bc, 0)
  v_w_bp <- ifelse(den_bp > 0, (a_bp * v_n_b * v_n_p) / den_bp, 0)
  v_w_cp <- ifelse(den_cp > 0, (a_cp * v_n_c * v_n_p) / den_cp, 0)

  ## v_f_xx: number of prey eaten per predator
  v_f_bp <- (a_bp * v_n_b) / (1 + a_bp * h_bp * v_n_b)
  v_f_cp <- (a_cp * v_n_c) / (1 + a_cp * h_cp * v_n_c)

  ## v_delta: preference to basal over intraguild prey
  zeta <- e_bp * v_f_bp + e_cp * v_f_cp
  den_delta <- ifelse(zeta > 0, zeta, 1)
  v_delta <- (e_bp * v_f_bp) / den_delta

  ## v_z_xx: total proportional loss = preference X proportion eaten
  ## xi_xx: indicator
  v_w_bcp <- v_w_bc + v_delta * v_w_bp
  xi_bcp <- ifelse(v_w_bcp > v_n_b, 1, 0)

  v_z_bc <- xi_bcp * (v_w_bc / v_w_bcp) + (1 - xi_bcp) * (v_w_bc / v_n_b)
  v_z_bp <- xi_bcp * ((v_delta * v_w_bp) / v_w_bcp) + (1 - xi_bcp) * ((v_delta * v_w_bp) / v_n_b)
  v_z_cp <- (1 - v_delta) * (v_w_cp / v_n_c)

  v_z_bc[is.nan(v_z_bc)] <- 0
  v_z_bp[is.nan(v_z_bp)] <- 0
  v_z_cp[is.nan(v_z_cp)] <- 0


  # trophic dynamics --------------------------------------------------------

  ## v_n_x_bar: density after background + predation-induced mortality
  ## v_n_x_hat: density after reproduction

  # basal species v_n_b ####
  ## survived = basal survival x predation-induced mortality
  nu <- (r_b - 1) / k
  v_n_b_bar <- s0_b * (1 - v_z_bc - v_z_bp) * v_n_b
  v_n_b_hat <- (v_n_b_bar * r_b) / (1 + nu * v_n_b_bar)

  # intraguild prey v_n_c ####
  ## survived = basal survival x predation-induced mortality
  v_n_c_bar <- s0_c * (1 - v_z_cp) * v_n_c
  v_n_c_hat <- e_bc * v_w_bc

  # intraguild predator v_n_p ####
  v_n_p_bar <- s0_p * v_n_p
  v_n_p_hat <- v_delta * e_bp * (v_z_bp * v_n_b) + (1 - v_delta) * e_cp * v_w_cp
  v_n_p_hat[is.nan(v_n_p_hat)] <- 0


  # export ------------------------------------------------------------------

  m_z <- rbind(v_z_bc,
               v_z_bp,
               v_z_cp)

  dimnames(m_z) <- NULL

  m_n_bar <- rbind(v_n_b_bar,
                   v_n_c_bar,
                   v_n_p_bar)

  dimnames(m_n_bar) <- NULL

  m_n_hat <- rbind(v_n_b_hat,
                   v_n_c_hat,
                   v_n_p_hat)

  dimnames(m_n_hat) <- NULL

  return(list(delta = v_delta,
              z = m_z,
              m_n_bar = m_n_bar,
              m_n_hat = m_n_hat)
  )

}

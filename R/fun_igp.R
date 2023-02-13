#' Internal function: igp function
#'
#' @param x Community matrix for basal, intraguild prey (ig-prey), and intraguild predator (ig-predator)
#' @param r_b Maximum reproductive rate of basal species
#' @param e Energetic conversion efficiency; must be given by the order of basal to ig-prey, basal to ig-predator, and ig-prey to ig-predator
#' @param k Carrying capacity of basal species
#' @param a Attack rate. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param h Handling time. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param s Strength of prey switching from ig-prey to basal
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_igp <- function(x,
                    r_b = 5,
                    k = 100,
                    e = c(0.8, 0.8, 0.8),
                    a = c(2, 2, 2),
                    h = c(0.1, 0.1, 0.1),
                    s = 0
) {

  # check inputs ------------------------------------------------------------

  if (!is.matrix(x)) stop("error in x; x must be matrix")
  if (dim(x)[1] != 3) stop("error in x; dim(x)[1] must be 3")
  if (length(s) != 1) stop("error in s; s must be a scalar")
  if (s < 0 | s > 1) stop("error in s; s must be 0 - 1")

  if (any(x < 0) |
      any(r_b < 0) |
      any(k < 0) |
      any(a < 0) |
      any(e < 0) |
      any(h < 0)
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

  ## attack rate x handling time
  h_bc <- a[1] * h[1]
  h_bp <- a[2] * h[2]
  h_cp <- a[3] * h[3]


  # trophic relationship ----------------------------------------------------

  ## v_p_bc: fraction of prey survived after predation by consumer
  v_p_bc <- exp(-((a_bc * v_n_c) / (1 + h_bc * v_n_b)))

  ## predation C on B & conversion
  v_n_b_minus_c <- v_n_b * v_p_bc
  v_n_c_plus <- e_bc * v_n_b * (1 - v_p_bc)

  ## preference function P on B & C
  ## phi: switching function
  ## s: strength of switching
  phi <- s * (v_n_b_minus_c - v_n_c_plus) / (v_n_b_minus_c + v_n_c_plus)
  phi <- ifelse(is.nan(phi), 0, phi) # NaN produced when n_b = n_c = 0
  om_b <- 1 + phi # preference to basal over ig-prey
  om_c <- 1 - phi # preference to ig-prey over basal

  ## v_p_xp: fraction of prey survived after predation by predator
  v_p_bp <- exp(-((om_b * a_bp * v_n_p) / (1 + h_bp * v_n_b_minus_c)))
  v_p_cp <- exp(-((om_c * a_cp * v_n_p) / (1 + h_bp * v_n_c_plus)))


  # trophic dynamics --------------------------------------------------------

  ## Basal species
  ### predation P on B
  v_n_b_minus_cp <- v_n_b_minus_c * v_p_bp

  ### beverton-holt growth
  nu <- (r_b - 1) / k
  v_n_b_hat <- v_n_b_minus_cp * (r_b / (1 + nu * v_n_b_minus_cp))

  ## IG prey
  ### conversion x basal n x faction captured x predation P on C
  v_n_c_hat <- v_n_c_plus * v_p_cp

  ## IG predator
  ### conversion x (basal n x fraction remained) x fraction captured
  ### conversion x consumer n x fraction captured
  v_n_p_hat <- e_bp * v_n_b_minus_c * (1 - v_p_bp) + e_cp * v_n_c_plus * (1 - v_p_cp)

  # omnivory level
  v_delta <- (e_bp * v_n_b_minus_c * (1 - v_p_bp)) / v_n_p_hat
  v_delta[is.nan(v_delta)] <- 0

  # export ------------------------------------------------------------------

  m_n_hat <- rbind(v_n_b_hat,
                   v_n_c_hat,
                   v_n_p_hat)

  dimnames(m_n_hat) <- NULL

  return(list(m_n_hat = m_n_hat,
              delta = v_delta))

}

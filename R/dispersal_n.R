#' Calculate number of individuals dispersing from patch x to y
#'
#' @description Internal function used by `igp_sim()` to calculate the number of individuals dispersing from patch x, and distributing them to patches y based on distanced using an exponential decay kernel. This function is called at the beginning of each timestep in the dynamic simulation. Inherits all arguments from `igp_sim()` or other variables calculated internally.
#'
#' @param N Abundance matrix with `nrow = n_sp = 3` and `ncol = n_patch`. Rows 1:3 correspond to species B, C, and P respectively.
#' @param v_p_dispersal A vector of dispersal probabilities of length = 3. in `igp_sim()`, this is controlled with the `p_dispersal` argument. Values should be from 0 to 1
#' @param v_theta Vector of length = 3, Parameters controlling the distance decay function. Larger values of theta result in a faster decline of dispersal distance. `v_theta` is caluclated internally in `igp_sim` based on the `theta` argument.
#' @param dist_mat a distance matrix with `nrow = ncol = n_patch` and diagonal == 0.
#' @param m_b_dispersal dispersal matrix for species B. This and the following 2 matrices are calculated internally in `igp_sim()`
#' @param m_c_dispersal dispersal matrix for species C.
#' @param m_p_dispersaldispersal matrix for species P.
#'
#' @details This function is used internally in `igp_sim()` to calculate the number of individuals of each dispersing from patch i to patch j. Emigrants from patch i are more likely to become Immigrants to patch j if the patches are closer together.
#'
#' @return `m_n_prime` A matrix of abundances after accounting for dispersal with `nrow = n_sp = 3` and `ncol = n_patch`. Values may be returned with decimal values, but will be converted to integers using `rpois()` within the simulation model.
#'
#' @export
#'
#' @examples
#' n_patch = 5
#' n_sp = 3
#' N = matrix(rpois(n_sp * n_patch, lambda = 10), nrow = n_sp, ncol = n_patch)
#' # sort coordinates so you can easily see increasing distances
#'  x_coord = sort(runif(5, 0, 5))
#'  y_coord = sort(runif(5, 0, 5))
#'  dist_mat = data.matrix(dist(cbind(x_coord, y_coord)))
#'  v_theta = c(1, 1, 1)
#'  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
#'  diag(m_b_dispersal) <- 0
#' # species C
#' m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
#' diag(m_c_dispersal) <- 0
#' # species P
#' m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
#' diag(m_p_dispersal) <- 0
#'
#' disperal_n(N = N,
#'  v_p_dispersal = 0.25,
#'  v_theta = 1,
#'  dist_mat = dist_mat,
#'  m_b_dispersal = m_b_dispersal,
#'  m_c_dispersal = m_c_dispersal,
#'  m_p_dispersal = m_p_dispersal)
#'
#'
dispersal_n <- function(N,
                        v_p_dispersal,
                        v_theta,
                        dist_mat,
                        m_b_dispersal,
                        m_c_dispersal,
                        m_p_dispersal){

  test1 = all(v_p_dispersal <=1) & all(v_p_dispersal>=0)
  if (test1 == FALSE){
    stop("in dispersal_n, v_p_dispersal must be between 0 and 1")
  }


  # estimated Emmigrants
  m_e_hat <- N * v_p_dispersal
  v_e_sum <- rowSums(m_e_hat)
  # distribute Emigrants to new patches as Immigrants
  # Separate calculation for each species
  m_i_b_raw <- m_e_hat[1,] %*% m_b_dispersal
  m_i_c_raw <- m_e_hat[2,] %*% m_c_dispersal
  m_i_p_raw <- m_e_hat[3,] %*% m_p_dispersal
  # combine immigrant totals
  m_i_raw <- rbind(m_i_b_raw, m_i_c_raw, m_i_p_raw)


  v_i_sum <- rowSums(m_i_raw)
  v_i_sum[v_i_sum == 0] <- 1
  m_i_prob <- m_i_raw / v_i_sum
  m_i_hat <- m_i_prob * v_e_sum
  # abundance matrix, N, after dispersal occurs
  m_n_prime <- N + m_i_hat - m_e_hat

  # make sure all values are numeric and >=0
  m_n_prime[is.nan(m_n_prime)] <- 0
  m_n_prime[is.na(m_n_prime)] <- 0
  m_n_prime[m_n_prime <0 ] <- 0

  m_n_prime
}

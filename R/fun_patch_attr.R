#' Internal function: weighted-mean patch attributes
#'
#' @param x adjacency matrix to be converted
#' @param n_branch number of branch
#' @param mean_source mean value of source attribute
#' @param sd_source SD of source attribute
#' @param sd_lon SD of longitudinal change
#' @param m_distance distance matrix
#' @param rho longitudinal autocorrelation
#' @param v_wa vector of watershed area
#' @param logit logit transformation of mean source attribute
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_patch_attr <- function(x,
                           n_branch,
                           mean_source,
                           sd_source,
                           sd_lon,
                           m_distance,
                           rho = 1,
                           v_wa,
                           logit = FALSE) {

  # extract upstream adjacency
  m_adj_up <- x
  m_adj_up[lower.tri(m_adj_up)] <- 0
  n_patch <- nrow(x)

  # upstream tributary proportion
  m_wa_prop <- t(apply(X = m_adj_up,
                       MARGIN = 1,
                       function(x) {
                         (x * v_wa) / ifelse(sum(x) == 0,
                                             yes = 1,
                                             no = sum(x * v_wa))
                         }
                       )
                 )

  n_source <- 0.5 * (n_branch + 1)
  source <- which(rowSums(m_adj_up) == 0)

  v_z_dummy <- v_z <- v_attr <- rep(0, n_patch)
  v_z[source] <- v_attr[source] <- rnorm(n = n_source,
                                         mean = ifelse(logit == TRUE,
                                                       yes = boot::logit(mean_source),
                                                       no = mean_source),
                                         sd = sd_source)
  v_z_dummy[source] <- 1

  if (!(rho <= 1 & rho >= 0)) stop("rho must be between 0 and 1")

  for (i in seq_len(max(m_distance[1, ]))) {

    v_eps <- rep(0, n_patch)
    v_eps[v_z_dummy != 0] <- rnorm(n = length(v_eps[v_z_dummy != 0]),
                                   mean = 0,
                                   sd = sd_lon)

    v_z <- m_wa_prop %*% ((rho * v_z) + v_eps)
    v_z_dummy <- m_wa_prop %*% v_z_dummy
    v_attr <- v_z + v_attr

  }

  return(c(v_attr))
}

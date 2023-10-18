#' Preferential prey model
#'
#' @param n_species Integer. Number of species
#' @param n_basal Integer. Number of basal species
#' @param l Interger. Expected number of links in the upper triangle
#' @param theta Numeric. Scale parameter of an exponential distribution. Smaller values indicate greater trophic specialization.
#' @param cannibal Logical. If \code{TRUE}, cannibalism allowed
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

ppm <- function(n_species,
                n_basal,
                l,
                theta,
                cannibal = FALSE) {
  # n_species: number of species in the pool
  # n_basal: number of basal species
  # l: expected number of links
  # theta: scale parameter

  # verify inputs
  if (n_basal < 1) stop("At least one basal species needed to construct a food web")
  if (n_basal >= n_species) stop("n_basal must be smaller than n_species")

  # number of consumers
  n_c <- n_species - n_basal

  # tp: trophic position
  # assign trophic position for basal species
  # assign -1 for consumers as initial values (to be updated)
  tp <- rep(-1, n_species)
  tp[seq_len(n_basal)] <- 1

  # beta parameter for beta distribution
  # determined so that E(L) = L
  if (l < n_c)
    stop("l must be at least equal to the number of consumers (n_species - n_basal)")

  if (cannibal) {
    max_l <- sum((n_basal + 1):n_species)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 1)
      stop("n_species + n_basal must be equal to or greater than 1")

    b <- ((n_species + n_basal - 1) * n_c) / (2 * (l - n_c)) - 1
  } else {
    max_l <- sum((n_basal + 1):n_species - 1)

    if (l > max_l)
      stop(paste("maximum l is", max_l))

    if (n_species + n_basal < 3)
      stop("n_species + n_basal must be equal to or greater than 3")

    b <- ((n_species + n_basal - 3) * n_c) / (2 * (l - n_c)) - 1
  }

  # vector for link proportions for all consumers
  v_xi <- stats::rbeta(n_c,
                       shape1 = 1,
                       shape2 = b)

  # vector for initial prey choice for all consumers
  v_i0 <- c(rep(-1, n_basal),
            sapply(X = (n_basal + 1):n_species,
                   FUN = function(j) resample(seq_len(j - 1), size = 1)))

  # realized number of prey nodes for all consumers
  # kappa follows a beta-binomial distribution
  # size parameter is the possible maximum number of prey nodes
  if (cannibal) {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = n_basal:(n_species - 1),
                               prob = v_xi))
  } else {
    v_kappa <- c(rep(-1, n_basal),
                 stats::rbinom(n = n_c,
                               size = (n_basal - 1):(n_species - 2),
                               prob = v_xi))
  }

  # alpha0: S x S interaction binary matrix
  # initialized with all zero
  # update the initial consumer's first prey
  alpha0 <- matrix(0, n_species, n_species)
  alpha0[v_i0[n_basal + 1], n_basal + 1] <- 1

  for (j in (n_basal + 1):n_species) {
    alpha0 <- extra_prey(alpha0 = alpha0,
                         j = j,
                         i0 = v_i0[j],
                         tp = tp,
                         kappa = v_kappa[j],
                         theta = theta,
                         cannibal = cannibal)

    v_n_prey <- colSums(alpha0)
    v_n_prey[v_n_prey == 0] <- 1

    tp_new <- (tp %*% alpha0) / v_n_prey + 1
    tp[j] <- tp_new[j]
  }

  attr(alpha0, "tp") <- tp

  return(alpha0)
}

#' Extra prey function
#'
#' @inheritParams ppm
#' @param alpha0 Interaction matrix
#' @param j Integer. Consumer's index
#' @param i0 Integer. Index for the first prey
#' @param tp Numeric. Initial trophic positions
#' @param kappa Integer. Number of extra prey items
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

extra_prey <- function(alpha0,
                       j,
                       i0,
                       tp,
                       theta,
                       kappa,
                       cannibal = FALSE) {
  # alpha0:      adjacency matrix defining trophic interactions
  # j:      consumer's index (j - 1 is the number of possible prey)
  # i0:     index of the first prey chosen (1 =< i0 < j)
  # tp:     vector of trophic positions for possible prey nodes (length = S)
  # kappa:  number of extra prey nodes in addition to the first prey i0
  # theta:  scale parameter for an exponential decay of prey preference

  # verify inputs
  if (length(tp) != ncol(alpha0)) stop(paste("'tp' length seems incorrect;",
                                             "the length must be",
                                             ncol(alpha0)))

  if (!(i0 < j && 1 <= i0)) stop("i0 must be a non-zero intger smaller than j")

  if (kappa >= j) stop("kappa must be an integer smaller than j")

  # probability of consumer j picks prey i
  lambda <- 1 / theta
  tp[j] <- tp[i0] + 1
  p_ij <- exp(- lambda * abs(tp[i0] - tp))

  # pick kappa prey species out of j - 1 nodes (without cannibalism) or j nodes
  if (cannibal) {
    i_index <- 1:j
    p_ij <- p_ij[1:j]
  } else {
    i_index <- 1:(j - 1)
    p_ij <- p_ij[1:(j - 1)]
  }

  # exclude i0 from resample() because i0 is the first pick
  # pick 'kappa' samples from index[-i0]
  # weighted by p_ij
  if (length(i_index) > 1) {
    i_pick <- resample(i_index[-i0],
                       size = kappa,
                       prob = p_ij[-i0])

    i <- c(i0, i_pick)
  } else {
    # when j = 2 with no cannibalism
    # no choice but i0 available as prey
    i <- i0
  }

  alpha0[i, j] <- 1

  return(alpha0)
}

#' Apply conversion efficiency and attack rate
#'
#' @inheritParams extra_prey
#' @param attack List for attack rates. Specify minimum and maximum values for a uniform distribution.
#' @param convert List for conversion efficiency. Specify minimum and maximum values for a uniform distribution.
#' @param mortal List for mortality (or intraspecific competition). Specify minimum and maximum values for a uniform distribution.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

to_alpha <- function(alpha0,
                     attack = list(min = 0,
                                   max = 1),
                     convert = list(min = 0,
                                    max = 1),
                     mortal = list(min = 0,
                                   max = 1)) {

  # identify basal species
  u_alpha0 <- alpha0
  basal <- which(colSums(alpha0) == 0)

  # scale by # of resources
  su_alpha0 <- t(t(alpha0[, -basal]) / colSums(alpha0)[-basal])

  # generate random parameters
  ## matrix for attack rates
  a <- with(attack, matrix(stats::runif(ncol(alpha0)^2, min = min, max = max),
                           nrow = nrow(alpha0),
                           ncol = ncol(alpha0)))

  ## matrix for conversion efficiency
  b <- with(convert, matrix(stats::runif(ncol(alpha0)^2, min = min, max = max),
                            nrow = nrow(alpha0),
                            ncol = ncol(alpha0)))

  ## vector for intraspecific competition
  m <- with(mortal, stats::runif(ncol(alpha0), min = min, max = max))

  # interaction matrix
  u_alpha0[, -basal] <- a[, -basal] * su_alpha0
  alpha <- b * u_alpha0 - t(u_alpha0)

  # intraspecific interaction
  diag(alpha) <- -m

  return(alpha)
}

#' Generate a food web based on the preferential prey model
#'
#' @inheritParams ppm
#' @inheritParams to_alpha
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodweb <- function(n_species,
                    n_basal,
                    l,
                    theta,
                    cannibal = FALSE,
                    attack = list(min = 0,
                                  max = 1),
                    convert = list(min = 0,
                                   max = 1),
                    mortal = list(min = 0,
                                  max = 1)) {

  alpha0 <- ppm(n_species = n_species,
                n_basal = n_basal,
                theta = theta,
                l = l,
                cannibal = cannibal)

  alpha <- to_alpha(alpha0 = alpha0,
                    attack = attack,
                    convert = convert,
                    mortal = mortal)

  return(alpha)
}

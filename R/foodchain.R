#' Maximum trophic position
#'
#' @param n_species Integer. Number of species
#' @param n Numeric. Vector of species density
#' @param alpha Numeric. Interaction matrix
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

max_tp <- function(n_species, n, alpha) {

  alpha[lower.tri(alpha, diag = TRUE)] <- 0

  # initialize trophic positions for basal species
  tp <- rep(-1, n_species)

  id_basal <- which(colSums(alpha) == 0)
  n_basal <- length(id_basal)
  tp[id_basal] <- 1

  if (all(n[id_basal] == 0)) {

    # if basal species are all absent
    return(0)

  } else {

    # total gain from prey (prey converted to predator density)
    g <- drop(n %*% alpha)

    # fractional contribution of each prey
    ## frac excludes basal species
    frac <- do.call(rbind,
                    lapply(1:n_species,
                           function(i) {
                             alpha[i, (n_basal + 1):n_species] * n[i] / g[(n_basal + 1):n_species]
                           }))

    ## add columns for basal species (all zero)
    frac <- cbind(matrix(0, nrow = nrow(frac), ncol = n_basal), frac)

    # update trophic positions recursively
    for (i in (n_basal + 1):n_species) {
      value <- tp %*% frac + 1
      tp[i] <- value[i]
    }

    max_tp <- max(tp[n > 0])

    return(max_tp)

  }

}

#' Spatiotemporal average of maximum trophic positions using \code{\link{sglv}} output
#'
#' @param n Numeric. Matrix output from \code{\link{sglv}}
#' @param n_species Integer. Number of species
#' @param n_patch Integer. Number of habitat patches
#' @param alpha Numeric. Interaction matrix
#' @param start Integer. Initial time step to be included in the calculation
#' @param end Integer. Last time step to be included in the calculation
#' @param full Logical. If \code{TRUE}, a full spatiotemporal matrix of food chain length returned as attributes.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodchain <- function(n,
                      n_species,
                      n_patch,
                      alpha,
                      start = max(n[, 1]) - 99,
                      end = max(n[, 1]),
                      full = FALSE) {

  # verify input ------------------------------------------------------------
  if (!inherits(n, "deSolve"))
    stop("n must be class 'deSolve'")

  if (ncol(n) != n_species * n_patch + 1)
    stop(paste("Dimension mismatch; n_species x n_patch must be",
               ncol(n) - 1))

  if (!all(dim(alpha) == n_species))
    stop(paste("Dimension mismatch; alpha's dimentions must be",
               n_species,
               "x",
               n_species))

  # subset data -------------------------------------------------------------
  times <- seq(start, end, by = 1)
  index_t <- which(n[, "time"] %in% times)
  x <- n[index_t, -which(colnames(n) == "time")]

  # tidy format -------------------------------------------------------------
  ## time vector
  v_t <- rep(times, times = n_species * n_patch)

  ## species vector
  v_sp <- rep(seq_len(n_species), each = length(times) * n_patch)

  ## patch vector
  v_patch <- rep(rep(seq_len(n_patch), each = length(times)), times = n_species)

  ## combine: g specifies a unique combination of time and patch
  y <- data.table::data.table(t = v_t,
                              species = v_sp,
                              patch = v_patch,
                              g = as.numeric(paste0(v_t, v_patch)),
                              n = c(x))

  # food chain length -------------------------------------------------------
  ## vector - time (t) and patch-specific (j) food chains
  v_fcl <- with(y, tapply(y, INDEX = g,
                          FUN = function(d) {
                            with(d, max_tp(n_species = n_species,
                                           n = n,
                                           alpha = alpha)
                                 )
                          }))

  ## matrix - time (t) and patch-specific (j) food chains
  ## row: time, col: patch
  m_fcl <- matrix(v_fcl,
                  nrow = length(times),
                  ncol = n_patch,
                  byrow = TRUE)

  colnames(m_fcl) <- paste0("patch", seq_len(n_patch))
  rownames(m_fcl) <- times

  fcl <- mean(v_fcl)

  if (full) {
    attr(fcl, "fcl_matrix") <- m_fcl
  }

  return(fcl)
}

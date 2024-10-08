#' Maximum trophic position
#'
#' @param n_species Integer. Number of species
#' @param n Numeric vector. Vector of species density
#' @param alpha Numeric matrix. \code{n_species x n_species} interaction matrix
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

max_tp <- function(n_species, n, alpha) {

  # remove upper.tri()
  # lower.tri() represents energetic conversion from prey to predator
  alpha[upper.tri(alpha, diag = TRUE)] <- 0

  # initialize trophic positions for basal species
  tp <- rep(-1, n_species)

  id_basal <- which(rowSums(alpha) == 0)
  n_basal <- length(id_basal)
  tp[id_basal] <- 1

  if (all(n[id_basal] == 0)) {

    # if basal species are all absent
    max_tp <- 0

  } else {

    # total gain from prey (prey converted to predator density)
    g <- drop(alpha %*% n)

    # fractional contribution of each prey
    ## frac excludes basal species
    frac <- do.call(rbind,
                    lapply(1:n_species,
                           function(i) {
                             ## num: contribution of each prey
                             ## den: total prey converted into consumer
                             ## frac0: factional contribution of each prey
                             num <- alpha[(n_basal + 1):n_species, i] * n[i]
                             den <- g[(n_basal + 1):n_species]
                             frac0 <- num / den

                             ## 0 / 0 returns NaN
                             ## replace with 0 if num == 0 && den == 0
                             k <- intersect(which(num == 0), which(den == 0))
                             frac0[k] <- 0

                             return(frac0)
                           }))

    ## add columns for basal species (all zero)
    frac <- cbind(matrix(0, nrow = nrow(frac), ncol = n_basal), frac)

    # update trophic positions recursively
    for (i in (n_basal + 1):n_species) {
      value <- tp %*% frac + 1
      tp[i] <- ifelse(g[i] == 0, 0, value[i])
    }

    max_tp <- max(tp[n > 0])
  }

  return(max_tp)
}

#' Spatiotemporal average of maximum trophic positions using \code{\link{sglv}} output
#'
#' @param n Numeric. Matrix output from \code{\link{sglv}}
#' @param n_species Integer. Number of species
#' @param n_patch Integer. Number of habitat patches
#' @param alpha Numeric. Interaction matrix
#' @param start Integer. Initial time step to be included in the calculation. \code{floor} applied internally
#' @param end Integer. Last time step to be included in the calculation. \code{floor} applied internally
#' @param full Logical. If \code{TRUE}, a full spatiotemporal matrix of food chain length returned as attributes.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

foodchain <- function(n,
                      n_species,
                      n_patch,
                      alpha,
                      start = ceiling(max(n[, 1]) - max(n[, 1]) * 0.5) + 1,
                      end = floor(max(n[, 1])),
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
  times <- seq(floor(start), floor(end), by = 1)
  index_t <- which(n[, "time"] %in% times)
  x <- n[index_t, -which(colnames(n) == "time")]

  # tidy format -------------------------------------------------------------
  ## time vector
  v_t <- rep(times, times = n_species * n_patch)

  ## species vector
  v_sp <- rep(rep(seq_len(n_species), each = length(times)),
              times = n_patch)

  ## patch vector
  v_patch <- rep(seq_len(n_patch), each = length(times) * n_species)

  ## combine: g specifies a unique combination of time and patch
  y <- data.table::data.table(t = v_t,
                              species = v_sp,
                              patch = v_patch,
                              g = paste0("t",
                                         sprintf("%05d", v_t),
                                         "p",
                                         sprintf("%05d", v_patch)),
                              n = c(x))

  # food chain length -------------------------------------------------------
  ## vector - time (t) and patch-specific (j) food chains
  v_fcl <- tapply(y, INDEX = y$g,
                  FUN = function(d) {
                    max_tp(n_species = n_species,
                           n = d$n,
                           alpha = alpha)
                  })

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


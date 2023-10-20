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
    for(i in (n_basal + 1):n_species) {
      value <- tp %*% frac + 1
      tp[i] <- value[i]
    }

    max_tp <- max(tp[n > 0])

    return(max_tp)

  }

}


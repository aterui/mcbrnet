#' Internal distance matrix function
#'
#' @description Function used internally by `igp_sim()`. Checks that supplied distance matrix is in correct format and returns it. If no matrix is supplied, it calculates one and returns it.
#' @param dist_mat Distance matrix describing the spatial arrangement of patches. should be a square matrix with `ncol = nrow = n_patch` and `diagonal = 0`. If `dist_mat = NULL`, a square landscape with dimension `landscape_size` and `n_patch` randomly distributed will be calculated.
#' @param landscape_size Size of each side of randomly generated 2D habitat. Only used if `dist_mat = NULL`
#' @param n_patch the number of patches, inherited from the `igp_sim()` function argument.
#'
#' @return `dist_mat` distance matrix. Returns value supplied in `igp_sim()`, or if no value supplied, calculates one by randomly distributing `n_patch` in a square landscape of size = `landscape_size`
#' `river_network_structure` Logical. `TRUE` if habitat architecture is a branching river network, or `FALSE` if a square landscape.
#'
#' @export
#'
#' @examples
#' n_patch = 5
#' landscape_size = 10
#' x_coord = runif(n_patch, 0, landscape_size)
#' y_coord = runif(n_patch, 0, landscape_size)
#' dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
#'                             diag = TRUE, upper = TRUE))
#'
#' dist_mat_internal(dist_mat = dist_mat, landscape_size = 10, n_patch = n_patch)
#'
dist_mat_internal <- function(dist_mat,
                              landscape_size,
                              n_patch){
  if(!is.null(dist_mat)){
    if (!is.matrix(dist_mat))
      stop("distance matrix should be provided as matrix")
    if (nrow(dist_mat) != n_patch)
      stop(
        "invalid dimension: distance matrix must have a dimension of
    n_patch * n_patch")
    if (any(diag(dist_mat) != 0))
      stop(
        "invalid distance matrix: diagonal elements must be zero")
    # df_xy_coord <- NULL
    out = list(dist_mat = dist_mat,
         river_network_structure = TRUE)
    # add / make adj_matrix here
  }
  if(is.null(dist_mat)){
    if (is.null(landscape_size)){
      landscape_size = 10
    }
    message(
      paste("No distance matrix supplied, assuming a ",
            landscape_size, "x", landscape_size, " square landscape", sep = ""))
    x_coord = runif(n_patch, 0, landscape_size)
    y_coord = runif(n_patch, 0, landscape_size)
    dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
                                diag = TRUE, upper = TRUE))
    out = list(dist_mat = dist_mat,
         river_network_structure = FALSE)
  }
  out
}


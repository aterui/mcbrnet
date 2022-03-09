#' Internal function: to dispersal matrix
#'
#' @param n_patch number of patches
#' @param landscape_size landscape size
#' @param theta rate parameter for a dispersal kernel
#' @param xy_coord xy coordinate
#' @param distance_matrix distance matrix
#' @param dispersal_matrix dispersal matrix
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_disp_mat <- function(n_patch,
                         landscape_size,
                         theta,
                         xy_coord = NULL,
                         distance_matrix = NULL,
                         dispersal_matrix = NULL) {

  if (is.null(xy_coord) &
      is.null(distance_matrix) &
      is.null(dispersal_matrix)) {

    message("neither xy_coord nor distance_matrix is given: generate a square landscape with landscape_size (default: 10) on a side")

    v_x_coord <- runif(n = n_patch,
                       min = 0,
                       max = landscape_size)

    v_y_coord <- runif(n = n_patch,
                       min = 0,
                       max = landscape_size)

    df_xy_coord <- dplyr::tibble(x_coord = v_x_coord,
                                 y_coord = v_y_coord)

    m_distance <- data.matrix(dist(df_xy_coord,
                                   diag = TRUE,
                                   upper = TRUE))

    m_dispersal <- data.matrix(exp(-theta * m_distance))

  } else {

    if (!is.null(xy_coord) &
        is.null(distance_matrix) &
        is.null(dispersal_matrix)) {

      if (nrow(xy_coord) != n_patch) stop("row numbers must match n_patch")
      if (ncol(xy_coord) != 2) stop("the number of columns must be two, describing x- and y-cooridnates")

      colnames(xy_coord) <- c("x_coord", "y_coord")

      df_xy_coord <- dplyr::as_tibble(xy_coord)

      m_distance <- data.matrix(dist(df_xy_coord,
                                     diag = TRUE,
                                     upper = TRUE))

      m_dispersal <- data.matrix(exp(-theta * m_distance))

    } else {

      if (!is.null(xy_coord)) message("both xy_coord and distance matrix are given: argument xy_coord is ignored")
      if (!is.matrix(distance_matrix)) stop("distance matrix must be provided as matrix")
      if (nrow(distance_matrix) != n_patch) stop("invalid dimension: distance matrix must have a dimension of n_patch * n_patch")
      if (any(diag(distance_matrix) != 0)) stop("invalid distance matrix: diagonal elements must be zero")

      df_xy_coord <- NULL

      m_distance <- distance_matrix

      if (!is.null(dispersal_matrix)) {

        message("dispersal_matrix is provided: dispersal_matrix is used to simulate dispersal process")
        if (!is.matrix(dispersal_matrix)) stop("dispersal_matrix must be provided as matrix")
        if (nrow(dispersal_matrix) != n_patch) stop("invalid dimension: dispersal_matrix must have a dimension of n_patch * n_patch")
        if (any(diag(dispersal_matrix) != 0)) stop("invalid dispersal_matrix: diagonal elements must be zero")

        m_dispersal <- data.matrix(dispersal_matrix)

      } else {

        m_dispersal <- data.matrix(exp(-theta * m_distance))

      }

    }
  }

  diag(m_dispersal) <- 0

  return(list(m_distance = m_distance,
              m_dispersal = m_dispersal,
              df_xy_coord = df_xy_coord))

}

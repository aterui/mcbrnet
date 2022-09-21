#' Visualize a branching network
#'
#' @param x `brnet()` object
#' @param patch_color Type of patch (vertex) label (either \code{"env"}, \code{"disturbance"}, \code{"other"} or any color). Default \code{"env"}.
#' @param edge_color Type of edge label (\code{"passability"} or any color). Default \code{NULL}.
#' @param value_col Patch values. Must be specified if \code{patch_color = "other"}.
#' @param color_label Color legend title
#' @param patch_label Type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{"none"}, no label will be shown on patches in the plot. Default \code{"none"}.
#' @param patch_size Patch (vertex) size in the plot.
#' @param ... Arguments passed to \code{ggraph::geom_node_label()}.
#'
#' @importFrom rlang .data
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

ggbrnet <- function(x,
                    patch_color = "env",
                    edge_color = NULL,
                    value_col = NULL,
                    color_label = NULL,
                    patch_label = "none",
                    patch_size = 3,
                    ...) {

  if (!inherits(x, what = "brnet")) stop("x must be a 'brnet' object")

  # patch attributes --------------------------------------------------------

  adj <- igraph::graph.adjacency(x$adjacency_matrix,
                                 mode = "undirected")

  patch_attr <- x$df_patch

  ## patch color
  if (patch_color == "env") {

    igraph::V(adj)$patch_value <- patch_attr$environment
    if (is.null(color_label)) color_label <- "Environment"

  }

  if (patch_color == "disturb") {

    igraph::V(adj)$patch_value <- patch_attr$disturbance
    if (is.null(color_label)) color_label <- "Disturbance"

  }

  if (patch_color == "other") {

    if (is.null(value_col)) stop("Provide 'value_col' if patch_color = 'other'")
    igraph::V(adj)$patch_value <- unlist(patch_attr[, value_col])
    if (is.null(color_label)) color_label <- "Value"

  }

  ## edge color
  if (edge_color == "passability") {

    edge_attr <- x$df_edge
    igraph::E(adj)$weight <- patch_attr$passability

  } else {

    if (is.null(edge_color)) edge_color <- grey(0.5)

  }

  ## patch label
  if (patch_label == "patch") {

    igraph::V(adj)$patch_id <- patch_attr$patch_id

  } else {

    if (patch_label == "branch") {

      igraph::V(adj)$patch_id <- patch_attr$branch_id

    } else {

      if (patch_label == "n_upstream") {

        igraph::V(adj)$patch_id <- patch_attr$n_patch_upstream

      } else {

        if (patch_label != "none") message("patch_label should be 'patch', 'branch', or 'n_upstream' to display")
        igraph::V(adj)$patch_id <- ""

      }

    }

  }


  # plot --------------------------------------------------------------------

  g <- ggraph::ggraph(adj,
                      layout = igraph::layout_as_tree(adj,
                                                      root = 1,
                                                      flip.y = FALSE))

  if (!(patch_color %in% c("env", "disturb", "other"))) {
  ## single patch color
    g <- g +
      ggraph::geom_edge_link(color = grey(0.5)) +
      ggraph::geom_node_point(size = patch_size,
                              color = patch_color) +
      ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                              fill = NA,
                              size = 3,
                              label.size = 0,
                              ...) +
      ggplot2::theme_void()

  } else {

    if (is.null(edge_color)) {
      ## variable patch color with single edge color
      g <- g +
        ggraph::geom_edge_link(color = edge_color) +
        ggraph::geom_node_point(ggplot2::aes(color = .data$patch_value),
                                size = patch_size) +
        ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                                fill = NA,
                                size = 3,
                                label.size = 0,
                                ...) +
        MetBrewer::scale_color_met_c("Hiroshige",
                                     direction = -1) +
        ggplot2::labs(color = color_label) +
        ggplot2::theme_void()

    } else {
      ## variable patch & edge color
      g <- g +
        ggraph::geom_edge_link(ggplot2::aes(color = .data$weight)) +
        ggraph::geom_node_point(ggplot2::aes(color = .data$patch_value),
                                size = patch_size) +
        ggraph::geom_node_label(ggplot2::aes(label = .data$patch_id),
                                fill = NA,
                                size = 3,
                                label.size = 0,
                                ...) +
        MetBrewer::scale_color_met_c("Hiroshige",
                                     direction = -1) +
        ggplot2::labs(color = color_label) +
        ggplot2::theme_void()

    }

  }

  return(g)

}

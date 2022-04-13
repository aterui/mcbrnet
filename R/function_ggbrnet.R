#' Visualize a branching network
#'
#' @param adjacency_matrix Adjacency matrix.
#' @param patch_attr Data frame for patch attributes.
#' @param patch_color Type of patch (vertex) label (either \code{"env"}, \code{"disturbance"} or any color code). Default \code{"env"}.
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

ggbrnet <- function(adjacency_matrix,
                    patch_attr,
                    patch_color = "env",
                    patch_label = "none",
                    patch_size = 3,
                    ...) {

  # patch attributes --------------------------------------------------------

  adj <- igraph::graph.adjacency(adjacency_matrix,
                                 mode = "undirected")

  ## patch color
  if (patch_color == "env") {

    igraph::V(adj)$patch_value <- patch_attr$environment
    color_label <- "Environment"

  } else {

    if (patch_color == "disturb") {

      igraph::V(adj)$patch_value <- patch_attr$disturbance
      color_label <- "Disturbance"

    }

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

  if (!(patch_color %in% c("env", "disturb"))) {

    g <- ggraph::ggraph(adj,
                        layout = igraph::layout_as_tree(adj,
                                                        root = 1,
                                                        flip.y = FALSE)) +
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

    g <- ggraph::ggraph(adj,
                        layout = igraph::layout_as_tree(adj,
                                                        root = 1,
                                                        flip.y = FALSE)) +
      ggraph::geom_edge_link(color = grey(0.5)) +
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

  return(g)

}

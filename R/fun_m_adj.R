#' Internal function: adjacency matrix in branching networks
#'
#' @param n_patch number of patches
#' @param p_branch branching probability
#' @param n_branch number of branches
#' @param n_patch_free constraint on number of patches
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_m_adj <- function(n_patch,
                      p_branch,
                      n_branch,
                      n_patch_free) {

  # adjacency matrix: linear network ----------------------------------------

  if (p_branch == 0 | n_branch == 1) {

    v_n_patch_branch <- n_patch
    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    x <- seq_len(n_patch - 1)
    y <- 2:n_patch
    m_adj[cbind(x, y)] <- 1
    m_adj[cbind(y, x)] <- 1

  }


  # adjacency matrix: branched network --------------------------------------

  if (p_branch > 0 & n_branch > 1) {

    # vector of the number of patches in each branch
    if (n_patch_free == FALSE) {
      repeat{

        repeat{

          v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
          if (sum(v_n_patch_branch) >= n_patch) break

        }

        if (sum(v_n_patch_branch) == n_patch) break

      }
    } else {

      v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
      n_patch <- sum(v_n_patch_branch)

    }

    # start_id, end_id, and neighbor list for each branch
    v_end_id <- cumsum(v_n_patch_branch)
    v_start_id <- v_end_id - (v_n_patch_branch - 1)

    list_neighbor_inbranch <- lapply(seq_len(n_branch),
                                     function(i) {
                                       fun_adj(n = v_n_patch_branch[i],
                                               start_id = v_start_id[i])
                                     }
    )

    m_neighbor_inbranch <- do.call(rbind, list_neighbor_inbranch)

    # combine parent and offspring branches at confluences
    if (n_branch == 3) {

      parent <- c(1, 1)
      offspg <- c(2, 3)
      m_po <- cbind(parent, offspg)

      list_confluence <- lapply(seq_len(nrow(m_po)),
                                function(x) cbind(v_end_id[m_po[x, 1]],
                                                  v_start_id[m_po[x, 2]]))

      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])

    } else {

      n_confluence <- 0.5 * (n_branch - 1)
      v_parent_branch <- seq_len(n_confluence)
      v_offspg_branch <- 2:n_branch

      m_offspg <- matrix(NA,
                         nrow = 2,
                         ncol = n_confluence)

      for (i in n_confluence:1) {

        v_y <- resample(v_offspg_branch[v_offspg_branch > v_parent_branch[i]],
                        size = 2)
        v_offspg_branch <- setdiff(v_offspg_branch, v_y)
        m_offspg[, i] <- v_y

      }

      parent <- rep(v_parent_branch, each = 2)
      offspg <- c(m_offspg)
      m_po <- cbind(parent, offspg)

      list_confluence <- lapply(seq_len(nrow(m_po)),
                                function(x) cbind(v_end_id[m_po[x, 1]],
                                                  v_start_id[m_po[x, 2]]))

      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])

    }

    m_neighbor_patch <- rbind(m_neighbor_inbranch, m_confluence)
    m_neighbor_patch <- m_neighbor_patch[complete.cases(m_neighbor_patch), ]

    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    m_adj[m_neighbor_patch] <- 1

  }

  return(list(v_n_patch_branch = v_n_patch_branch,
              m_adj = m_adj))
}

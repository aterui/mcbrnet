
# setup -------------------------------------------------------------------

n_patch <- round(runif(1, 20, 50))
p_branch <- runif(1, 0.2, 0.8)

n_branch <- rbinom(1, n_patch, p_branch)
while(n_branch %% 2 == 0) {
  n_branch <- rbinom(1, n_patch, p_branch)
}

list_adj <- fun_m_adj(n_patch = n_patch,
                      p_branch = p_branch,
                      n_branch = n_branch,
                      n_patch_free = FALSE)

n_conf_raw <- sum(rowSums(list_adj$m_adj) == 3)
n_conf <- ifelse(sum(list_adj$m_adj[1,]) == 2,
                 n_conf_raw + 1,
                 n_conf_raw)

n_source_raw <- sum(rowSums(list_adj$m_adj) == 1)
n_source <- ifelse(sum(list_adj$m_adj[1,]) == 2,
                   n_source_raw,
                   n_source_raw - 1)


# test --------------------------------------------------------------------

test_that("adjacency degree is less than 4 for all patches", {
  expect_true(any(rowSums(list_adj$m_adj) < 4))
})

test_that("number of confluences", {
  expect_equal(n_conf,
               (n_branch - 1) * 0.5)
})

test_that("number of source streams", {
  expect_equal(n_source,
               (n_branch + 1) * 0.5)
})


# setup -------------------------------------------------------------------

## param
n_species <- round(runif(1, 5, 15))
alpha <- runif(1)
min_alpha <- runif(1, 0, 0.4)
max_alpha <- runif(1, 0.5, 1)

## sample
m_int_c <- fun_int_mat(n_species = n_species,
                       alpha = alpha,
                       min_alpha = min_alpha,
                       max_alpha = max_alpha,
                       interaction_type = "constant")

diag(m_int_c) <- NA
v_alpha_c <- m_int_c[!is.na(m_int_c)]

m_int_r <- fun_int_mat(n_species = n_species,
                       alpha = alpha,
                       min_alpha = min_alpha,
                       max_alpha = max_alpha,
                       interaction_type = "random")

diag(m_int_r) <- NA
v_alpha_r <- m_int_r[!is.na(m_int_r)]


# test --------------------------------------------------------------------

test_that("multiplication works", {

  expect_equal(unique(max(v_alpha_c)),
               alpha)
  expect_lt(max(v_alpha_r),
            max_alpha)
  expect_gt(min(v_alpha_r),
            min_alpha)

})

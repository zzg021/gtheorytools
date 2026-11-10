# REML estimates and their observed-information covariance against a REML fit
# from n x n matrices.

# REML by Fisher scoring on the variance scale, and the inverse observed
# information at the estimates. With V_k = Z_k Z_k', V = sum_k sigma^2_k V_k,
# and P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1 for the grand mean X = 1, the
# REML score is (y'P V_k P y - tr(P V_k)) / 2, the expected information
# tr(P V_k P V_l) / 2, and the observed information
# y'P V_k P V_l P y - tr(P V_k P V_l) / 2. The fit starts with every
# component at the variance of the scores over the number of components.
# `effects` lists the data columns of each effect, as in dense_anova().
dense_reml <- function(dat, effects) {
  dat <- dat[!is.na(dat$score), ]
  n <- nrow(dat)
  q <- length(effects)
  indicator <- function(columns) {
    cells <- interaction(dat[columns], drop = TRUE)
    return(outer(as.integer(cells), seq_len(nlevels(cells)), `==`) * 1)
  }
  v_matrices <- lapply(effects, function(columns) tcrossprod(indicator(columns)))
  x <- matrix(1, n, 1)
  projection <- function(theta) {
    v_inverse <- solve(Reduce(`+`, Map(`*`, theta, v_matrices)))
    v_inverse_x <- v_inverse %*% x
    return(v_inverse - v_inverse_x %*% solve(t(x) %*% v_inverse_x, t(v_inverse_x)))
  }

  theta <- rep(var(dat$score) / q, q)
  for (iteration in 1:100) {
    p_matrix <- projection(theta)
    py <- p_matrix %*% dat$score
    pv <- lapply(v_matrices, function(v) p_matrix %*% v)
    score <- numeric(q)
    expected_information <- matrix(0, q, q)
    for (k in 1:q) {
      score[k] <- (sum(py * (v_matrices[[k]] %*% py)) - sum(diag(pv[[k]]))) / 2
      for (l in 1:q) {
        expected_information[k, l] <- sum(pv[[k]] * t(pv[[l]])) / 2
      }
    }
    step <- solve(expected_information, score)
    theta <- theta + step
    if (max(abs(step)) < 1e-12 * max(theta)) {
      break
    }
  }

  p_matrix <- projection(theta)
  py <- p_matrix %*% dat$score
  pv <- lapply(v_matrices, function(v) p_matrix %*% v)
  observed_information <- matrix(0, q, q)
  for (k in 1:q) {
    for (l in 1:q) {
      observed_information[k, l] <- as.numeric(t(py) %*% v_matrices[[k]] %*% pv[[l]] %*% py) -
        sum(pv[[k]] * t(pv[[l]])) / 2
    }
  }
  return(list(estimates = theta, vcov = solve(observed_information)))
}

# lmer's estimates agree to its convergence tolerance. The covariance comes
# from a numerical Hessian on the SD scale; at the REML estimates its delta-
# method map equals the inverse observed information on the variance scale,
# which differs from the inverse expected information by 1 to 9 percent in
# these data.
expect_dense_reml <- function(g, dat, effects) {
  expect_identical(names(g$variance_component), names(effects))
  dense <- dense_reml(dat, effects)
  expect_equal(unname(g$variance_component), dense$estimates, tolerance = 1e-5)
  expect_equal(unname(g$variance_component_vcov), dense$vcov, tolerance = 1e-4)
}

test_that("REML agrees with the dense calculation for Brennan Table 7.2 (i:p)", {
  dat <- brennan_i_in_p_unbalanced
  g <- g_study_general(dat, design = "i:p")
  expect_dense_reml(g, dat, effects_i_in_p)
})

test_that("REML agrees with the dense calculation for Brennan Table 7.7 (p x i, missing scores)", {
  dat <- brennan_p_x_i_missing
  g <- g_study_general(dat, design = "p x i")
  expect_dense_reml(g, dat, effects_p_x_i)
})

test_that("REML agrees with the dense calculation for urGENOVA's p x (i:h) with empty cells", {
  dat <- urgenova_p_x_i_in_h()
  g <- g_study_general(dat, design = "p x (i:h)")
  expect_dense_reml(g, dat, effects_p_x_i_in_h)
})

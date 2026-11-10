# ANOVA estimates against values printed by Brennan (2001) and GENOVA. Each
# tolerance is half a unit in the last printed decimal. Expected values are
# named by effect, in the package's order.

expect_printed <- function(object, expected, tolerance) {
  if (!is.null(names(expected))) {
    expect_identical(names(object), names(expected))
  }
  expect_lt(max(abs(unname(object) - unname(expected))), tolerance)
}

# A column of the ANOVA table, named by effect.
anova_column <- function(g, column) {
  result <- g$anova_table[[column]]
  names(result) <- g$anova_table$effect
  return(result)
}

# Analogous ANOVA estimates and their unbiased normal-theory covariance from
# n x n matrices. `effects` lists the data columns of each effect in the
# package's order, the residual last. SS(a) = y'Q_a y, where Q_a combines the
# projections on cell means; df(a) = tr(Q_a); E SS(a) = sum_k sigma^2_k
# tr(Q_a Z_k Z_k'); and Cov(SS(a), SS(b)) = 2 tr(Q_a V Q_b V).
dense_anova <- function(dat, effects) {
  dat <- dat[!is.na(dat$score), ]
  n <- nrow(dat)
  q <- length(effects)
  indicator <- function(columns) {
    cells <- interaction(dat[columns], drop = TRUE)
    return(outer(as.integer(cells), seq_len(nlevels(cells)), `==`) * 1)
  }
  projection <- function(columns) {
    cell_indicator <- indicator(columns)
    return(cell_indicator %*% diag(1 / colSums(cell_indicator), ncol(cell_indicator)) %*%
      t(cell_indicator))
  }
  q_matrices <- list()
  for (a in order(lengths(effects))) {
    lower <- names(q_matrices)[vapply(
      names(q_matrices),
      function(b) all(effects[[b]] %in% effects[[a]]),
      logical(1)
    )]
    q_matrices[[names(effects)[a]]] <- projection(effects[[a]]) - matrix(1 / n, n, n) -
      Reduce(`+`, q_matrices[lower], matrix(0, n, n))
  }
  q_matrices <- q_matrices[names(effects)]
  z_matrices <- lapply(effects, function(columns) tcrossprod(indicator(columns)))

  df <- vapply(q_matrices, function(m) sum(diag(m)), numeric(1))
  ems <- matrix(0, q, q)
  qz <- list()
  for (a in 1:q) {
    for (k in 1:q) {
      qz[[paste(a, k)]] <- q_matrices[[a]] %*% z_matrices[[k]]
      ems[a, k] <- sum(diag(qz[[paste(a, k)]])) / df[[a]]
    }
  }
  ss <- vapply(q_matrices, function(m) as.numeric(t(dat$score) %*% m %*% dat$score), numeric(1))
  weights <- solve(ems) %*% diag(1 / df)
  estimates <- as.vector(weights %*% ss)

  g_matrix <- matrix(0, q^2, q^2)
  for (k in 1:q) {
    for (l in 1:q) {
      ss_covariance <- matrix(0, q, q)
      for (a in 1:q) {
        for (b in 1:q) {
          ss_covariance[a, b] <- 2 * sum(qz[[paste(a, k)]] * t(qz[[paste(b, l)]]))
        }
      }
      g_matrix[, (l - 1) * q + k] <- as.vector(weights %*% ss_covariance %*% t(weights))
    }
  }
  products <- solve(diag(q^2) + g_matrix, as.vector(estimates %o% estimates))
  vcov <- matrix(g_matrix %*% products, q, q)
  return(list(df = df, ems = ems, estimates = estimates, vcov = (vcov + t(vcov)) / 2))
}

expect_dense_agreement <- function(g, dat, effects) {
  expect_identical(g$anova_table$effect, names(effects))
  dense <- dense_anova(dat, effects)
  expect_equal(g$anova_table$df, unname(dense$df))
  expect_equal(unname(g$ems_coefficients), dense$ems)
  expect_equal(unname(g$anova_table$estimate), dense$estimates)
  expect_equal(unname(g$variance_component_vcov), dense$vcov)
}

# Error variances of the D study with the observed pattern, from n x n
# matrices: with M the matrix of the means of the objects of measurement,
# V = sum_k sigma^2_k Z_k Z_k', and C the centering matrix over objects,
# sigma^2(Delta) is the mean over objects of Var(Xbar_o - mu_o) =
# (M V M')_oo - sigma^2(o), and sigma^2(delta) is
# E S^2(o) - sigma^2(o) = tr(C M V M' C) / (n_o - 1) - sigma^2(o). `effects`
# is as in dense_anova(), with the object of measurement first.
dense_errors <- function(dat, effects, components) {
  dat <- dat[!is.na(dat$score), ]
  indicator <- function(columns) {
    cells <- interaction(dat[columns], drop = TRUE)
    return(outer(as.integer(cells), seq_len(nlevels(cells)), `==`) * 1)
  }
  objects <- indicator(effects[[1]])
  object_mean <- t(objects) / colSums(objects)
  score_covariance <- matrix(0, nrow(dat), nrow(dat))
  for (k in seq_along(effects)) {
    score_covariance <- score_covariance + components[[k]] * tcrossprod(indicator(effects[[k]]))
  }
  mean_covariance <- object_mean %*% score_covariance %*% t(object_mean)
  n_objects <- ncol(objects)
  centering <- diag(n_objects) - 1 / n_objects
  Delta <- mean(diag(mean_covariance)) - components[[1]]
  delta <- sum(diag(centering %*% mean_covariance %*% centering)) / (n_objects - 1) -
    components[[1]]
  return(c(delta = delta, Delta = Delta))
}

d_estimates <- function(d) {
  result <- d$coefficient_results$estimate
  names(result) <- d$coefficient_results$quantity
  return(result)
}

expect_dense_d_study <- function(g, dat, effects) {
  facets <- setdiff(names(g$parents), "p")
  universe_sizes <- rep(Inf, length(facets))
  names(universe_sizes) <- facets
  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = universe_sizes)
  dense <- dense_errors(dat, effects, unname(g$variance_component[names(effects)]))
  expect_equal(d_estimates(d)[c("delta", "Delta")], dense)
}

test_that("ANOVA reproduces GENOVA for p x (r:t) (Brennan Appendix F)", {
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", method = "anova")

  # The T terms are in Table 3.4.
  expect_printed(
    anova_column(g, "t_term"),
    c(p = 2800.1667, t = 2755.7000, "r:t" = 2835.4000, pt = 2931.5000, "pr:t" = 3204.0000),
    5e-5
  )
  expect_printed(
    anova_column(g, "sum_of_squares"),
    c(p = 92.66667, t = 48.20000, "r:t" = 79.70000, pt = 83.13333, "pr:t" = 192.80000),
    5e-6
  )
  expect_printed(
    anova_column(g, "mean_square"),
    c(p = 10.29630, t = 24.10000, "r:t" = 8.85556, pt = 4.61852, "pr:t" = 2.38025),
    5e-6
  )
  expect_equal(
    g$ems_coefficients["t", ], c(p = 0, t = 40, "r:t" = 10, pt = 4, "pr:t" = 1)
  )
  expect_printed(
    g$variance_component,
    c(p = 0.4731481, t = 0.3251543, "r:t" = 0.6475309, pt = 0.5595679, "pr:t" = 2.3802469),
    5e-8
  )
  expect_printed(
    g$variance_component_se,
    c(p = 0.3855758, t = 0.4379875, "r:t" = 0.3794056, pt = 0.3766291, "pr:t" = 0.3694860),
    5e-8
  )

  # GENOVA's covariance matrix, in its order p, t, r:t, pt, pr:t.
  genova_order <- c("p", "t", "r:t", "pt", "pr:t")
  genova_vcov <- matrix(0, 5, 5, dimnames = list(genova_order, genova_order))
  genova_vcov[lower.tri(genova_vcov, diag = TRUE)] <- c(
    0.1486687, 0.0044439, 0.0000000, -0.0444390, 0.0000000,
    0.1918331, -0.0359871, -0.0141849, 0.0034130,
    0.1439486, 0.0034130, -0.0136520,
    0.1418495, -0.0341300,
    0.1365199
  )
  genova_vcov <- genova_vcov + t(genova_vcov) - diag(diag(genova_vcov))
  package_order <- names(g$variance_component)
  expect_printed(
    g$variance_component_vcov, genova_vcov[package_order, package_order], 5e-8
  )

  # On balanced data, the covariance used for unbalanced data is Eq. 6.1.
  dat <- .prepare_g_study_data(
    brennan_p_x_r_in_t, .build_g_study_design("p x (r:t)", NULL, "score")
  )
  general_vcov <- .anova_covariance_unbalanced(
    dat, package_order, g$anova_table$df, g$ems_coefficients, g$anova_table$estimate
  )
  expect_equal(unname(general_vcov), unname(g$variance_component_vcov))

  # D study #1: 3 tasks and 4 raters per task, both random. GENOVA prints SEs
  # for tau, delta, and Delta.
  d <- d_study_general(
    g, object_of_measurement = "p",
    universe_sizes = c(r = Inf, t = Inf), d_sizes = c(r = 4, t = 3)
  )
  expect_printed(
    d$coefficient_results$estimate, c(0.47315, 0.38488, 0.54722, 0.55144, 0.46370), 5e-6
  )
  expect_printed(d$coefficient_results$standard_error[1:3], c(0.38558, 0.12171, 0.17935), 5e-6)
})

test_that("ANOVA estimates below zero are set to zero (Brennan Table 3.5)", {
  # Synthetic Data Set No. 3 analyzed as p x (i:o), as if different items were
  # used on the two occasions.
  g <- g_study_general(brennan_p_x_i_x_o, design = "p x (i:o)", method = "anova")

  expect_printed(
    anova_column(g, "mean_square"),
    c(p = 6.9111, o = 3.2000, "i:o" = 7.4917, po = 1.3389, "pi:o" = 1.5102),
    5e-5
  )
  expect_printed(
    anova_column(g, "estimate"),
    c(p = 0.6965, o = -0.1030, "i:o" = 0.5981, po = -0.0428, "pi:o" = 1.5102),
    5e-5
  )
  expect_printed(
    g$variance_component,
    c(p = 0.6965, o = 0, "i:o" = 0.5981, po = 0, "pi:o" = 1.5102),
    5e-5
  )
})

test_that("ANOVA reproduces Brennan Tables 2.3, 2.9, and 3.3", {
  g_p_x_i <- g_study_general(brennan_p_x_i, design = "p x i", method = "anova")
  expect_printed(g_p_x_i$variance_component, c(p = 0.0574, i = 0.0754, pi = 0.1269), 5e-5)

  g_i_in_p <- g_study_general(brennan_i_in_p, design = "i:p", method = "anova")
  expect_printed(g_i_in_p$variance_component, c(p = 0.6108, "i:p" = 2.0250), 5e-5)

  g_p_x_i_x_o <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o", method = "anova")
  expect_printed(
    g_p_x_i_x_o$variance_component,
    c(p = 0.5528, i = 0.4417, o = 0.0074, pi = 0.5750, po = 0.1009, io = 0.1565, pio = 0.9352),
    5e-5
  )
})

test_that("ANOVA with tasks fixed reproduces GENOVA D study #3 (Brennan Appendix F)", {
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t", method = "anova")
  reported <- g$variance_components$variance_estimate
  names(reported) <- g$variance_components$effect
  expect_printed(
    reported[c("p", "r:t", "pt", "pr:t")],
    c(p = 0.65967, "r:t" = 0.64753, pt = 0.55957, "pr:t" = 2.38025),
    5e-6
  )

  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = c(r = Inf), d_sizes = c(r = 4))
  expect_printed(
    d$coefficient_results$estimate, c(0.65967, 0.19835, 0.25231, 0.76882, 0.72333), 5e-6
  )
  expect_printed(d$coefficient_results$standard_error[1:3], c(0.36716, 0.03079, 0.04193), 5e-6)
})

test_that("ANOVA with b and c fixed reproduces GENOVA's test2", {
  # GENOVA manual, pp. B-8 to B-13: S:(A x B x C) with B and C fixed, 2 levels
  # of a, b, and c, and 4 s in each abc cell, read with s varying fastest.
  # Eq. 3.34 uses sigma^2(a) = -.625 as computed (p. B-12).
  dat <- expand.grid(s = 1:4, c = 1:2, b = 1:2, a = 1:2)
  dat$score <- c(
    3, 6, 3, 3, 4, 5, 4, 3, 7, 8, 7, 6, 7, 8, 9, 8,
    1, 2, 2, 2, 2, 3, 4, 3, 5, 6, 5, 6, 10, 10, 9, 11
  )
  g <- g_study_general(dat, design = "s:(a x b x c)", fixed = c("b", "c"), method = "anova")

  expect_printed(anova_column(g, "estimate")[["a"]], -0.625, 5e-8)
  reported <- g$variance_components$variance_estimate
  names(reported) <- g$variance_components$effect
  expect_printed(
    reported[c("a", "ab", "ac", "abc", "s:abc")],
    c(a = 0.1471354, ab = 0.6692708, ac = 1.1692708, abc = 0.5885417, "s:abc" = 0.7708333),
    5e-8
  )
})

test_that("ANOVA reproduces Brennan Table 7.2 (unbalanced i:p)", {
  dat <- brennan_i_in_p_unbalanced
  g <- g_study_general(dat, design = "i:p", method = "anova")

  expect_printed(anova_column(g, "t_term"), c(p = 1069.40, "i:p" = 1120.00), 5e-3)
  expect_printed(anova_column(g, "mean_square"), c(p = 3.8912, "i:p" = 1.8071), 5e-5)
  expect_printed(g$variance_component, c(p = 0.5098, "i:p" = 1.8071), 5e-5)
  expect_dense_agreement(g, dat, effects_i_in_p)

  # The D study with the G-study pattern (Sec. 7.2.1): n~_i = 3.9130, the
  # harmonic mean of the n_i:p (Eq. 7.19), sigma^2(Delta) = sigma^2(delta) =
  # 1.8071 / 3.9130 = .4618 (Eq. 7.18), and E rho^2 = .525.
  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf))
  expect_printed(d$d_mean_results$divisor[[2]], 3.9130, 5e-5)
  estimates <- d_estimates(d)
  expect_printed(estimates[c("delta", "Delta")], c(delta = 0.4618, Delta = 0.4618), 5e-5)
  expect_printed(estimates[["G"]], 0.525, 5e-4)
  expect_dense_d_study(g, dat, effects_i_in_p)
})

test_that("ANOVA reproduces Brennan Table 7.5 and Eq. 7.13 (unbalanced p x (i:h))", {
  dat <- brennan_p_x_i_in_h_unbalanced
  g <- g_study_general(dat, design = "p x (i:h)", method = "anova")

  expect_printed(
    anova_column(g, "t_term"),
    c(p = 1390.0, h = 1424.0, "i:h" = 1440.0, ph = 1564.5, "pi:h" = 1610.0),
    0.05
  )
  expect_printed(
    anova_column(g, "mean_square"),
    c(p = 13.4286, h = 64.0000, "i:h" = 3.2000, ph = 3.3214, "pi:h" = 0.8429),
    5e-5
  )
  # Brennan prints .9913 for sigma^2(ph); Eq. 7.14 with his mean squares gives
  # (3.3214 - .8429) / 2.5 = .9914.
  expect_printed(
    g$variance_component,
    c(p = 1.2014, h = 2.9161, "i:h" = 0.2946, ph = 0.9914, "pi:h" = 0.8429),
    5e-5
  )

  # Eq. 7.13 with n_p = 8, n_i+ = 8, r_i = 24 / 8 = 3, and t_i = (8 - 3) / 2.
  effects <- c("p", "h", "i:h", "ph", "pi:h")
  expected_ems <- rbind(
    c(8, 0, 0, 3, 1),
    c(0, 8 * 2.5, 8, 2.5, 1),
    c(0, 0, 8, 0, 1),
    c(0, 0, 0, 2.5, 1),
    c(0, 0, 0, 0, 1)
  )
  dimnames(expected_ems) <- list(effects, effects)
  expect_equal(g$ems_coefficients, expected_ems)
  expect_dense_agreement(g, dat, effects_p_x_i_in_h)

  # The D study with the G-study pattern (Sec. 7.2.2): n~_h = 64 / 24 and
  # n_i+ = 8 (Eq. 7.28), sigma^2(Delta) = 1.6075 (Eq. 7.27), sigma^2(delta)
  # = .4771 (Eq. 7.31), and E rho^2 = .716 (Eq. 7.32).
  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf, h = Inf))
  divisor <- d$d_mean_results$divisor
  names(divisor) <- d$d_mean_results$effect
  expect_equal(divisor, c(p = 1, h = 64 / 24, "i:h" = 8, ph = 64 / 24, "pi:h" = 8))
  estimates <- d_estimates(d)
  expect_printed(estimates[c("delta", "Delta")], c(delta = 0.4771, Delta = 1.6075), 5e-5)
  expect_printed(estimates[["G"]], 0.716, 5e-4)
  # The standard error of sigma^2(Delta) from the G-study covariance with the
  # weights 1 / n~_h and 1 / n_i+.
  weights <- c(0, 24 / 64, 1 / 8, 24 / 64, 1 / 8)
  expect_equal(
    d$coefficient_results$standard_error[[3]],
    sqrt(as.numeric(t(weights) %*% g$variance_component_vcov %*% weights))
  )
  expect_dense_d_study(g, dat, effects_p_x_i_in_h)
})

test_that("ANOVA reproduces urGENOVA for p x i with missing scores (Brennan Table 7.7)", {
  dat <- brennan_p_x_i_missing
  g <- g_study_general(dat, design = "p x i", method = "anova")

  # urGENOVA manual, Appendix F (Run 5).
  expect_equal(g$anova_table$df, c(11, 5, 42))
  expect_printed(
    anova_column(g, "t_term"), c(p = 27.80000, i = 24.74318, pi = 37.00000), 5e-6
  )
  expect_printed(
    anova_column(g, "mean_square"), c(p = 0.41787, i = 0.30796, pi = 0.18239), 5e-6
  )
  expect_printed(g$variance_component, c(p = 0.04725, i = 0.01174, pi = 0.18401), 5e-6)
  urgenova_ems <- rbind(
    c(4.903, 0.186, 1),
    c(0.186, 9.810, 1),
    c(-0.022, -0.049, 1)
  )
  expect_printed(g$ems_coefficients, urgenova_ems, 5e-4)
  expect_dense_agreement(g, dat, effects_p_x_i)

  # The D study with the G-study pattern (Sec. 7.2.4): n~ = 4.7682, the
  # harmonic mean of the n_p (Eq. 7.42), and sigma^2(Delta) = (sigma^2(i) +
  # sigma^2(pi)) / n~ (Eq. 7.41). Brennan's .0410 uses the components rounded
  # to four decimals; urGENOVA's give .04105.
  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf))
  expect_printed(d$d_mean_results$divisor[2:3], c(4.7682, 4.7682), 5e-5)
  observed <- dat[!is.na(dat$score), ]
  took_item <- unclass(table(observed$p, observed$i))
  n_items <- rowSums(took_item)
  harmonic_mean <- length(n_items) / sum(1 / n_items)
  vc <- g$variance_component
  estimates <- d_estimates(d)
  expect_equal(estimates[["Delta"]], (vc[["i"]] + vc[["pi"]]) / harmonic_mean)
  # sigma^2(delta) from the last two terms of Eq. 7.45: the coefficient of
  # sigma^2(i) is 1 / n~ less the mean, over pairs of persons, of
  # n_pp' / (n_p n_p'), with n_pp' the number of items both took. This gives
  # .0391 and E rho^2 = .547.
  pair_share <- tcrossprod(took_item) / outer(n_items, n_items)
  diag(pair_share) <- 0
  n_persons <- length(n_items)
  coefficient_i <- 1 / harmonic_mean - sum(pair_share) / (n_persons * (n_persons - 1))
  expect_equal(
    estimates[["delta"]], coefficient_i * vc[["i"]] + vc[["pi"]] / harmonic_mean
  )
  expect_dense_d_study(g, dat, effects_p_x_i)

  # With items as the object of measurement, the pattern weights each item's
  # scores by person.
  d_items <- d_study_general(g, object_of_measurement = "i", universe_sizes = c(p = Inf))
  dense <- dense_errors(dat, effects_p_x_i[c("i", "p", "pi")], unname(vc[c("i", "p", "pi")]))
  expect_equal(d_estimates(d_items)[c("delta", "Delta")], dense)
})

test_that("ANOVA reproduces urGENOVA for p x (i:h) with empty cells", {
  dat <- urgenova_p_x_i_in_h()
  g <- g_study_general(dat, design = "p x (i:h)", method = "anova")

  # urGENOVA manual, Appendix G (Run 6).
  expect_equal(g$anova_table$df, c(9, 2, 9, 18, 75))
  expect_printed(
    anova_column(g, "t_term"),
    c(p = 2390.53788, h = 2360.44089, "i:h" = 2436.26667, ph = 2512.00000, "pi:h" = 2787.00000),
    5e-6
  )
  expect_printed(
    anova_column(g, "mean_square"),
    c(p = 9.11532, h = 25.97044, "i:h" = 8.42509, ph = 3.86229, "pi:h" = 2.65566),
    5e-6
  )
  expect_printed(
    g$variance_component,
    c(p = 0.43859, h = 0.45415, "i:h" = 0.60321, ph = 0.32937, "pi:h" = 2.67060),
    5e-6
  )
  # urGENOVA's order p, h, ph, i:h, pi:h.
  urgenova_order <- c("p", "h", "ph", "i:h", "pi:h")
  urgenova_ems <- rbind(
    c(11.398, 0.031, 4.254, 0.051, 1),
    c(0.042, 35.868, 3.621, 9.614, 1),
    c(-0.005, -0.016, 3.567, 0.043, 1),
    c(0.055, 0, 0.055, 9.470, 1),
    c(-0.007, 0, -0.007, -0.016, 1)
  )
  dimnames(urgenova_ems) <- list(urgenova_order, urgenova_order)
  package_order <- names(g$variance_component)
  expect_printed(g$ems_coefficients, urgenova_ems[package_order, package_order], 5e-4)
  expect_dense_agreement(g, dat, effects_p_x_i_in_h)
  expect_dense_d_study(g, dat, effects_p_x_i_in_h)
})

test_that("ANOVA agrees with the closed form for unbalanced i:h:p (Brennan Appendix E)", {
  set.seed(7)
  n_h_per_p <- c(2, 3, 2, 3, 3, 2)
  dat <- data.frame(p = rep(seq_along(n_h_per_p), n_h_per_p))
  dat$h <- sequence(n_h_per_p)
  n_i_per_ph <- sample(1:4, nrow(dat), replace = TRUE)
  dat <- dat[rep(seq_len(nrow(dat)), n_i_per_ph), ]
  dat$i <- sequence(n_i_per_ph)
  dat$score <- round(rnorm(nrow(dat), mean = 5, sd = 2))
  g <- g_study_general(dat, design = "i:h:p", method = "anova")

  t_term <- function(cells) {
    return(sum(tapply(dat$score, cells, sum)^2 / tapply(dat$score, cells, length)))
  }
  n_plus <- nrow(dat)
  n_h_plus <- sum(n_h_per_p)
  n_p <- length(n_h_per_p)
  n_ph <- as.vector(table(dat$h, dat$p)[table(dat$h, dat$p) > 0])
  n_p_of_ph <- as.vector(table(dat$p)[rep(seq_along(n_h_per_p), n_h_per_p)])
  r_pi <- sum(n_ph^2 / n_p_of_ph)
  r_i <- sum(n_ph^2) / n_plus
  r_ih <- sum(table(dat$p)^2) / n_plus
  t_mean <- sum(dat$score)^2 / n_plus
  t_p <- t_term(dat$p)
  t_hp <- t_term(paste(dat$p, dat$h))
  t_ihp <- sum(dat$score^2)
  vc_ihp <- (t_ihp - t_hp) / (n_plus - n_h_plus)
  vc_hp <- (t_hp - t_p - (n_h_plus - n_p) * vc_ihp) / (n_plus - r_pi)
  vc_p <- (t_p - t_mean - (n_p - 1) * vc_ihp - (r_pi - r_i) * vc_hp) / (n_plus - r_ih)

  expect_equal(unname(g$anova_table$estimate), c(vc_p, vc_hp, vc_ihp))
  effects <- list(p = "p", "h:p" = c("h", "p"), "i:h:p" = c("i", "h", "p"))
  expect_dense_agreement(g, dat, effects)
  expect_dense_d_study(g, dat, effects)
})

test_that("ANOVA agrees with the dense calculation for unbalanced p x (i:h) x o", {
  # Five persons, 2 and 3 items in the two levels of h, two occasions, and
  # four missing scores: 11 components.
  set.seed(11)
  conditions <- expand.grid(i = 1:3, h = 1:2, o = 1:2)
  conditions <- conditions[!(conditions$h == 1 & conditions$i == 3), ]
  scores <- matrix(round(rnorm(5 * nrow(conditions), mean = 5, sd = 2)), nrow = 5)
  scores[cbind(c(1, 2, 4, 5), c(2, 7, 4, 10))] <- NA
  dat <- long_scores(scores, conditions)
  g <- g_study_general(dat, design = "p x (i:h) x o", method = "anova")
  effects <- list(
    p = "p", h = "h", o = "o", "i:h" = c("i", "h"), ph = c("p", "h"), po = c("p", "o"),
    ho = c("h", "o"), "pi:h" = c("p", "i", "h"), "io:h" = c("i", "o", "h"),
    pho = c("p", "h", "o"), "pio:h" = c("p", "i", "o", "h")
  )
  expect_dense_agreement(g, dat, effects)
  expect_dense_d_study(g, dat, effects)
})

test_that("ANOVA and REML estimates agree for balanced p x (a:b) x (c:d:b)", {
  # Every ANOVA estimate is above zero, so the REML estimates of these
  # balanced data solve the same equations; they agree to lmer's convergence
  # tolerance.
  design <- "p x (a:b) x (c:d:b)"
  set.seed(15)
  dat <- expand.grid(p = 1:12, a = 1:2, b = 1:3, c = 1:3, d = 1:2)
  dat$score <- 0
  for (members in .design_effects(.parse_design(design, "design")$parents)) {
    cells <- interaction(dat[members], drop = TRUE)
    dat$score <- dat$score + rnorm(nlevels(cells))[cells]
  }
  g_anova <- g_study_general(dat, design = design, method = "anova")
  g_reml <- g_study_general(dat, design = design)
  expect_equal(g_anova$variance_component, g_reml$variance_component, tolerance = 1e-4)
})

test_that("ANOVA requires positive degrees of freedom for every effect", {
  # Each person takes different items, so p x i has negative df.
  dat <- data.frame(p = rep(1:3, each = 2), i = 1:6, score = c(3, 5, 4, 6, 2, 7))
  expect_error(
    g_study_general(dat, design = "p x i", method = "anova"),
    "positive degrees of freedom"
  )
})

test_that("A sampling variance estimated below zero gives an NA standard error", {
  scores <- rbind(
    c(4, -1, 7, 5, 3, 7), c(6, 5, 3, 5, 6, 1), c(7, 5, 0, 5, 3, 7), c(5, 3, 6, 7, 5, 4)
  )
  conditions <- data.frame(h = c(1, 1, 1, 2, 3, 3), i = c(1, 2, 3, 1, 1, 2))
  expect_warning(
    g <- g_study_general(long_scores(scores, conditions), design = "p x (i:h)", method = "anova"),
    "sampling variance of p is below zero"
  )
  expect_true(is.na(g$variance_component_se[["p"]]))
  expect_false(anyNA(g$variance_component_se[-1]))
})

test_that("A D-study sampling variance estimated below zero gives an NA standard error", {
  # Every G-study sampling variance is above zero, but the estimated covariance
  # has a negative eigenvalue.
  scores <- rbind(c(7, 6, 5, 8, 6), c(7, 3, 5, 9, 6), c(6, 7, 7, 9, 7), c(6, 5, 9, 3, 7))
  conditions <- data.frame(h = c(1, 2, 3, 3, 3), i = c(1, 1, 1, 2, 3))
  g <- g_study_general(long_scores(scores, conditions), design = "p x (i:h)", method = "anova")
  expect_warning(
    d <- d_study_general(
      g, object_of_measurement = "p",
      universe_sizes = c(i = Inf, h = Inf), d_sizes = c(i = 1, h = 6)
    ),
    "sampling variance of G is below zero"
  )
  standard_errors <- d$coefficient_results$standard_error
  expect_true(is.na(standard_errors[[4]]))
  expect_false(anyNA(standard_errors[-4]))
})

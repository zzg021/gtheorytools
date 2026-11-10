# ANOVA estimates against values printed by Brennan (2001) and GENOVA. Each
# tolerance is half a unit in the last printed decimal.

expect_printed <- function(object, expected, tolerance) {
  expect_lt(max(abs(unname(object) - unname(expected))), tolerance)
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
  dense <- dense_anova(dat, effects)
  expect_equal(g$anova_table$df, unname(dense$df))
  expect_equal(unname(g$ems_coefficients), dense$ems)
  expect_equal(unname(g$anova_table$estimate), dense$estimates)
  expect_equal(unname(g$variance_component_vcov), dense$vcov)
}

# Error variances of the D study with the observed pattern, from n x n
# matrices: with M the matrix of person means, V = sum_k sigma^2_k Z_k Z_k',
# and C the centering matrix over persons, sigma^2(Delta) is the mean over
# persons of Var(Xbar_p - mu_p) = (M V M')_pp - sigma^2(p), and
# sigma^2(delta) is E S^2(p) - sigma^2(p) = tr(C M V M' C) / (n_p - 1) -
# sigma^2(p). `effects` is as in dense_anova(), persons first.
dense_errors <- function(dat, effects, components) {
  dat <- dat[!is.na(dat$score), ]
  indicator <- function(columns) {
    cells <- interaction(dat[columns], drop = TRUE)
    return(outer(as.integer(cells), seq_len(nlevels(cells)), `==`) * 1)
  }
  persons <- indicator(effects[[1]])
  person_mean <- t(persons) / colSums(persons)
  score_covariance <- matrix(0, nrow(dat), nrow(dat))
  for (k in seq_along(effects)) {
    score_covariance <- score_covariance + components[[k]] * tcrossprod(indicator(effects[[k]]))
  }
  mean_covariance <- person_mean %*% score_covariance %*% t(person_mean)
  n_persons <- ncol(persons)
  centering <- diag(n_persons) - 1 / n_persons
  Delta <- mean(diag(mean_covariance)) - components[[1]]
  delta <- sum(diag(centering %*% mean_covariance %*% centering)) / (n_persons - 1) -
    components[[1]]
  return(c(delta = delta, Delta = Delta))
}

d_estimates <- function(d) {
  result <- d$coefficient_results$estimate
  names(result) <- d$coefficient_results$quantity
  return(result)
}

expect_dense_d_study <- function(g, dat, effects) {
  universe_sizes <- rep(Inf, g$n_facets)
  names(universe_sizes) <- c("F1", "F2")[seq_len(g$n_facets)]
  d <- d_study(g, universe_sizes = universe_sizes)
  dense <- dense_errors(dat, effects, unname(g$variance_component))
  expect_equal(d_estimates(d)[c("delta", "Delta")], dense)
}

test_that("ANOVA reproduces GENOVA for p x (r:t) (Brennan Appendix F)", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    method = "anova"
  )

  # Package order p, t, pt, r:t, pr:t; the T terms are in Table 3.4.
  expect_printed(
    g$anova_table$t_term, c(2800.1667, 2755.7000, 2931.5000, 2835.4000, 3204.0000), 5e-5
  )
  expect_printed(
    g$anova_table$sum_of_squares, c(92.66667, 48.20000, 83.13333, 79.70000, 192.80000), 5e-6
  )
  expect_printed(
    g$anova_table$mean_square, c(10.29630, 24.10000, 4.61852, 8.85556, 2.38025), 5e-6
  )
  expect_equal(
    g$ems_coefficients["F1", ], c(P = 0, F1 = 40, PxF1 = 4, "F2:F1" = 10, "F2:F1xP" = 1)
  )
  expect_printed(
    g$variance_component, c(0.4731481, 0.3251543, 0.5595679, 0.6475309, 2.3802469), 5e-8
  )
  expect_printed(
    g$variance_component_se, c(0.3855758, 0.4379875, 0.3766291, 0.3794056, 0.3694860), 5e-8
  )

  # GENOVA's covariance matrix, in its order p, t, r:t, pt, pr:t.
  genova_order <- c("P", "F1", "F2:F1", "PxF1", "F2:F1xP")
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
    brennan_p_x_r_in_t, .build_g_study_design("score", "p", "t", "r", "partially_nested")
  )
  general_vcov <- .anova_covariance_unbalanced(
    dat, package_order, g$anova_table$df, g$ems_coefficients, g$anova_table$estimate
  )
  expect_equal(unname(general_vcov), unname(g$variance_component_vcov))

  # D study #1: 3 tasks and 4 raters per task, both random. GENOVA prints SEs
  # for tau, delta, and Delta.
  d <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
  expect_printed(
    d$coefficient_results$estimate, c(0.47315, 0.38488, 0.54722, 0.55144, 0.46370), 5e-6
  )
  expect_printed(d$coefficient_results$standard_error[1:3], c(0.38558, 0.12171, 0.17935), 5e-6)
})

test_that("ANOVA estimates below zero are set to zero (Brennan Table 3.5)", {
  # Synthetic Data Set No. 3 analyzed as p x (i:o), as if different items were
  # used on the two occasions. Package order p, o, po, i:o, pi:o.
  g <- g_study(
    brennan_p_x_i_x_o, p = "p", f1 = "o", f2 = "i", design = "partially_nested",
    method = "anova"
  )

  expect_printed(g$anova_table$mean_square, c(6.9111, 3.2000, 1.3389, 7.4917, 1.5102), 5e-5)
  expect_printed(g$anova_table$estimate, c(0.6965, -0.1030, -0.0428, 0.5981, 1.5102), 5e-5)
  expect_printed(g$variance_component, c(0.6965, 0, 0, 0.5981, 1.5102), 5e-5)
})

test_that("ANOVA reproduces Brennan Tables 2.3, 2.9, and 3.3", {
  g_p_x_i <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed", method = "anova")
  expect_printed(g_p_x_i$variance_component, c(0.0574, 0.0754, 0.1269), 5e-5)

  g_i_in_p <- g_study(brennan_i_in_p, p = "p", f1 = "i", design = "nested", method = "anova")
  expect_printed(g_i_in_p$variance_component, c(0.6108, 2.0250), 5e-5)

  g_p_x_i_x_o <- g_study(
    brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed",
    method = "anova"
  )
  expect_printed(
    g_p_x_i_x_o$variance_component,
    c(0.5528, 0.4417, 0.0074, 0.5750, 0.1009, 0.1565, 0.9352),
    5e-5
  )
})

test_that("ANOVA with tasks fixed reproduces GENOVA D study #3 (Brennan Appendix F)", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    fixed = "t", method = "anova"
  )
  reported <- g$variance_components$variance_estimate
  names(reported) <- g$variance_components$effect
  expect_printed(
    reported[c("P", "PxF1", "F2:F1", "F2:F1xP")], c(0.65967, 0.55957, 0.64753, 2.38025), 5e-6
  )

  d <- d_study(g, universe_sizes = c(F2 = Inf), d_sizes = c(F2 = 4))
  expect_printed(
    d$coefficient_results$estimate, c(0.65967, 0.19835, 0.25231, 0.76882, 0.72333), 5e-6
  )
  expect_printed(d$coefficient_results$standard_error[1:3], c(0.36716, 0.03079, 0.04193), 5e-6)
})

test_that("ANOVA reproduces Brennan Table 7.2 (unbalanced i:p)", {
  dat <- brennan_i_in_p_unbalanced
  g <- g_study(dat, p = "p", f1 = "i", design = "nested", method = "anova")

  expect_printed(g$anova_table$t_term, c(1069.40, 1120.00), 5e-3)
  expect_printed(g$anova_table$mean_square, c(3.8912, 1.8071), 5e-5)
  expect_printed(g$variance_component, c(0.5098, 1.8071), 5e-5)
  expect_dense_agreement(g, dat, effects_i_in_p)

  # The D study with the G-study pattern (Sec. 7.2.1): n~_i = 3.9130, the
  # harmonic mean of the n_i:p (Eq. 7.19), sigma^2(Delta) = sigma^2(delta) =
  # 1.8071 / 3.9130 = .4618 (Eq. 7.18), and E rho^2 = .525.
  d <- d_study(g, universe_sizes = c(F1 = Inf))
  expect_printed(d$d_mean_results$divisor[[2]], 3.9130, 5e-5)
  estimates <- d_estimates(d)
  expect_printed(estimates[c("delta", "Delta")], c(0.4618, 0.4618), 5e-5)
  expect_printed(estimates[["G"]], 0.525, 5e-4)
  expect_dense_d_study(g, dat, effects_i_in_p)
})

test_that("ANOVA reproduces Brennan Table 7.5 and Eq. 7.13 (unbalanced p x (i:h))", {
  dat <- brennan_p_x_i_in_h_unbalanced
  g <- g_study(dat, p = "p", f1 = "h", f2 = "i", design = "partially_nested", method = "anova")

  # Package order p, h, ph, i:h, pi:h.
  expect_printed(g$anova_table$t_term, c(1390.0, 1424.0, 1564.5, 1440.0, 1610.0), 0.05)
  expect_printed(g$anova_table$mean_square, c(13.4286, 64.0000, 3.3214, 3.2000, 0.8429), 5e-5)
  # Brennan prints .9913 for sigma^2(ph); Eq. 7.14 with his mean squares gives
  # (3.3214 - .8429) / 2.5 = .9914.
  expect_printed(g$variance_component, c(1.2014, 2.9161, 0.9914, 0.2946, 0.8429), 5e-5)

  # Eq. 7.13 with n_p = 8, n_i+ = 8, r_i = 24 / 8 = 3, and t_i = (8 - 3) / 2.
  expected_ems <- rbind(
    c(8, 0, 3, 0, 1),
    c(0, 8 * 2.5, 2.5, 8, 1),
    c(0, 0, 2.5, 0, 1),
    c(0, 0, 0, 8, 1),
    c(0, 0, 0, 0, 1)
  )
  expect_equal(unname(g$ems_coefficients), expected_ems)
  expect_dense_agreement(g, dat, effects_p_x_i_in_h)

  # The D study with the G-study pattern (Sec. 7.2.2): n~_h = 64 / 24 and
  # n_i+ = 8 (Eq. 7.28), sigma^2(Delta) = 1.6075 (Eq. 7.27), sigma^2(delta)
  # = .4771 (Eq. 7.31), and E rho^2 = .716 (Eq. 7.32).
  d <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf))
  expect_equal(d$d_mean_results$divisor, c(1, 64 / 24, 64 / 24, 8, 8))
  estimates <- d_estimates(d)
  expect_printed(estimates[c("delta", "Delta")], c(0.4771, 1.6075), 5e-5)
  expect_printed(estimates[["G"]], 0.716, 5e-4)
  # The standard error of sigma^2(Delta) from the G-study covariance with the
  # weights 1 / n~_h and 1 / n_i+.
  weights <- c(0, 24 / 64, 24 / 64, 1 / 8, 1 / 8)
  expect_equal(
    d$coefficient_results$standard_error[[3]],
    sqrt(as.numeric(t(weights) %*% g$variance_component_vcov %*% weights))
  )
  expect_dense_d_study(g, dat, effects_p_x_i_in_h)
})

test_that("ANOVA reproduces urGENOVA for p x i with missing scores (Brennan Table 7.7)", {
  dat <- brennan_p_x_i_missing
  g <- g_study(dat, p = "p", f1 = "i", design = "crossed", method = "anova")

  # urGENOVA manual, Appendix F (Run 5).
  expect_equal(g$anova_table$df, c(11, 5, 42))
  expect_printed(g$anova_table$t_term, c(27.80000, 24.74318, 37.00000), 5e-6)
  expect_printed(g$anova_table$mean_square, c(0.41787, 0.30796, 0.18239), 5e-6)
  expect_printed(g$variance_component, c(0.04725, 0.01174, 0.18401), 5e-6)
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
  d <- d_study(g, universe_sizes = c(F1 = Inf))
  expect_printed(d$d_mean_results$divisor[2:3], c(4.7682, 4.7682), 5e-5)
  observed <- dat[!is.na(dat$score), ]
  took_item <- unclass(table(observed$p, observed$i))
  n_items <- rowSums(took_item)
  harmonic_mean <- length(n_items) / sum(1 / n_items)
  vc <- g$variance_component
  estimates <- d_estimates(d)
  expect_equal(estimates[["Delta"]], (vc[["F1"]] + vc[["PxF1"]]) / harmonic_mean)
  # sigma^2(delta) from the last two terms of Eq. 7.45: the coefficient of
  # sigma^2(i) is 1 / n~ less the mean, over pairs of persons, of
  # n_pp' / (n_p n_p'), with n_pp' the number of items both took. This gives
  # .0391 and E rho^2 = .547.
  pair_share <- tcrossprod(took_item) / outer(n_items, n_items)
  diag(pair_share) <- 0
  n_persons <- length(n_items)
  coefficient_i <- 1 / harmonic_mean - sum(pair_share) / (n_persons * (n_persons - 1))
  expect_equal(
    estimates[["delta"]], coefficient_i * vc[["F1"]] + vc[["PxF1"]] / harmonic_mean
  )
  expect_dense_d_study(g, dat, effects_p_x_i)
})

test_that("ANOVA reproduces urGENOVA for p x (i:h) with empty cells", {
  dat <- urgenova_p_x_i_in_h()
  g <- g_study(dat, p = "p", f1 = "h", f2 = "i", design = "partially_nested", method = "anova")

  # urGENOVA manual, Appendix G (Run 6). Package order p, h, ph, i:h, pi:h.
  expect_equal(g$anova_table$df, c(9, 2, 18, 9, 75))
  expect_printed(
    g$anova_table$t_term, c(2390.53788, 2360.44089, 2512.00000, 2436.26667, 2787.00000), 5e-6
  )
  expect_printed(
    g$anova_table$mean_square, c(9.11532, 25.97044, 3.86229, 8.42509, 2.65566), 5e-6
  )
  expect_printed(g$variance_component, c(0.43859, 0.45415, 0.32937, 0.60321, 2.67060), 5e-6)
  urgenova_ems <- rbind(
    c(11.398, 0.031, 4.254, 0.051, 1),
    c(0.042, 35.868, 3.621, 9.614, 1),
    c(-0.005, -0.016, 3.567, 0.043, 1),
    c(0.055, 0, 0.055, 9.470, 1),
    c(-0.007, 0, -0.007, -0.016, 1)
  )
  expect_printed(g$ems_coefficients, urgenova_ems, 5e-4)
  expect_dense_agreement(g, dat, effects_p_x_i_in_h)
  expect_dense_d_study(g, dat, effects_p_x_i_in_h)
})

test_that("Expected mean squares give Brennan Table 8.1 from the T terms of Exercise 8.1", {
  # ACT Science Reasoning Form 1: 3458 persons and seven passages with
  # 5, 5, 5, 6, 6, 6, and 7 items.
  items_per_passage <- c(5, 5, 5, 6, 6, 6, 7)
  layout <- expand.grid(
    P = factor(seq_len(3458)),
    F2 = factor(seq_len(sum(items_per_passage)))
  )
  layout$F1 <- factor(rep(seq_along(items_per_passage), items_per_passage))[layout$F2]
  layout$score <- 0
  effect_names <- c("P", "F1", "PxF1", "F2:F1", "F2:F1xP")
  df <- .anova_table(layout, effect_names)$df
  ems_coefficients <- .ems_coefficients(layout, effect_names, df)

  t_p <- 61977.95000
  t_h <- 60279.46591
  t_ih <- 62005.88433
  t_ph <- 69160.59524
  t_pih <- 89416.00000
  t_mean <- 57802.35003
  sum_of_squares <- c(
    t_p - t_mean, t_h - t_mean, t_ph - t_p - t_h + t_mean, t_ih - t_h,
    t_pih - t_ph - t_ih + t_h
  )
  estimates <- solve(ems_coefficients, sum_of_squares / df)
  expect_printed(estimates, c(0.0245, 0.0183, 0.0113, 0.0151, 0.1624), 5e-5)
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
  g <- g_study(dat, p = "p", f1 = "h", f2 = "i", design = "fully_nested", method = "anova")

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
  effects <- list(p = "p", hp = c("h", "p"), ihp = c("i", "h", "p"))
  expect_dense_agreement(g, dat, effects)
  expect_dense_d_study(g, dat, effects)
})

test_that("ANOVA covariance agrees with the dense calculation for p x i x o with missing scores", {
  dat <- brennan_p_x_i_x_o
  dat$score[c(3, 17, 22, 41, 58, 76)] <- NA
  g <- g_study(dat, p = "p", f1 = "i", f2 = "o", design = "fully_crossed", method = "anova")
  effects <- list(
    p = "p", i = "i", o = "o", pi = c("p", "i"), po = c("p", "o"), io = c("i", "o"),
    pio = c("p", "i", "o")
  )
  expect_dense_agreement(g, dat, effects)
  expect_dense_d_study(g, dat, effects)
})

test_that("ANOVA requires positive degrees of freedom for every effect", {
  # Each person takes different items, so p x i has negative df.
  dat <- data.frame(p = rep(1:3, each = 2), i = 1:6, score = c(3, 5, 4, 6, 2, 7))
  expect_error(
    g_study(dat, p = "p", f1 = "i", design = "crossed", method = "anova"),
    "positive degrees of freedom"
  )
})

test_that("A sampling variance estimated below zero gives an NA standard error", {
  scores <- rbind(
    c(4, -1, 7, 5, 3, 7), c(6, 5, 3, 5, 6, 1), c(7, 5, 0, 5, 3, 7), c(5, 3, 6, 7, 5, 4)
  )
  conditions <- data.frame(h = c(1, 1, 1, 2, 3, 3), i = c(1, 2, 3, 1, 1, 2))
  expect_warning(
    g <- g_study(
      long_scores(scores, conditions), p = "p", f1 = "h", f2 = "i",
      design = "partially_nested", method = "anova"
    ),
    "sampling variance of P is below zero"
  )
  expect_true(is.na(g$variance_component_se[["P"]]))
  expect_false(anyNA(g$variance_component_se[-1]))
})

test_that("A D-study sampling variance estimated below zero gives an NA standard error", {
  # Every G-study sampling variance is above zero, but the estimated covariance
  # has a negative eigenvalue.
  scores <- rbind(c(7, 6, 5, 8, 6), c(7, 3, 5, 9, 6), c(6, 7, 7, 9, 7), c(6, 5, 9, 3, 7))
  conditions <- data.frame(h = c(1, 2, 3, 3, 3), i = c(1, 1, 1, 2, 3))
  g <- g_study(
    long_scores(scores, conditions), p = "p", f1 = "h", f2 = "i",
    design = "partially_nested", method = "anova"
  )
  expect_false(anyNA(g$variance_component_se))

  expect_warning(
    d <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 6, F2 = 1)),
    "sampling variance of G is below zero"
  )
  standard_errors <- d$coefficient_results$standard_error
  expect_true(is.na(standard_errors[[4]]))
  expect_false(anyNA(standard_errors[-4]))
})

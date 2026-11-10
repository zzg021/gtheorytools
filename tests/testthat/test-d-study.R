# D studies without d_sizes, which repeat the G study.

test_that("For balanced data a D study without d_sizes uses the G-study sizes", {
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  expect_equal(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(r = Inf, t = 3)),
    d_study_general(
      g, object_of_measurement = "p",
      universe_sizes = c(r = Inf, t = 3), d_sizes = c(r = 4, t = 3)
    )
  )
  expect_error(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(r = 2, t = 3)),
    "r has n = 4 in the G study, more than its universe size, N' = 2"
  )
})

test_that("The G-study pattern gives the same divisors with REML estimates", {
  dat <- brennan_i_in_p_unbalanced
  g <- g_study_general(dat, design = "i:p")
  d <- d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf))

  n_items <- as.vector(table(dat$p))
  harmonic_mean <- length(n_items) / sum(1 / n_items)
  Delta <- g$variance_component[["i:p"]] / harmonic_mean
  expect_equal(d$coefficient_results$estimate[2:3], c(Delta, Delta))
  expect_equal(
    d$coefficient_results$standard_error[[3]],
    g$variance_component_se[["i:p"]] / harmonic_mean
  )
})

test_that("The G-study pattern requires infinite universes and the G-study design", {
  g <- g_study_general(brennan_i_in_p_unbalanced, design = "i:p", method = "anova")
  expect_error(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = 20)),
    "every universe size to be Inf"
  )

  dat <- brennan_p_x_i
  dat$score[5] <- NA
  g <- g_study_general(dat, design = "p x i", method = "anova")
  expect_error(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf), d_design = "i:p"),
    "d_design must be the G-study design"
  )

  # The G-study design written in another order is the same design.
  g <- g_study_general(brennan_p_x_i_in_h_unbalanced, design = "p x (i:h)", method = "anova")
  expect_equal(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf, h = Inf),
            d_design = "(i:h) x p")$coefficient_results,
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf, h = Inf))$coefficient_results
  )
})

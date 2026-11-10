# D studies without d_sizes, which repeat the G study.

test_that("For balanced data a D study without d_sizes uses the G-study sizes", {
  g <- g_study(brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested")
  expect_equal(
    d_study(g, universe_sizes = c(F1 = 3, F2 = Inf)),
    d_study(g, universe_sizes = c(F1 = 3, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
  )
  expect_error(
    d_study(g, universe_sizes = c(F1 = 3, F2 = 2)),
    "r has n = 4 in the G study, more than universe_sizes\\[\"F2\"\\] = 2"
  )
})

test_that("The G-study pattern gives the same divisors with REML estimates", {
  dat <- brennan_i_in_p_unbalanced
  g <- g_study(dat, p = "p", f1 = "i", design = "nested")
  d <- d_study(g, universe_sizes = c(F1 = Inf))

  n_items <- as.vector(table(dat$p))
  harmonic_mean <- length(n_items) / sum(1 / n_items)
  Delta <- g$variance_component[["F1:P"]] / harmonic_mean
  expect_equal(d$coefficient_results$estimate[2:3], c(Delta, Delta))
  expect_equal(
    d$coefficient_results$standard_error[[3]],
    g$variance_component_se[["F1:P"]] / harmonic_mean
  )
})

test_that("The G-study pattern requires infinite universes and the G-study design", {
  g <- g_study(brennan_i_in_p_unbalanced, p = "p", f1 = "i", design = "nested", method = "anova")
  expect_error(d_study(g, universe_sizes = c(F1 = 20)), "every universe size to be Inf")

  dat <- brennan_p_x_i
  dat$score[5] <- NA
  g <- g_study(dat, p = "p", f1 = "i", design = "crossed", method = "anova")
  expect_error(
    d_study(g, universe_sizes = c(F1 = Inf), d_design = "nested"),
    "d_design must be the G-study design"
  )
})

# ANOVA estimates against values printed by Brennan (2001) and GENOVA. Each
# tolerance is half a unit in the last printed decimal.

expect_printed <- function(object, expected, tolerance) {
  expect_lt(max(abs(unname(object) - unname(expected))), tolerance)
}

test_that("ANOVA reproduces GENOVA for p x (r:t) (Brennan Appendix F)", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    method = "anova"
  )

  # Package order p, t, pt, r:t, pr:t.
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

# g_study(), d_study(), and compare_d_studies() compute with the design string
# of each keyword, so they agree with the general functions, with the effects
# named as in version 0.1.0.

# brennan_p_x_r_in_t with different tasks and raters for each person.
r_in_t_in_p <- function() {
  dat <- brennan_p_x_r_in_t
  dat$t <- paste(dat$p, dat$t)
  dat$r <- paste(dat$t, dat$r)
  return(dat)
}

test_that("Each design keyword gives the components of its design string", {
  cases <- list(
    list(brennan_p_x_i, "i", NULL, "crossed", "p x i", c(P = "p", F1 = "i", PxF1 = "pi")),
    list(brennan_i_in_p, "i", NULL, "nested", "i:p", c(P = "p", "F1:P" = "i:p")),
    list(
      brennan_p_x_i_x_o, "i", "o", "fully_crossed", "p x i x o",
      c(P = "p", F1 = "i", F2 = "o", PxF1 = "pi", PxF2 = "po", F1xF2 = "io",
        PxF1xF2 = "pio")
    ),
    list(
      brennan_p_x_r_in_t, "t", "r", "partially_nested", "p x (r:t)",
      c(P = "p", F1 = "t", PxF1 = "pt", "F2:F1" = "r:t", "F2:F1xP" = "pr:t")
    ),
    list(
      r_in_t_in_p(), "t", "r", "fully_nested", "r:t:p",
      c(P = "p", "F1:P" = "t:p", "F2:F1:P" = "r:t:p")
    )
  )

  for (case in cases) {
    names(case) <- c("data", "f1", "f2", "keyword", "string", "effects")
    for (method in c("reml", "anova")) {
      keyword <- g_study(
        case$data, p = "p", f1 = case$f1, f2 = case$f2, design = case$keyword,
        method = method
      )
      general <- g_study_general(case$data, design = case$string, method = method)
      effects <- case$effects
      expect_equal(keyword$variance_components$effect, names(effects))
      expect_equal(
        unname(keyword$variance_component), unname(general$variance_component[effects]),
        tolerance = 1e-6
      )
      expect_equal(
        unname(keyword$variance_component_vcov),
        unname(general$variance_component_vcov[effects, effects]),
        tolerance = 1e-6
      )
    }
  }
})

test_that("Each D-study design keyword gives the D study of its design string", {
  g_pio <- g_study(
    brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed",
    method = "anova"
  )
  g_pio_general <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o", method = "anova")
  g_prt <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    method = "anova"
  )
  g_prt_general <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", method = "anova")
  g_pi <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed", method = "anova")
  g_pi_general <- g_study_general(brennan_p_x_i, design = "p x i", method = "anova")

  cases <- list(
    list(g_pio, g_pio_general, NULL, NULL, "p x i x o",
         c(P = "p", F1 = "i", F2 = "o", PxF1 = "pi", PxF2 = "po", F1xF2 = "io",
           PxF1xF2 = "pio")),
    list(g_pio, g_pio_general, "partially_nested", "F2:F1", "p x (o:i)",
         c(P = "p", F1 = "i", PxF1 = "pi", "F2:F1" = "o:i", "F2:F1xP" = "po:i")),
    list(g_pio, g_pio_general, "partially_nested", "F1:F2", "p x (i:o)",
         c(P = "p", F2 = "o", PxF2 = "po", "F1:F2" = "i:o", "F1:F2xP" = "pi:o")),
    list(g_pio, g_pio_general, "fully_nested", "F2:F1", "o:i:p",
         c(P = "p", "F1:P" = "i:p", "F2:F1:P" = "o:i:p")),
    list(g_pio, g_pio_general, "fully_nested", "F1:F2", "i:o:p",
         c(P = "p", "F2:P" = "o:p", "F1:F2:P" = "i:o:p")),
    list(g_prt, g_prt_general, "fully_nested", NULL, "r:t:p",
         c(P = "p", "F1:P" = "t:p", "F2:F1:P" = "r:t:p")),
    list(g_pi, g_pi_general, "nested", NULL, "i:p", c(P = "p", "F1:P" = "i:p"))
  )

  for (case in cases) {
    names(case) <- c("g", "g_general", "d_design", "d_nesting", "string", "effects")
    n_facets <- case$g$n_facets
    facet_names <- c("F1", "F2")[seq_len(n_facets)]
    columns <- case$g$facet_columns[facet_names]
    d_sizes <- c(4, 2)[seq_len(n_facets)]
    names(d_sizes) <- facet_names
    universe_sizes <- c(F1 = Inf, F2 = Inf)[facet_names]
    keyword <- d_study(
      case$g, universe_sizes = universe_sizes, d_sizes = d_sizes,
      d_design = case$d_design, d_nesting = case$d_nesting
    )
    general_sizes <- d_sizes
    names(general_sizes) <- columns
    general_universe <- universe_sizes
    names(general_universe) <- columns
    general <- d_study_general(
      case$g_general, object_of_measurement = "p", universe_sizes = general_universe,
      d_sizes = general_sizes, d_design = case$string
    )

    effects <- case$effects
    expect_equal(keyword$d_mean_results$effect, names(effects))
    mean_scores <- general$d_mean_results$mean_score_estimate
    names(mean_scores) <- general$d_mean_results$effect
    expect_equal(keyword$d_mean_results$mean_score_estimate, unname(mean_scores[effects]))
    expect_equal(keyword$coefficient_results, general$coefficient_results)
    expect_equal(keyword$d_sizes, d_sizes)
  }

  # Brennan's (2001, p. 63) rule for p x (i:o) from p x i x o.
  recoded <- d_study(
    g_pio, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 4, F2 = 2),
    d_design = "partially_nested", d_nesting = "F1:F2"
  )
  expect_equal(
    recoded$component_mapping$source_effects,
    c("P", "F2", "PxF2", "F1 + F1xF2", "PxF1 + PxF1xF2")
  )
})

test_that("A facet fixed in g_study() and finite universes reach the general functions", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    fixed = "t", method = "anova"
  )
  g_general <- g_study_general(
    brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t", method = "anova"
  )
  effects <- c(P = "p", F1 = "t", PxF1 = "pt", "F2:F1" = "r:t", "F2:F1xP" = "pr:t")
  reported <- g_general$variance_components
  expect_equal(
    g$variance_components$variance_estimate,
    reported$variance_estimate[match(effects, reported$effect)]
  )
  expect_equal(g$fixed_facet, "F1")

  keyword <- d_study(g, universe_sizes = c(F2 = 20), d_sizes = c(F2 = 4))
  general <- d_study_general(
    g_general, object_of_measurement = "p", universe_sizes = c(r = 20), d_sizes = c(r = 4)
  )
  expect_equal(keyword$coefficient_results, general$coefficient_results, tolerance = 1e-6)
  expect_equal(keyword$universe_sizes, c(F1 = 3, F2 = 20))
  expect_equal(keyword$facet_status, c(F1 = "fixed", F2 = "finite"))
})

test_that("d_study() repeats the G-study pattern of unbalanced data", {
  g <- g_study(
    brennan_p_x_i_in_h_unbalanced, p = "p", f1 = "h", f2 = "i",
    design = "partially_nested", method = "anova"
  )
  g_general <- g_study_general(
    brennan_p_x_i_in_h_unbalanced, design = "p x (i:h)", method = "anova"
  )
  keyword <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf))
  general <- d_study_general(
    g_general, object_of_measurement = "p", universe_sizes = c(i = Inf, h = Inf)
  )
  effects <- c(P = "p", F1 = "h", PxF1 = "ph", "F2:F1" = "i:h", "F2:F1xP" = "pi:h")
  divisors <- general$d_mean_results$divisor
  names(divisors) <- general$d_mean_results$effect
  expect_equal(keyword$d_mean_results$divisor, unname(divisors[effects]))
  expect_equal(keyword$coefficient_results, general$coefficient_results)
  expect_null(keyword$d_sizes)
})

test_that("compare_d_studies() varies F1 fastest and works with find_d_designs()", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
    method = "anova"
  )
  g_general <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", method = "anova")
  keyword <- compare_d_studies(
    g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = list(F1 = 1:3, F2 = c(2, 4))
  )
  general <- compare_d_studies_general(
    g_general, object_of_measurement = "p", universe_sizes = c(r = Inf, t = Inf),
    d_sizes = data.frame(t = rep(1:3, 2), r = rep(c(2, 4), each = 3))
  )
  expect_equal(keyword$scenarios$F1, rep(1:3, 2))
  expect_equal(keyword$scenarios$F2, rep(c(2, 4), each = 3))
  expect_equal(keyword$scenarios[c("G", "Phi", "se_G")], general$scenarios[c("G", "Phi", "se_G")])
  expect_equal(
    find_d_designs(keyword, "G", 0.6)$designs$scenario,
    find_d_designs(general, "G", 0.6)$designs$scenario
  )
})

test_that("The keyword functions check their own arguments", {
  expect_error(
    g_study(brennan_p_x_i, p = "p", f1 = "i", design = "p x i"),
    "g_study_general\\(\\) takes a design string"
  )
  expect_error(
    g_study(brennan_p_x_i, p = "p", f1 = "p", design = "crossed"),
    "p and f1 must name different columns"
  )
  expect_error(
    g_study(
      brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
      fixed = "r"
    ),
    "only f1 can be fixed"
  )

  g <- g_study(brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed")
  expect_error(
    d_study(
      g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 4, F2 = 2),
      d_design = "partially_nested"
    ),
    "d_nesting must be specified"
  )
  expect_error(
    d_study(g, universe_sizes = c(i = Inf, o = Inf), d_sizes = c(i = 4, o = 2)),
    "containing exactly: F1, F2"
  )
  expect_error(
    compare_d_studies(g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = list(F1 = 4)),
    "d_sizes must name each facet \\(F1, F2\\)"
  )
  expect_error(
    d_study_general(g, object_of_measurement = "p", universe_sizes = c(i = Inf, o = Inf)),
    "result of g_study\\(\\) goes to d_study\\(\\)"
  )
  g_general <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
  expect_error(
    d_study(g_general, universe_sizes = c(F1 = Inf, F2 = Inf)),
    "result of g_study_general\\(\\) goes to d_study_general\\(\\)"
  )
})

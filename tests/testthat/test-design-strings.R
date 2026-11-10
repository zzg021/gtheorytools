# Design strings: the nesting relation and Brennan's effect names of each,
# the errors for strings that are not designs, and the facet symbols as data
# columns.

test_that("A design string gives its nesting relation and Brennan's effect names", {
  none <- character(0)
  cases <- list(
    list(
      design = "p x (r:t) x o",
      parents = list(p = none, r = "t", t = none, o = none),
      effects = c("p", "t", "o", "r:t", "pt", "po", "to", "pr:t", "ro:t", "pto", "pro:t")
    ),
    list(
      design = "p x (i:h:g)",
      parents = list(p = none, i = c("h", "g"), h = "g", g = none),
      effects = c("p", "g", "h:g", "i:h:g", "pg", "ph:g", "pi:h:g")
    ),
    list(
      design = "p x (i:(r x h))",
      parents = list(p = none, i = c("r", "h"), r = none, h = none),
      effects = c("p", "r", "h", "i:rh", "pr", "ph", "rh", "prh", "pi:rh")
    ),
    list(
      design = "(i:p) x o",
      parents = list(i = "p", p = none, o = none),
      effects = c("p", "o", "i:p", "po", "io:p")
    ),
    list(
      design = "p x (a:b) x (c:d:b)",
      parents = list(p = none, a = "b", b = none, c = c("b", "d"), d = "b"),
      effects = c("p", "b", "a:b", "d:b", "c:d:b", "pb", "pa:b", "pd:b", "ad:b",
                  "pc:d:b", "ac:d:b", "pad:b", "pac:d:b")
    )
  )
  for (case in cases) {
    parsed <- .parse_design(case$design, "design")
    expect_identical(parsed$parents, case$parents)
    expect_identical(names(.design_effects(parsed$parents)), case$effects)
  }
  expect_identical(.parse_design("pxi", "design")$design, "p x i")
})

test_that("Strings that are not designs are rejected", {
  for (bad_design in c("p x (r:t", "p x", "p i", "p x (x:t)")) {
    expect_error(
      g_study_general(brennan_p_x_r_in_t, design = bad_design),
      "Could not read design"
    )
  }
  expect_error(
    g_study_general(brennan_p_x_r_in_t, design = "p x (r:t) x (t:r)"),
    "nests r within itself"
  )
  expect_error(
    g_study_general(brennan_p_x_i_x_o, design = "fully_crossed"),
    "Here design is a design string, such as \"p x i x o\"",
    fixed = TRUE
  )
  expect_error(g_study_general(brennan_p_x_i, design = "p"), "at least two facets")
})

test_that("Facet symbols name the data columns in their case", {
  upper <- brennan_p_x_r_in_t
  names(upper) <- c("P", "T", "R", "score")
  for (method in c("reml", "anova")) {
    g_upper <- g_study_general(upper, design = "P X (R:T)", method = method)
    g_lower <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", method = method)
    expect_equal(g_upper$g_design, "P x (R:T)")
    expect_equal(g_upper$variance_components$effect, c("P", "T", "R:T", "PT", "PR:T"))
    expect_equal(unname(g_upper$variance_component), unname(g_lower$variance_component))

    d_upper <- d_study_general(g_upper, "P", c(R = Inf, T = 3), c(R = 4, T = 3))
    d_lower <- d_study_general(g_lower, "p", c(r = Inf, t = 3), c(r = 4, t = 3))
    expect_equal(d_upper$coefficient_results, d_lower$coefficient_results)
  }

  expect_error(
    g_study_general(brennan_p_x_r_in_t, design = "P x (R:T)"),
    paste0(
      "missing from data: P, R, T. .*facets = c\\(P = \"P_column\"\\). ",
      "Names are case-sensitive: P is not p, R is not r, T is not t."
    )
  )
  g_lower <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  expect_error(
    d_study_general(g_lower, "p", c(r = Inf, t = Inf), d_design = "p x (R:T)"),
    paste0(
      "d_design must have the facets of the G-study design: p, r, t. ",
      "Names are case-sensitive: R is not r, T is not t."
    ),
    fixed = TRUE
  )
})

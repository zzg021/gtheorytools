# g_study(), d_study(), and compare_d_studies() name the one- and two-facet
# designs by keyword, with persons as the object of measurement. They compute
# with design strings in which p, a, and b stand for the columns p, f1, and
# f2, and report the components as P, F1, and F2 joined by x for crossing and
# by : for nesting.

.keyword_symbols <- c(P = "p", F1 = "a", F2 = "b")

.keyword_design_labels <- c(
  crossed = "One-facet crossed G-study: P x F1",
  nested = "One-facet nested G-study: F1:P",
  fully_crossed = "Two-facet fully crossed G-study: P x F1 x F2",
  partially_nested = "Two-facet partially nested G-study: P x F1 x (F2:F1)",
  fully_nested = "Two-facet fully nested G-study: F2:F1:P"
)

# The design string of each keyword and the keyword name of each of its
# effects, named by the effect of the design string and in reported order.
# The entries ending in F1:F2 are the nested D-study designs of a fully crossed
# G study with F1 nested within F2.
.keyword_designs <- list(
  crossed = list(string = "p x a", effects = c(p = "P", a = "F1", pa = "PxF1")),
  nested = list(string = "a:p", effects = c(p = "P", "a:p" = "F1:P")),
  fully_crossed = list(
    string = "p x a x b",
    effects = c(
      p = "P", a = "F1", b = "F2", pa = "PxF1", pb = "PxF2", ab = "F1xF2",
      pab = "PxF1xF2"
    )
  ),
  partially_nested = list(
    string = "p x (b:a)",
    effects = c(p = "P", a = "F1", pa = "PxF1", "b:a" = "F2:F1", "pb:a" = "F2:F1xP")
  ),
  "partially_nested F1:F2" = list(
    string = "p x (a:b)",
    effects = c(p = "P", b = "F2", pb = "PxF2", "a:b" = "F1:F2", "pa:b" = "F1:F2xP")
  ),
  fully_nested = list(
    string = "b:a:p",
    effects = c(p = "P", "a:p" = "F1:P", "b:a:p" = "F2:F1:P")
  ),
  "fully_nested F1:F2" = list(
    string = "a:b:p",
    effects = c(p = "P", "b:p" = "F2:P", "a:b:p" = "F1:F2:P")
  )
)

.keyword_design <- function(keyword, nesting = NULL) {
  if (identical(nesting, "F1:F2")) {
    keyword <- paste(keyword, nesting)
  }
  return(.keyword_designs[[keyword]])
}

# An effect name written with the data's column names: F2:F1xP as "r:t x p".
.keyword_effect_label <- function(effect_name, facet_columns) {
  tokens <- regmatches(effect_name, gregexpr("F1|F2|P|[:x]", effect_name))[[1]]
  parts <- tokens
  is_facet <- tokens %in% names(facet_columns)
  parts[is_facet] <- facet_columns[tokens[is_facet]]
  parts[tokens == "x"] <- " x "
  return(paste0(parts, collapse = ""))
}

# Positions that put effects of the design string in the order of
# `effect_map`, the keyword names of `.keyword_designs`.
.reported_order <- function(effect_names, effect_map) {
  return(order(match(effect_names, names(effect_map))))
}

.rename_effect_vector <- function(values, effect_map) {
  values <- values[.reported_order(names(values), effect_map)]
  names(values) <- unname(effect_map[names(values)])
  return(values)
}

# A matrix with effects as row and column names, or a one-row matrix with
# effects as column names.
.rename_effect_matrix <- function(values, effect_map) {
  if (is.null(values)) {
    return(NULL)
  }
  columns <- .reported_order(colnames(values), effect_map)
  if (is.null(rownames(values))) {
    values <- values[, columns, drop = FALSE]
    colnames(values) <- unname(effect_map[colnames(values)])
    return(values)
  }
  values <- values[.reported_order(rownames(values), effect_map), columns, drop = FALSE]
  dimnames(values) <- list(
    unname(effect_map[rownames(values)]), unname(effect_map[colnames(values)])
  )
  return(values)
}

# A table with an effect column: keyword names and labels, in reported order.
.rename_effect_table <- function(table, effect_map, facet_columns) {
  table <- table[.reported_order(table$effect, effect_map), , drop = FALSE]
  table$effect <- unname(effect_map[table$effect])
  table$effect_label <- vapply(
    table$effect,
    .keyword_effect_label,
    character(1),
    facet_columns = facet_columns,
    USE.NAMES = FALSE
  )
  rownames(table) <- NULL
  return(table)
}

# Sizes named F1 and F2 renamed to the facet symbols a and b.
.symbol_sizes <- function(sizes) {
  if (is.null(sizes)) {
    return(NULL)
  }
  names(sizes) <- unname(.keyword_symbols[names(sizes)])
  return(sizes)
}

# Sizes named by facet symbol renamed to F1 and F2, in that order.
.keyword_sizes <- function(sizes, facet_names) {
  if (is.null(sizes)) {
    return(NULL)
  }
  result <- sizes[unname(.keyword_symbols[facet_names])]
  names(result) <- facet_names
  return(result)
}

.check_keyword_design <- function(score_column, p, f1, f2, design) {
  .check_single_string(score_column, "score_column")
  .check_single_string(p, "p")
  .check_single_string(f1, "f1")
  .check_single_string(design, "design")

  if (is.null(f2)) {
    if (!(design %in% c("crossed", "nested"))) {
      stop(
        "For one-facet G-studies with f2 = NULL, design must be one of: crossed, ",
        "nested. g_study_general() takes a design string such as \"p x i\".",
        call. = FALSE
      )
    }
    if (p == f1) {
      stop("p and f1 must name different columns.", call. = FALSE)
    }
    return(invisible(NULL))
  }

  .check_single_string(f2, "f2")
  if (!(design %in% c("fully_crossed", "partially_nested", "fully_nested"))) {
    stop(
      "For two-facet G-studies, design must be one of: fully_crossed, ",
      "partially_nested, fully_nested. g_study_general() takes a design string ",
      "such as \"p x i x o\".",
      call. = FALSE
    )
  }
  if (anyDuplicated(c(p, f1, f2)) > 0) {
    stop("p, f1, and f2 must name different columns.", call. = FALSE)
  }
  return(invisible(NULL))
}

# Resolve `fixed` to F1 or F2 and check that the design allows it.
.resolve_keyword_fixed <- function(fixed, keyword, facet_columns) {
  if (is.null(fixed)) {
    return(NULL)
  }

  .check_single_string(fixed, "fixed")
  facet <- names(facet_columns)[facet_columns == fixed]
  if (length(facet) != 1 || facet == "P") {
    stop("fixed must name the f1 or f2 column.", call. = FALSE)
  }

  if (keyword == "fully_crossed") {
    return(facet)
  }
  if (keyword == "partially_nested") {
    if (facet == "F2") {
      stop(
        "In the partially nested design only f1 can be fixed: the levels of ",
        fixed, " differ across the levels of ", facet_columns[["F1"]], ".",
        call. = FALSE
      )
    }
    return(facet)
  }
  if (keyword == "crossed") {
    stop(
      "fixed is not available for the crossed design: fixing its only facet ",
      "leaves no random facet.",
      call. = FALSE
    )
  }

  stop(
    "fixed is not available for the ", keyword, " design: a facet nested ",
    "within persons has no common set of conditions.",
    call. = FALSE
  )
}

#' Generalizability (G) Study
#'
#' Estimates the variance components of a one- or two-facet G study, with
#' persons as the object of measurement, by REML with [lme4::lmer()] or from
#' the ANOVA expected mean squares.
#'
#' @inheritParams g_study_general
#' @param p Name of the person (object of measurement) column.
#' @param f1,f2 Names of the facet columns. `f2 = NULL` gives a one-facet
#'   design. In the nested designs `f2` is nested within `f1`.
#' @param design One of `"crossed"` or `"nested"` for one facet, and
#'   `"fully_crossed"`, `"partially_nested"`, or `"fully_nested"` for two
#'   facets (see Details).
#' @param fixed Optional name of the `f1` or `f2` column of a facet fixed in
#'   the universe of admissible observations. `fully_crossed` allows either
#'   facet and `partially_nested` only `f1`. Requires complete, balanced data.
#'
#' @details
#' P, F1, and F2 correspond to the columns named by `p`, `f1`, and `f2`.
#' In the design notation, `x` denotes crossing and `F2:F1` means F2 is
#' nested within F1. The last component listed is the residual.
#'
#' | `design` | Design | Components | Design string |
#' |---|---|---|---|
#' | `crossed` | P x F1 | P, F1, PxF1 | `"p x i"` |
#' | `nested` | F1:P | P, F1:P | `"i:p"` |
#' | `fully_crossed` | P x F1 x F2 | P, F1, F2, PxF1, PxF2, F1xF2, PxF1xF2 | `"p x i x o"` |
#' | `partially_nested` | P x (F2:F1) | P, F1, PxF1, F2:F1, F2:F1xP | `"p x (r:t)"` |
#' | `fully_nested` | F2:F1:P | P, F1:P, F2:F1:P | `"r:t:p"` |
#'
#' The last column gives the same design for [g_study_general()], with i, o,
#' r, and t standing for the facet columns.
#'
#' A facet declared with `fixed` keeps its observed conditions in every
#' D study. To fix a facet only in the D study, leave it random here and
#' set its sample size equal to its universe size in [d_study()].
#'
#' @inheritSection g_study_general ANOVA estimates
#'
#' @returns An object of class `g_study`, a list with the following elements:
#' \describe{
#'   \item{\code{variance_components}}{A data frame of reported variance
#'     components and standard errors, with effect labels built from the
#'     data's column names. With `fixed`, these are components for the universe
#'     with that facet fixed; the fixed facet's main-effect row is `NA`.}
#'   \item{\code{variance_component}, \code{variance_component_se},
#'     \code{variance_component_vcov}}{The fitted variance components,
#'     their standard errors, and the covariance matrix used by [d_study()].}
#'   \item{\code{method}}{The estimation method.}
#'   \item{\code{anova_table}}{With `method = "anova"`, a data frame of
#'     degrees of freedom, T terms (Brennan, 2001, Eq. 7.1), sums of squares,
#'     mean squares, and the random effects estimates before those below zero
#'     are set to zero.}
#'   \item{\code{ems_coefficients}}{With `method = "anova"`, the matrix whose
#'     row alpha and column beta give the coefficient of sigma^2(beta) in the
#'     expected mean square of alpha.}
#'   \item{\code{fit}}{The fitted `lmer` model.}
#'   \item{\code{fit_diagnostics}}{Singularity and convergence information,
#'     with messages and warnings from `lmer`.}
#'   \item{Design information}{\code{g_design}, \code{design_label},
#'     \code{formula}, \code{n_facets}, \code{facet_columns},
#'     \code{fixed_facet}, and \code{sample_sizes}, which is `NA` for a
#'     nested facet with unequal numbers of conditions.}
#'   \item{\code{observed_design}}{The P, F1, and F2 columns of the observed
#'     scores.}
#'   \item{\code{balanced}}{Whether every cell holds one score and each nested
#'     facet has the same number of conditions in every level of its parent.}
#'   \item{\code{general}}{The same G study as an object of class
#'     `g_study_general`, with p, a, and b for P, F1, and F2, from which
#'     [d_study()] computes.}
#' }
#' Elements for the other method are `NULL`.
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [d_study()]; [g_study_general()] for designs written as design
#'   strings.
#'
#' @examples
#' # p x i (Brennan, 2001, Table 2.3).
#' g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed")
#'
#' # i:p, a different set of items for each person (Table 2.9).
#' g_study(brennan_i_in_p, p = "p", f1 = "i", design = "nested")
#'
#' # p x i x o (Table 3.3).
#' g_study(brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed")
#'
#' # p x (r:t), raters nested within tasks (Table 3.4).
#' g_study(brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested")
#'
#' # The same data with tasks fixed in the universe of admissible observations
#' # (Appendix F).
#' g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
#'   fixed = "t"
#' )
#'
#' # ANOVA estimates, with GENOVA's standard errors (Appendix F).
#' g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
#'   method = "anova"
#' )
#'
#' # p x (i:h) with 2, 4, and 2 items in the levels of h (Table 7.5), and
#' # the analogous ANOVA estimates (Sec. 7.1.3).
#' g_study(
#'   brennan_p_x_i_in_h_unbalanced, p = "p", f1 = "h", f2 = "i",
#'   design = "partially_nested", method = "anova"
#' )
#'
#' @export
g_study <- function(data, score_column = "score", p, f1, f2 = NULL, design,
                    fixed = NULL, method = c("reml", "anova")) {
  method <- match.arg(method)
  .check_keyword_design(score_column, p, f1, f2, design)
  facet_columns <- c(P = p, F1 = f1, F2 = f2)
  fixed_facet <- .resolve_keyword_fixed(fixed, design, facet_columns)
  keyword_design <- .keyword_design(design)
  effect_map <- keyword_design$effects

  # The design string writes the facets as p, a, and b, and the lmer model
  # names them P, F1, and F2.
  facets <- facet_columns
  names(facets) <- .keyword_symbols[names(facet_columns)]
  model_names <- names(facet_columns)
  names(model_names) <- names(facets)
  fixed_symbol <- if (is.null(fixed_facet)) NULL else .keyword_symbols[[fixed_facet]]
  general <- .fit_g_study(
    data, keyword_design$string, facets, score_column, fixed_symbol, method, model_names
  )

  anova_table <- NULL
  ems_coefficients <- NULL
  if (method == "anova") {
    anova_table <- .rename_effect_table(general$anova_table, effect_map, facet_columns)
    ems_coefficients <- .rename_effect_matrix(general$ems_coefficients, effect_map)
  }
  observed_design <- general$observed_design[names(facets)]
  names(observed_design) <- names(facet_columns)
  facet_names <- names(facet_columns)[-1]

  result <- list(
    fit = general$fit,
    fit_diagnostics = general$fit_diagnostics,
    variance_component = .rename_effect_vector(general$variance_component, effect_map),
    variance_component_se = .rename_effect_vector(general$variance_component_se, effect_map),
    variance_component_vcov = .rename_effect_matrix(
      general$variance_component_vcov, effect_map
    ),
    method = method,
    anova_table = anova_table,
    ems_coefficients = ems_coefficients,
    design_label = .keyword_design_labels[[design]],
    formula = general$formula,
    g_design = design,
    n_facets = length(facet_names),
    facet_columns = facet_columns,
    fixed_facet = fixed_facet,
    sample_sizes = .keyword_sizes(general$sample_sizes, facet_names),
    observed_design = observed_design,
    balanced = general$balanced,
    variance_components = .rename_effect_table(
      general$variance_components, effect_map, facet_columns
    ),
    general = general
  )
  class(result) <- "g_study"
  return(result)
}

#' @export
print.g_study <- function(x, ..., digits = 7) {
  cat("\n")
  cat("========================================================\n")
  cat(x$design_label, "\n", sep = "")
  cat("========================================================\n\n")
  general <- x$general
  symbols <- unname(.keyword_symbols[names(x$facet_columns)])
  sample_size_lines <- .sample_size_lines(
    general$observed_design, general$parents[symbols], general$facet_columns
  )
  # With missing scores, the number of scores per person.
  if (length(sample_size_lines) > 1) {
    scores_per_person <- range(table(x$observed_design$P))
    per_person <- paste(unique(scores_per_person), collapse = " to ")
    sample_size_lines[2] <- paste0(sample_size_lines[2], ", ", per_person, " per person")
  }
  cat(sample_size_lines, sep = "\n")
  cat("\n")
  .print_estimation(x, digits)
  cat("\n")
  if (is.null(x$fixed_facet)) {
    cat("G-study variance components:\n")
  } else {
    fixed <- x$facet_columns[[x$fixed_facet]]
    cat(
      "G-study variance components with ", fixed, " fixed (N = n = ",
      x$sample_sizes[[x$fixed_facet]], "):\n",
      sep = ""
    )
  }
  print(x$variance_components, digits = digits, row.names = FALSE)
  if (!is.null(x$fixed_facet)) {
    cat(
      "The ", fixed, " main effect is a fixed effect, so its row is NA.\n",
      sep = ""
    )
  }
  .print_estimation_notes(x, x$fixed_facet)
  return(invisible(x))
}

# The D-study design of a keyword G study: the keyword, its nesting order, and
# its design string and effects.
.resolve_keyword_d_design <- function(g_design, d_design, d_nesting, n_facets) {
  if (n_facets == 1) {
    allowed_designs <- c("crossed", "nested")
  } else {
    allowed_designs <- c("fully_crossed", "partially_nested", "fully_nested")
  }

  if (is.null(d_design)) {
    d_design <- g_design
  } else if (length(d_design) != 1 || !(d_design %in% allowed_designs)) {
    stop("d_design must be one of: ", paste(allowed_designs, collapse = ", "), ".", call. = FALSE)
  }

  nesting_level <- c(
    crossed = 1,
    nested = 2,
    fully_crossed = 1,
    partially_nested = 2,
    fully_nested = 3
  )

  if (nesting_level[[d_design]] < nesting_level[[g_design]]) {
    stop(
      "d_design cannot be less nested than g_study_result$g_design. ",
      "Variance components can be summed into a more nested D-study design, ",
      "but nested components cannot be divided into a less nested design.",
      call. = FALSE
    )
  }

  target_is_nested <- d_design %in% c("partially_nested", "fully_nested")

  if (!is.null(d_nesting)) {
    if (!target_is_nested) {
      stop(
        "d_nesting is only used when d_design is partially_nested or fully_nested.",
        call. = FALSE
      )
    }

    if (length(d_nesting) != 1 || !(d_nesting %in% c("F2:F1", "F1:F2"))) {
      stop("d_nesting must be one of: F2:F1, F1:F2.", call. = FALSE)
    }
  }

  if (g_design == "fully_crossed" && target_is_nested && is.null(d_nesting)) {
    stop(
      "d_nesting must be specified when g_study_result$g_design is fully_crossed ",
      "and d_design is partially_nested or fully_nested. Use d_nesting = \"F2:F1\" ",
      "or d_nesting = \"F1:F2\".",
      call. = FALSE
    )
  }

  # A nested two-facet G study already sets the nesting order as F2:F1.
  if (g_design %in% c("partially_nested", "fully_nested")) {
    if (!is.null(d_nesting) && d_nesting != "F2:F1") {
      stop(
        "d_nesting conflicts with the nesting order in g_study_result$g_design. ",
        "This G-study object uses d_nesting = \"F2:F1\".",
        call. = FALSE
      )
    }

    d_nesting <- "F2:F1"
  }

  keyword_design <- .keyword_design(d_design, d_nesting)
  return(list(
    d_design = d_design,
    d_nesting = d_nesting,
    string = keyword_design$string,
    effects = keyword_design$effects
  ))
}

# The checks of d_study() and compare_d_studies() that do not depend on
# d_sizes, with the sizes named F1 and F2.
.prepare_keyword_d_study <- function(g_study_result, universe_sizes, d_design, d_nesting) {
  .validate_g_study_result(g_study_result, "g_study")
  facet_names <- names(g_study_result$facet_columns)[-1]
  fixed_sizes <- g_study_result$sample_sizes[g_study_result$fixed_facet]
  universe_sizes <- .fill_fixed_facet_sizes(
    universe_sizes, "universe_sizes", fixed_sizes, g_study_result$facet_columns
  )
  universe_values <- .check_size_vector(
    universe_sizes, "universe_sizes", facet_names, allow_inf = TRUE
  )[facet_names]
  design <- .resolve_keyword_d_design(
    g_study_result$g_design, d_design, d_nesting, g_study_result$n_facets
  )
  return(list(
    facet_names = facet_names,
    fixed_sizes = fixed_sizes,
    universe_sizes = universe_values,
    design = design
  ))
}

# The result of d_study_general() with the keyword names of the D-study and
# G-study effects.
.keyword_d_study_result <- function(general, prepared, g_study_result) {
  facet_columns <- g_study_result$facet_columns
  facet_names <- prepared$facet_names
  design <- prepared$design
  d_map <- design$effects
  g_map <- .keyword_design(g_study_result$g_design)$effects

  # The G-study components summed into each D-study component.
  general_mapping <- general$component_mapping
  general_mapping <- general_mapping[.reported_order(general_mapping$effect, d_map), ]
  source_names <- lapply(
    strsplit(general_mapping$source_effects, " + ", fixed = TRUE),
    function(sources) unname(g_map[sources[.reported_order(sources, g_map)]])
  )
  names(source_names) <- unname(d_map[general_mapping$effect])
  source_labels <- vapply(
    source_names,
    function(sources) {
      labels <- vapply(sources, .keyword_effect_label, character(1), facet_columns = facet_columns)
      return(paste(labels, collapse = " + "))
    },
    character(1)
  )
  component_mapping <- data.frame(
    effect = names(source_names),
    source_effects = vapply(source_names, paste, character(1), collapse = " + "),
    row.names = NULL
  )

  single_obs_results <- .rename_effect_table(general$single_obs_results, d_map, facet_columns)
  single_obs_results$source_effects <- unname(source_labels[single_obs_results$effect])
  d_mean_results <- .rename_effect_table(general$d_mean_results, d_map, facet_columns)
  d_mean_results$source_effects <- unname(source_labels[d_mean_results$effect])

  facet_status <- .keyword_sizes(general$facet_status, facet_names)
  universe_sizes <- .keyword_sizes(general$universe_sizes, facet_names)
  d_sizes <- .keyword_sizes(general$d_sizes, facet_names)

  result <- list(
    g_design = g_study_result$g_design,
    d_design = design$d_design,
    d_nesting = design$d_nesting,
    n_facets = g_study_result$n_facets,
    component_mapping = component_mapping,
    d_variance_component = .rename_effect_vector(general$d_variance_component, d_map),
    d_variance_component_vcov = .rename_effect_matrix(
      general$d_variance_component_vcov, d_map
    ),
    single_obs = .rename_effect_matrix(general$single_obs, d_map),
    single_obs_vcov = .rename_effect_matrix(general$single_obs_vcov, d_map),
    single_obs_results = single_obs_results,
    d_mean = .rename_effect_matrix(general$d_mean, d_map),
    d_mean_vcov = .rename_effect_matrix(general$d_mean_vcov, d_map),
    d_mean_results = d_mean_results,
    coefficient_results = general$coefficient_results,
    universe_sizes = universe_sizes,
    d_sizes = d_sizes,
    facet_status = facet_status,
    facet_labels = .d_facet_labels(
      facet_status, universe_sizes, d_sizes, g_study_result$fixed_facet, facet_columns
    ),
    se_available = general$se_available
  )
  class(result) <- "d_study"
  return(result)
}

#' Decision (D) Study from a G Study
#'
#' Computes D-study variance components, error variances, and the
#' generalizability (G) and dependability (Phi) coefficients for one set of
#' facet sample sizes. Standard errors use the G-study covariance matrix
#' and the delta method.
#'
#' @inheritParams d_study_general
#' @param g_study_result The object returned by [g_study()].
#' @param universe_sizes Universe sizes N', named `F1`, or `F1` and `F2`:
#'   whole numbers or `Inf` for an infinite universe. A facet fixed in the
#'   G study may be omitted.
#' @param d_design The D-study design, one of the `design` values of
#'   [g_study()]. Defaults to the G-study design and may be more nested.
#' @param d_nesting `"F2:F1"` or `"F1:F2"`, required when a fully crossed
#'   G study is recoded to a nested design. It is `"F2:F1"` for a nested
#'   G study.
#'
#' @inheritSection d_study_general Coefficients
#' @inheritSection d_study_general Unbalanced data
#'
#' @details
#' A facet is random when N' is infinite, finite when `n' < N' < Inf`, and
#' fixed in the D study when `n' = N'`.
#'
#' @returns An object of class `d_study`, a list with the following elements:
#' \describe{
#'   \item{\code{coefficient_results}}{A data frame with estimates and
#'     standard errors for `tau`, `delta`, `Delta`, `G`, and `Phi`.}
#'   \item{\code{single_obs_results}}{A data frame of single-observation
#'     variance components and standard errors.}
#'   \item{\code{d_mean_results}}{A data frame of mean-score variance
#'     components and standard errors, with `divisor` and `delta_weight` for
#'     the G-study pattern of unbalanced data.}
#'   \item{\code{d_variance_component}}{Variance components after recoding
#'     to the D-study design.}
#'   \item{\code{single_obs}, \code{d_mean}}{The single-observation and
#'     mean-score component estimates as one-row matrices.}
#'   \item{\code{d_variance_component_vcov}, \code{single_obs_vcov},
#'     \code{d_mean_vcov}}{Covariance matrices for the corresponding component
#'     estimates. `single_obs_vcov` and `d_mean_vcov` are `NULL` when
#'     `se_available` is `FALSE`.}
#'   \item{Design information}{\code{g_design}, \code{d_design},
#'     \code{d_nesting}, \code{n_facets}, \code{component_mapping},
#'     \code{universe_sizes}, \code{d_sizes}, \code{facet_status}, and
#'     \code{facet_labels}. `d_sizes` is `NULL` for the G-study pattern of
#'     unbalanced data.}
#'   \item{\code{se_available}}{Whether the G-study covariance matrix can
#'     be used to calculate standard errors.}
#' }
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [compare_d_studies()] for several sets of sizes;
#'   [d_study_general()] for the result of [g_study_general()].
#'
#' @examples
#' # p x i with 12 items (Brennan, 2001, Table 2.4).
#' g_pi <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed")
#' d_study(g_pi, universe_sizes = c(F1 = Inf), d_sizes = c(F1 = 12))
#'
#' # p x i with missing scores (Table 7.7), and the D study that repeats the
#' # G-study pattern (Sec. 7.2.4).
#' g_missing <- g_study(
#'   brennan_p_x_i_missing, p = "p", f1 = "i", design = "crossed",
#'   method = "anova"
#' )
#' d_study(g_missing, universe_sizes = c(F1 = Inf))
#'
#' # p x i x o with 4 items and 2 occasions (Table 4.4).
#' g_pio <- g_study(
#'   brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed"
#' )
#' d_study(g_pio, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 4, F2 = 2))
#'
#' # The same sample sizes with different items on each occasion, p x (i:o)
#' # (Table 4.5).
#' d_study(
#'   g_pio,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = c(F1 = 4, F2 = 2),
#'   d_design = "partially_nested",
#'   d_nesting = "F1:F2"
#' )
#'
#' # p x (r:t) with the 3 tasks fixed and 4 raters per task (Table 4.8).
#' g_prt <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
#' )
#' d_study(g_prt, universe_sizes = c(F1 = 3, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
#'
#' @export
d_study <- function(g_study_result, universe_sizes, d_sizes = NULL,
                    d_design = NULL, d_nesting = NULL) {
  prepared <- .prepare_keyword_d_study(g_study_result, universe_sizes, d_design, d_nesting)
  design <- prepared$design
  if (!is.null(d_sizes)) {
    d_sizes <- .fill_fixed_facet_sizes(
      d_sizes, "d_sizes", prepared$fixed_sizes, g_study_result$facet_columns
    )
    d_sizes <- .check_d_study_sizes(
      prepared$universe_sizes, d_sizes, prepared$facet_names, g_study_result$facet_columns
    )
  } else if (!g_study_result$balanced && design$d_design != g_study_result$g_design) {
    stop(
      "Without d_sizes, the D study repeats the G-study pattern of these ",
      "unbalanced data, so d_design must be the G-study design (",
      g_study_result$g_design, "). Give d_sizes for a balanced D study in another design.",
      call. = FALSE
    )
  }

  general <- d_study_general(
    g_study_result$general,
    object_of_measurement = "p",
    universe_sizes = .symbol_sizes(prepared$universe_sizes),
    d_sizes = .symbol_sizes(d_sizes),
    d_design = design$string
  )
  return(.keyword_d_study_result(general, prepared, g_study_result))
}

#' @export
print.d_study <- function(x, ..., digits = 7) {
  cat("\n")
  cat("G-study design: ", x$g_design, "\n", sep = "")
  cat("D-study design: ", x$d_design, "\n", sep = "")
  if (!is.null(x$d_nesting)) {
    cat("D-study nesting: ", x$d_nesting, "\n", sep = "")
  }
  cat("Facets: ", paste(x$facet_labels, collapse = "; "), "\n", sep = "")
  .print_d_study_tables(x, digits)
  return(invisible(x))
}

#' Compare D Studies at Different Sample Sizes
#'
#' Evaluates combinations of facet sample sizes from one G study,
#' with one D-study design and one set of universe sizes.
#'
#' @inheritParams d_study
#' @param d_sizes D-study sample sizes: a named list of numeric vectors for
#'   all combinations (F1 varies fastest), a data frame with one combination
#'   per row, or a named numeric vector for one D study.
#'
#' @returns An object of class `d_study_comparison`, a list with the following
#'   elements:
#' \describe{
#'   \item{\code{scenarios}}{A data frame with one row per D study. Columns
#'     contain the row number (`scenario`), facet sample sizes, estimates of
#'     `tau`, `delta`, `Delta`, `G`, and `Phi`, their `se_` columns,
#'     `se_available`, and `status` (`"ok"`, `"warning"`, or `"error"`) with
#'     the `warning` and `error` messages. See [d_study()] for the
#'     definitions of the variances and coefficients.}
#'   \item{Design information}{\code{g_design}, \code{d_design},
#'     \code{d_nesting}, \code{universe_sizes},
#'     \code{facet_columns}, and \code{fixed_facet}, which apply to all rows.}
#' }
#' @seealso [find_d_designs()] to select D studies that meet a target;
#'   [compare_d_studies_general()] for the result of [g_study_general()].
#'
#' @examples
#' # p x i x o: 4, 8, or 16 items on 1 or 2 occasions (Brennan, 2001,
#' # Table 4.4).
#' g_pio <- g_study(
#'   brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed"
#' )
#' compare_d_studies(
#'   g_pio,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = list(F1 = c(4, 8, 16), F2 = c(1, 2))
#' )
#'
#' # p x (r:t): the listed numbers of tasks and raters per task, with at most
#' # 12 ratings per person (Table 4.6).
#' g_prt <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
#' )
#' compare_d_studies(
#'   g_prt,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = data.frame(F1 = 1:6, F2 = c(12, 6, 4, 3, 2, 2))
#' )
#'
#' # Tasks fixed in the G study: the task sizes can be omitted (Table 4.8).
#' g_prt_fixed <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
#'   fixed = "t"
#' )
#' compare_d_studies(g_prt_fixed, universe_sizes = c(F2 = Inf), d_sizes = list(F2 = 1:4))
#'
#' @export
compare_d_studies <- function(g_study_result, universe_sizes, d_sizes,
                              d_design = NULL, d_nesting = NULL) {
  prepared <- .prepare_keyword_d_study(g_study_result, universe_sizes, d_design, d_nesting)
  design <- prepared$design
  facet_names <- prepared$facet_names
  candidates <- .d_study_candidates(d_sizes, facet_names, prepared$fixed_sizes)
  general <- compare_d_studies_general(
    g_study_result$general,
    object_of_measurement = "p",
    universe_sizes = .symbol_sizes(prepared$universe_sizes),
    d_sizes = .symbol_sizes(candidates),
    d_design = design$string
  )

  scenarios <- general$scenarios
  size_columns <- match(unname(.keyword_symbols[facet_names]), names(scenarios))
  names(scenarios)[size_columns] <- facet_names
  other_columns <- setdiff(names(scenarios), c("scenario", facet_names))
  scenarios <- scenarios[c("scenario", facet_names, other_columns)]

  result <- list(
    scenarios = scenarios,
    g_design = g_study_result$g_design,
    d_design = design$d_design,
    d_nesting = design$d_nesting,
    universe_sizes = prepared$universe_sizes,
    facet_columns = g_study_result$facet_columns,
    fixed_facet = g_study_result$fixed_facet
  )
  class(result) <- "d_study_comparison"
  return(result)
}

#' @export
print.d_study_comparison <- function(x, ..., digits = 7, row.names = FALSE) {
  design_line <- x$d_design
  if (x$d_design != x$g_design) {
    design_line <- paste0(x$g_design, " -> ", x$d_design)
  }
  cat("D-study comparison: ", design_line, "\n", sep = "")
  if (!is.null(x$d_nesting)) {
    cat("Nesting: ", x$d_nesting, " (child sizes are per parent level)\n", sep = "")
  }
  facets <- names(x$universe_sizes)
  for (facet in facets) {
    cat(facet, " = ", x$facet_columns[[facet]], "; N' = ",
        x$universe_sizes[[facet]], "\n", sep = "")
  }
  .print_scenarios(x, facets, digits, row.names, ...)
  return(invisible(x))
}

# The facet columns of a design: `facets` maps each facet symbol to a data
# column, and NULL uses the symbols as the column names.
.resolve_facet_columns <- function(facets, symbols) {
  if (is.null(facets)) {
    result <- symbols
    names(result) <- symbols
    return(result)
  }

  if (!is.character(facets) || is.null(names(facets)) || anyNA(facets) ||
      length(facets) != length(symbols) ||
      !setequal(names(facets), symbols)) {
    stop(
      "facets must be a named character vector that maps each facet symbol of ",
      "the design (", paste(symbols, collapse = ", "), ") to a data column, ",
      "such as c(", symbols[[1]], " = \"", symbols[[1]], "_column\").",
      .case_note(names(facets), symbols),
      call. = FALSE
    )
  }
  if (anyDuplicated(facets) > 0) {
    stop("facets maps two facet symbols to the same data column.", call. = FALSE)
  }

  return(facets[symbols])
}

# The effects of the design, the columns of its facets, and the lmer term of
# each effect: the interaction of its facets. `model_names` names the facets
# in the lmer model, by default their symbols.
.build_g_study_design <- function(design, facets, score_column, model_names = NULL) {
  .check_single_string(score_column, "score_column")
  parsed <- .parse_design(design, "design")
  if (length(parsed$facets) < 2) {
    stop(
      "design must have at least two facets, such as \"p x i\": the residual ",
      "is the interaction of all of them.",
      call. = FALSE
    )
  }

  if (is.null(model_names)) {
    model_names <- parsed$facets
    names(model_names) <- parsed$facets
  }
  effects <- .design_effects(parsed$parents)
  effect_names <- names(effects)
  n_effects <- length(effect_names)
  lmer_terms <- vapply(
    effects,
    function(members) paste(model_names[members], collapse = ":"),
    character(1)
  )

  return(list(
    g_design = parsed$design,
    parents = parsed$parents,
    random_effects = effect_names[-n_effects],
    residual_effect = effect_names[[n_effects]],
    lmer_terms = lmer_terms,
    model_names = model_names[parsed$facets],
    score_column = score_column,
    facet_columns = .resolve_facet_columns(facets, parsed$facets)
  ))
}

# Check `fixed`, the symbols of nonnested facets fixed in the universe of
# admissible observations, and return them in the order of the design.
.resolve_fixed_facets <- function(fixed, design) {
  if (is.null(fixed)) {
    return(character(0))
  }

  facets <- names(design$parents)
  if (!is.character(fixed) || anyNA(fixed) || !all(fixed %in% facets)) {
    stop(
      "fixed must contain facet symbols of the design: ",
      paste(facets, collapse = ", "), ".",
      .case_note(fixed, facets),
      call. = FALSE
    )
  }

  for (facet in fixed) {
    parents <- design$parents[[facet]]
    if (length(parents) > 0) {
      stop(
        facet, " is nested within ", paste(parents, collapse = " and "),
        ", so its conditions differ across the levels of ",
        paste(parents, collapse = " and "), ". Only a nonnested facet can be fixed.",
        call. = FALSE
      )
    }
  }

  fixed_facets <- facets[facets %in% fixed]
  fixed_effects <- .effects_within(design$random_effects, fixed_facets)
  if (length(fixed_effects) == length(design$random_effects)) {
    stop(
      "With ", paste(fixed_facets, collapse = " and "), " fixed, every effect ",
      "other than the residual is a fixed effect, and the model needs at least ",
      "one random effect besides the residual.",
      call. = FALSE
    )
  }
  return(fixed_facets)
}

# The cell of each score among the observed combinations of levels of the
# factor columns, numbered from 1.
.cell_index <- function(factors) {
  result <- rep(1L, nrow(factors))
  n_cells <- 1
  for (column in factors) {
    key <- (as.integer(column) - 1) * n_cells + result
    distinct_keys <- sort(unique(key))
    result <- match(key, distinct_keys)
    n_cells <- length(distinct_keys)
  }
  return(result)
}

# The number of levels of a facet in each cell of the facets it is nested
# within among the observed scores or, for a nonnested facet, its number of
# levels.
.levels_per_parent_cell <- function(observed_design, facet, parent_columns) {
  if (length(parent_columns) == 0) {
    return(nlevels(observed_design[[facet]]))
  }

  # Each distinct pair of parent cell and nested condition counts once.
  pair_id <- .cell_index(observed_design[c(parent_columns, facet)])
  first_of_pair <- !duplicated(pair_id)
  parent_id <- .cell_index(observed_design[first_of_pair, parent_columns, drop = FALSE])
  result <- tabulate(parent_id, nbins = max(parent_id))
  return(result)
}

# G-study sample size of each facet: its number of levels, or for a nested
# facet its number of levels in each combination of levels of the facets it
# is nested within (NA if that varies).
.g_study_sample_sizes <- function(observed_design, parents) {
  result <- vapply(
    names(parents),
    function(facet) {
      counts <- .levels_per_parent_cell(observed_design, facet, parents[[facet]])
      sample_size <- if (length(unique(counts)) == 1) counts[[1]] else NA_real_
      return(sample_size)
    },
    FUN.VALUE = numeric(1)
  )

  return(result)
}

# The number of cells of the design among the observed conditions: the
# combinations of levels in which each facet has a level observed within the
# levels of the facets it is nested within. Facets linked by nesting are
# joined on the facets they share, and the groups so formed are crossed.
.n_design_cells <- function(observed_design, parents) {
  facets <- names(parents)
  group <- seq_along(facets)
  names(group) <- facets
  for (facet in facets) {
    linked <- c(facet, parents[[facet]])
    group[group %in% group[linked]] <- min(group[linked])
  }

  result <- 1
  for (group_id in unique(group)) {
    cells <- NULL
    for (facet in facets[group == group_id]) {
      observed_cells <- unique(observed_design[c(facet, parents[[facet]])])
      cells <- if (is.null(cells)) observed_cells else merge(cells, observed_cells)
    }
    result <- result * nrow(cells)
  }
  return(result)
}

# The sample sizes of the observed scores as printed lines: the number of
# levels of each facet, per level of the facets it is nested within (a range
# for unequal numbers), and, with missing scores, the number of scores
# observed out of the cells of the design.
.sample_size_lines <- function(observed_design, parents, facet_columns) {
  parts <- character(0)
  for (facet in names(parents)) {
    counts <- .levels_per_parent_cell(observed_design, facet, parents[[facet]])
    size <- if (min(counts) == max(counts)) {
      min(counts)
    } else {
      paste(min(counts), "to", max(counts))
    }
    effect_name <- .effect_name(c(facet, parents[[facet]]), parents)
    parts <- c(parts, paste(.effect_label(effect_name, facet_columns), size))
  }
  result <- paste0("Sample sizes: ", paste(parts, collapse = "; "))

  n_cells <- .n_design_cells(observed_design, parents)
  if (nrow(observed_design) < n_cells) {
    result <- c(
      result,
      paste0("Observed scores: ", nrow(observed_design), " of ", n_cells)
    )
  }
  return(result)
}

# One non-missing score in every cell, with the same number of levels of each
# nested facet in every cell of the facets it is nested within.
.is_balanced <- function(observed_design, sample_sizes) {
  result <- !anyNA(sample_sizes) && nrow(observed_design) == prod(sample_sizes)
  return(result)
}

# A fixed facet requires one non-missing score in every cell of a balanced
# design. Brennan's (2001) Eq. 3.34, used below to report the components with
# the facets fixed, rests on the expected mean squares of balanced designs, and
# unequal numbers within a fixed facet make a multivariate design (Sec. 9.1).
.check_fixed_facet_data <- function(dat_for_lmer, balanced) {
  if (!balanced || anyNA(dat_for_lmer$score)) {
    stop(
      "fixed requires complete, balanced data: one non-missing ",
      "score in every cell of the design, with the same number of nested ",
      "conditions in each level of their parents. Unequal numbers within a fixed ",
      "facet call for a multivariate design.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Components for the universe of admissible observations with the facets in M
# fixed (Brennan Eq. 3.34): sigma^2(A | M) sums sigma^2_S over the effects S
# that add facets of M to A, each divided by the product of n over the added
# facets. The covariance is W V W', where W holds these weights. The effects
# whose facets are all in M are quadratic forms of fixed effects, so W has no
# rows for them.
.mixed_model_components <- function(components, fixed_facets, sample_sizes) {
  fitted_effects <- names(components$vc)
  weights <- matrix(
    0,
    nrow = length(fitted_effects),
    ncol = length(fitted_effects),
    dimnames = list(fitted_effects, fitted_effects)
  )

  for (target in fitted_effects) {
    target_facets <- .parse_effect(target)
    for (source in fitted_effects) {
      source_facets <- .parse_effect(source)
      if (!all(target_facets %in% source_facets)) {
        next
      }
      extra_facets <- setdiff(source_facets, target_facets)
      if (all(extra_facets %in% fixed_facets)) {
        weights[target, source] <- 1 / prod(sample_sizes[extra_facets])
      }
    }
  }

  vc <- as.vector(weights %*% components$vc)
  names(vc) <- fitted_effects

  # Each covariance uses only the contributing components, so the NA row of a
  # zero-boundary component reaches only the components it enters.
  vcov <- matrix(
    NA_real_,
    nrow = length(fitted_effects),
    ncol = length(fitted_effects),
    dimnames = list(fitted_effects, fitted_effects)
  )
  for (row_effect in fitted_effects) {
    row_sources <- fitted_effects[weights[row_effect, ] != 0]
    for (column_effect in fitted_effects) {
      column_sources <- fitted_effects[weights[column_effect, ] != 0]
      vcov[row_effect, column_effect] <- as.numeric(
        weights[row_effect, row_sources] %*%
          components$vcov[row_sources, column_sources, drop = FALSE] %*%
          weights[column_effect, column_sources]
      )
    }
  }
  se <- sqrt(diag(vcov))
  names(se) <- fitted_effects

  return(list(vc = vc, se = se, vcov = vcov))
}

# Rename the design columns to score and the facet symbols, with the facets as
# factors.
.prepare_g_study_data <- function(data, design) {
  needed_columns <- c(design$score_column, unname(design$facet_columns))

  if (!all(needed_columns %in% names(data))) {
    missing_columns <- needed_columns[!(needed_columns %in% names(data))]
    # Without `facets`, each facet symbol is the name of a column.
    facets_hint <- ""
    missing_symbols <- intersect(missing_columns, names(design$facet_columns))
    if (all(design$facet_columns == names(design$facet_columns)) &&
        length(missing_symbols) > 0) {
      facets_hint <- paste0(
        ". Each facet symbol of the design is the name of a data column unless ",
        "facets maps it to another column, such as facets = c(",
        missing_symbols[[1]], " = \"", missing_symbols[[1]], "_column\")"
      )
    }
    stop(
      "These columns are missing from data: ", paste(missing_columns, collapse = ", "),
      facets_hint, ".", .case_note(missing_columns, names(data)),
      call. = FALSE
    )
  }

  result <- data[, needed_columns]
  names(result) <- c("score", names(design$facet_columns))

  facet_names <- names(design$facet_columns)
  result[facet_names] <- lapply(result[facet_names], factor)

  # The residual is interpreted as the highest-order interaction, which is
  # only identified with a single observation per full design cell.
  cell_ids <- .cell_index(result[names(design$facet_columns)])
  if (anyDuplicated(cell_ids) > 0) {
    n_replicated <- sum(duplicated(cell_ids))
    stop(
      "Data contain ", n_replicated, " replicated observation(s) within ",
      paste(unname(design$facet_columns), collapse = " x "),
      " cells. The G-study residual is the highest-order interaction, which ",
      "requires exactly one observation per cell. Aggregate replicated ",
      "observations before fitting, or restructure the design so replication ",
      "is an explicit facet.",
      call. = FALSE
    )
  }

  return(result)
}

#' Generalizability (G) Study from a Design String
#'
#' Estimates the variance components of a G study with any number of crossed
#' and nested facets, written as a design string, by REML with [lme4::lmer()]
#' or from the ANOVA expected mean squares.
#'
#' @param data A data frame in long format, with one row per score and no
#'   replicated design cells. Missing scores are dropped.
#' @param design The design, written with one-letter facet symbols, `x` for
#'   crossing, `:` for nesting, and parentheses, such as `"p x i"`, `"i:p"`,
#'   or `"p x (r:t)"` (see Details).
#' @param facets A named character vector that maps each facet symbol to a
#'   data column, such as `c(p = "person", i = "item")`. The default, `NULL`,
#'   uses the symbols as the column names.
#' @param score_column Name of the score column.
#' @param fixed Optional symbols of nonnested facets fixed in the universe of
#'   admissible observations, such as `"t"` or `c("b", "c")`. Requires
#'   complete, balanced data.
#' @param method `"reml"` for a random-intercept linear mixed model fitted by
#'   REML, with standard errors from the observed information, or `"anova"`
#'   for estimates from the expected mean squares (see ANOVA estimates).
#'
#' @details
#' In `design`, `A:B` nests every facet of A within every facet of B, and a
#' chain `i:h:g` nests i within h and h within g; the other pairs of facets
#' are crossed. A symbol may appear more than once and denotes the same facet
#' each time: in `"p x (a:b) x (c:d:b)"`, a and d are nested within b. Symbols
#' are letters other than x and X, which both mean crossing, and their case
#' matters, as in the column names: `"P X I"` uses the columns `P` and `I`.
#' The other arguments name facets by the same symbols.
#'
#' The components are those of the design's effects (Brennan, 2001, Sec. 3.2),
#' named in Brennan's notation: the primary indices, then the nesting indices
#' of each level after a colon, each level in the order of the facets in
#' `design`. The last component is the residual, which combines the
#' highest-order interaction and error.
#'
#' | `design` | Components |
#' |---|---|
#' | `"p x i"` | p, i, pi |
#' | `"i:p"` | p, i:p |
#' | `"p x i x o"` | p, i, o, pi, po, io, pio |
#' | `"p x (r:t)"` | p, t, r:t, pt, pr:t |
#' | `"r:t:p"` | p, t:p, r:t:p |
#' | `"(i:p) x o"` | p, o, i:p, po, io:p |
#'
#' The object of measurement is chosen in [d_study_general()], so one G study serves
#' any nonnested random facet as the object of measurement.
#'
#' The effects whose facets are all declared in `fixed` enter the model as
#' fixed effects, and the reported components are those of the universe of
#' admissible observations with these facets fixed (Brennan, 2001, Eq. 3.34).
#' The fixed facets keep their observed conditions in every D study. To fix a
#' facet only in the D study, leave it random here and set its sample size
#' equal to its universe size in [d_study_general()].
#'
#' @section ANOVA estimates:
#' With `method = "anova"`, the variance components solve the expected mean
#' square equations of the random effects model (Brennan, 2001, Eqs. 3.24 and
#' 3.25). For unbalanced data, with missing scores or unequal numbers of
#' nested conditions, the sums of squares and expected mean squares are the
#' analogous ones of Henderson's Method 1 (Sec. 7.1.1), as in urGENOVA.
#' Estimates below zero are set to zero, leaving the other estimates unchanged
#' (Sec. 3.4.6); with `fixed`, Eq. 3.34 is applied before that step. The
#' sampling variances and covariances of the estimates are estimated without
#' bias under normality, which for balanced data is Brennan's Eq. 6.1, with
#' df + 2, as in GENOVA. A component whose estimated sampling variance is
#' below zero has an `NA` standard error. For balanced data with every
#' estimate above zero, REML gives the same estimates, with standard errors
#' larger by a factor of up to sqrt((df + 2) / df).
#'
#' @returns An object of class `g_study_general`, a list with the following
#'   elements:
#' \describe{
#'   \item{\code{variance_components}}{A data frame of reported variance
#'     components and standard errors, with effect labels built from the
#'     data's column names. With `fixed`, the rows of the effects whose facets
#'     are all fixed are `NA`.}
#'   \item{\code{variance_component}, \code{variance_component_se},
#'     \code{variance_component_vcov}}{The fitted variance components,
#'     their standard errors, and the covariance matrix used by
#'     [d_study_general()].}
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
#'   \item{Design information}{\code{g_design}, the design string;
#'     \code{parents}, the facets each facet is nested within;
#'     \code{facet_columns}, \code{formula}, \code{fixed_facets}, and
#'     \code{sample_sizes}, the G-study n of each facet (for a nested facet,
#'     per level of its parents; `NA` if that varies).}
#'   \item{\code{observed_design}}{The facet columns of the observed scores,
#'     the pattern used by [d_study_general()] for unbalanced data.}
#'   \item{\code{balanced}}{Whether every cell holds one score and each nested
#'     facet has the same number of levels in every cell of the facets it is
#'     nested within.}
#' }
#' Elements for the other method are `NULL`.
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [d_study_general()]; [g_study()] for the one- and two-facet
#'   designs named by keyword.
#'
#' @examples
#' # p x i (Brennan, 2001, Table 2.3).
#' g_study_general(brennan_p_x_i, design = "p x i")
#'
#' # i:p, a different set of items for each person (Table 2.9).
#' g_study_general(brennan_i_in_p, design = "i:p")
#'
#' # p x i x o (Table 3.3).
#' g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
#'
#' # p x (r:t), raters nested within tasks (Table 3.4).
#' g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
#'
#' # The same data with other column names, mapped to the facet symbols.
#' ratings <- brennan_p_x_r_in_t
#' names(ratings) <- c("person", "task", "rater", "rating")
#' g_study_general(
#'   ratings, design = "p x (r:t)",
#'   facets = c(p = "person", r = "rater", t = "task"), score_column = "rating"
#' )
#'
#' # Tasks fixed in the universe of admissible observations (Appendix F).
#' g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t")
#'
#' # ANOVA estimates, with GENOVA's standard errors (Appendix F).
#' g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", method = "anova")
#'
#' # p x (i:h) with 2, 4, and 2 items in the levels of h (Table 7.5), and
#' # the analogous ANOVA estimates (Sec. 7.1.3).
#' g_study_general(brennan_p_x_i_in_h_unbalanced, design = "p x (i:h)", method = "anova")
#'
#' @export
g_study_general <- function(data, design, facets = NULL, score_column = "score",
                            fixed = NULL, method = c("reml", "anova")) {
  method <- match.arg(method)
  return(.fit_g_study(data, design, facets, score_column, fixed, method))
}

# The G study of g_study_general(). g_study() also gives `model_names`.
.fit_g_study <- function(data, design, facets, score_column, fixed, method,
                         model_names = NULL) {
  design <- .build_g_study_design(design, facets, score_column, model_names)
  fixed_facets <- .resolve_fixed_facets(fixed, design)
  dat_for_lmer <- .prepare_g_study_data(data, design)
  observed <- dat_for_lmer[!is.na(dat_for_lmer$score), ]
  observed_design <- droplevels(observed[names(design$parents)])
  rownames(observed_design) <- NULL
  sample_sizes <- .g_study_sample_sizes(observed_design, design$parents)
  balanced <- .is_balanced(observed_design, sample_sizes)
  effect_names <- c(design$random_effects, design$residual_effect)
  effect_labels <- vapply(
    effect_names,
    .effect_label,
    character(1),
    facet_columns = design$facet_columns,
    USE.NAMES = FALSE
  )

  # The effects whose facets are all fixed in the G study enter the model as
  # fixed effects; the other effects stay random.
  fixed_effects <- .effects_within(design$random_effects, fixed_facets)
  if (length(fixed_facets) > 0) {
    .check_fixed_facet_data(dat_for_lmer, balanced)
    design$random_effects <- setdiff(design$random_effects, fixed_effects)
  }
  fitted_effects <- c(design$random_effects, design$residual_effect)

  formula <- NULL
  fit <- NULL
  fit_diagnostics <- NULL
  anova_table <- NULL
  ems_coefficients <- NULL
  if (method == "reml") {
    formula <- .lmer_formula(
      "score", design$lmer_terms[design$random_effects], design$lmer_terms[fixed_effects]
    )
    model_data <- dat_for_lmer
    names(model_data) <- c("score", unname(design$model_names[names(design$facet_columns)]))
    fit_result <- .fit_lmer_with_diagnostics(formula, model_data)
    fit <- fit_result$fit
    fit_diagnostics <- fit_result$diagnostics
    components <- .extract_variance_components(fit, design)
  } else {
    anova <- .anova_table(observed, effect_names)
    .check_anova_df(anova, effect_labels)
    ems_coefficients <- .ems_coefficients(observed, effect_names, anova$df)
    estimates <- .anova_components(anova, ems_coefficients, observed, balanced, effect_labels)
    components <- list(
      vc = estimates$vc[fitted_effects],
      se = estimates$se[fitted_effects],
      vcov = estimates$vcov[fitted_effects, fitted_effects, drop = FALSE]
    )
    anova_table <- data.frame(
      effect = effect_names,
      effect_label = effect_labels,
      anova[c("df", "t_term", "sum_of_squares", "mean_square")],
      estimate = unname(estimates$vc[effect_names]),
      row.names = NULL
    )
  }

  reported <- components
  if (length(fixed_facets) > 0) {
    reported <- .mixed_model_components(components, fixed_facets, sample_sizes)
  }

  # ANOVA estimates below zero are set to zero after Eq. 3.34, which uses them
  # as computed (Brennan 2001, pp. 84 and 90). Their standard errors are those
  # of the estimates as computed.
  if (method == "anova") {
    components$vc <- pmax(components$vc, 0)
    reported$vc <- pmax(reported$vc, 0)
  }

  vc_table <- data.frame(
    effect = effect_names,
    effect_label = effect_labels,
    variance_estimate = unname(reported$vc[effect_names]),
    standard_error = unname(reported$se[effect_names]),
    row.names = NULL
  )

  # variance_component holds the fitted components, the input to
  # d_study_general(); variance_components is the reported table.
  result <- list(
    fit = fit,
    fit_diagnostics = fit_diagnostics,
    variance_component = components$vc,
    variance_component_se = components$se,
    variance_component_vcov = components$vcov,
    method = method,
    anova_table = anova_table,
    ems_coefficients = ems_coefficients,
    formula = formula,
    g_design = design$g_design,
    parents = design$parents,
    facet_columns = design$facet_columns,
    fixed_facets = fixed_facets,
    sample_sizes = sample_sizes,
    observed_design = observed_design,
    balanced = balanced,
    variance_components = vc_table
  )
  class(result) <- "g_study_general"
  return(result)
}

#' @export
print.g_study_general <- function(x, ..., digits = 7) {
  cat("\n")
  cat("========================================================\n")
  cat("G study: ", x$g_design, "\n", sep = "")
  cat("========================================================\n\n")
  facets <- names(x$facet_columns)
  if (any(unname(x$facet_columns) != facets)) {
    facet_labels <- vapply(facets, .facet_display, character(1), facet_columns = x$facet_columns)
    cat("Facets: ", paste(facet_labels, collapse = ", "), "\n", sep = "")
  }
  cat(.sample_size_lines(x$observed_design, x$parents, x$facet_columns), sep = "\n")
  cat("\n")
  .print_estimation(x, digits)
  cat("\n")
  if (length(x$fixed_facets) == 0) {
    cat("G-study variance components:\n")
  } else {
    fixed_labels <- paste0(x$fixed_facets, " (N = n = ", x$sample_sizes[x$fixed_facets], ")")
    cat(
      "G-study variance components with ", paste(fixed_labels, collapse = " and "),
      " fixed:\n",
      sep = ""
    )
  }
  print(x$variance_components, digits = digits, row.names = FALSE)
  fixed_effects <- .effects_within(x$variance_components$effect, x$fixed_facets)
  if (length(fixed_effects) > 0) {
    cat("Fixed effects, whose rows are NA: ", paste(fixed_effects, collapse = ", "), "\n",
        sep = "")
  }
  .print_estimation_notes(x, fixed_effects)
  return(invisible(x))
}

# The printed estimation step of a G study: the lmer formula, or the ANOVA
# table and expected mean squares.
.print_estimation <- function(x, digits) {
  if (x$method == "reml") {
    cat("Model formula:\n")
    cat(paste(deparse(x$formula, width.cutoff = 500), collapse = "\n"), "\n", sep = "")
    return(invisible(NULL))
  }
  cat("ANOVA table:\n")
  print(x$anova_table, digits = digits, row.names = FALSE)
  cat("\n")
  effect_labels <- x$anova_table$effect_label
  names(effect_labels) <- x$anova_table$effect
  .print_expected_mean_squares(x$ems_coefficients, effect_labels)
  return(invisible(NULL))
}

# The lmer fit diagnostics, or a note when an ANOVA estimate of a random
# effect, one not in `fixed_effects`, is below zero.
.print_estimation_notes <- function(x, fixed_effects) {
  if (x$method == "reml") {
    .print_lmer_fit_diagnostics(x$fit_diagnostics)
    return(invisible(NULL))
  }
  random_rows <- !(x$anova_table$effect %in% fixed_effects)
  if (any(x$anova_table$estimate[random_rows] < 0)) {
    cat("Estimates below zero are set to zero; the ANOVA table gives them as computed.\n")
  }
  return(invisible(NULL))
}

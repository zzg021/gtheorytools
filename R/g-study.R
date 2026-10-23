# The components of the chosen design and the columns mapped to P, F1, and F2.
.build_g_study_design <- function(score_column, p, f1, f2 = NULL, design) {
  .check_single_string(score_column, "score_column")
  .check_single_string(p, "p")
  .check_single_string(f1, "f1")
  .check_single_string(design, "design")

  one_facet_designs <- c("crossed", "nested")
  two_facet_designs <- c("fully_crossed", "partially_nested", "fully_nested")

  if (is.null(f2)) {
    if (!(design %in% one_facet_designs)) {
      stop(
        "For one-facet G-studies with f2 = NULL, design must be one of: crossed, nested.",
        call. = FALSE
      )
    }

    design_specs <- list(
      crossed = list(
        label = "One-facet crossed G-study: P x F1",
        random_effects = c("P", "F1"),
        residual_effect = "PxF1",
        n_facets = 1L
      ),
      nested = list(
        label = "One-facet nested G-study: F1:P",
        random_effects = c("P"),
        residual_effect = "F1:P",
        n_facets = 1L
      )
    )

    result <- design_specs[[design]]
    result$g_design <- design
    result$score_column <- score_column
    result$facet_columns <- c(P = p, F1 = f1)
    return(result)
  }

  .check_single_string(f2, "f2")

  if (!(design %in% two_facet_designs)) {
    stop(
      "For two-facet G-studies, design must be one of: ",
      "fully_crossed, partially_nested, fully_nested.",
      call. = FALSE
    )
  }

  # In the nested designs F2 is nested within F1.
  design_specs <- list(
    fully_crossed = list(
      label = "Two-facet fully crossed G-study: P x F1 x F2",
      random_effects = c("P", "F1", "F2", "PxF1", "PxF2", "F1xF2"),
      residual_effect = "PxF1xF2",
      n_facets = 2L
    ),
    partially_nested = list(
      label = "Two-facet partially nested G-study: P x F1 x (F2:F1)",
      random_effects = c("P", "F1", "PxF1", "F2:F1"),
      residual_effect = "F2:F1xP",
      n_facets = 2L
    ),
    fully_nested = list(
      label = "Two-facet fully nested G-study: F2:F1:P",
      random_effects = c("P", "F1:P"),
      residual_effect = "F2:F1:P",
      n_facets = 2L
    )
  )

  result <- design_specs[[design]]
  result$g_design <- design
  result$score_column <- score_column
  result$facet_columns <- c(P = p, F1 = f1, F2 = f2)
  return(result)
}

# Resolve `fixed` to F1 or F2 and check that the design supports it.
.resolve_fixed_facet <- function(fixed, design) {
  if (is.null(fixed)) {
    return(NULL)
  }

  .check_single_string(fixed, "fixed")
  facet <- names(design$facet_columns)[design$facet_columns == fixed]
  if (length(facet) != 1 || facet == "P") {
    stop("fixed must name the f1 or f2 column.", call. = FALSE)
  }

  if (design$g_design == "fully_crossed") {
    return(facet)
  }
  if (design$g_design == "partially_nested") {
    if (facet == "F2") {
      stop(
        "In the partially nested design only f1 can be fixed: the levels of ",
        fixed, " differ across the levels of ", design$facet_columns[["F1"]], ".",
        call. = FALSE
      )
    }
    return(facet)
  }
  if (design$g_design == "crossed") {
    stop(
      "fixed is not available for the crossed design: fixing its only facet ",
      "leaves no random facet.",
      call. = FALSE
    )
  }

  stop(
    "fixed is not available for the ", design$g_design, " design: a facet nested ",
    "within persons has no common set of conditions.",
    call. = FALSE
  )
}

# G-study sample size of each facet: its number of levels, or for a nested
# facet its number of levels in each level of the parent (NA if that varies).
.g_study_sample_sizes <- function(dat_for_lmer, design) {
  parents <- switch(
    design$g_design,
    crossed = list(F1 = character(0)),
    nested = list(F1 = "P"),
    fully_crossed = list(F1 = character(0), F2 = character(0)),
    partially_nested = list(F1 = character(0), F2 = "F1"),
    fully_nested = list(F1 = "P", F2 = c("P", "F1"))
  )
  observed <- dat_for_lmer[!is.na(dat_for_lmer$score), ]
  result <- vapply(
    names(parents),
    function(facet) {
      parent_columns <- parents[[facet]]
      if (length(parent_columns) == 0) {
        return(length(unique(observed[[facet]])))
      }

      parent_id <- interaction(observed[parent_columns], drop = TRUE)
      counts <- tapply(observed[[facet]], parent_id, function(x) length(unique(x)))
      sample_size <- if (length(unique(counts)) == 1) counts[[1]] else NA_real_
      return(sample_size)
    },
    FUN.VALUE = numeric(1)
  )

  return(result)
}

# A fixed facet requires one non-missing score in every cell of a balanced
# design. Brennan's (2001) Eq. 3.34, used below to report the components with
# the facet fixed, rests on the expected mean squares of balanced designs, and
# unequal numbers within a fixed facet make a multivariate design (Sec. 9.1).
.check_fixed_facet_data <- function(dat_for_lmer, sample_sizes, fixed) {
  n_persons <- length(unique(dat_for_lmer$P))
  expected_rows <- n_persons * prod(sample_sizes)

  if (anyNA(sample_sizes) || anyNA(dat_for_lmer$score) ||
      nrow(dat_for_lmer) != expected_rows) {
    stop(
      "fixed = \"", fixed, "\" requires complete, balanced data: one non-missing ",
      "score for every person in every condition, with the same number of nested ",
      "conditions in each level of their parent. Unequal numbers within a fixed ",
      "facet call for a multivariate design.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Components for the universe of admissible observations with facet F fixed
# (Brennan Eq. 3.34): sigma^2(A | M) is sigma^2_A plus, when A does not contain
# F, the component of the effect that adds F to A divided by n_F. The
# covariance is W V W', where W holds these weights. The main effect of F is a
# quadratic form of fixed effects, so W has no row for it.
.mixed_model_components <- function(components, fixed_facet, n_fixed) {
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
      if (length(extra_facets) == 0) {
        weights[target, source] <- 1
      } else if (identical(extra_facets, fixed_facet)) {
        weights[target, source] <- 1 / n_fixed
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

# Rename the design columns to score, P, F1, and F2, with the facets as factors.
.prepare_g_study_data <- function(data, design) {
  needed_columns <- c(design$score_column, unname(design$facet_columns))

  if (!all(needed_columns %in% names(data))) {
    missing_columns <- needed_columns[!(needed_columns %in% names(data))]
    stop(
      "These columns are missing from data: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  result <- data[, needed_columns]
  names(result) <- c("score", names(design$facet_columns))

  facet_names <- names(design$facet_columns)
  result[facet_names] <- lapply(result[facet_names], factor)

  # The residual is interpreted as the highest-order interaction, which is
  # only identified with a single observation per full design cell.
  cell_ids <- interaction(result[names(design$facet_columns)], drop = TRUE)
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

#' Generalizability (G) Study with a Linear Mixed Model
#'
#' Estimates the variance components of a one- or two-facet G study, with
#' persons as the object of measurement, by fitting a random-intercept linear
#' mixed model with [lme4::lmer()] (REML). Standard errors use the REML
#' observed information.
#'
#' @param data A data frame in long format, with one row per score and no
#'   replicated design cells. Missing scores are dropped when `fixed = NULL`.
#' @param score_column Name of the score column.
#' @param p Name of the person (object of measurement) column.
#' @param f1,f2 Names of the facet columns. `f2 = NULL` gives a one-facet
#'   design. In the nested designs `f2` is nested within `f1`; swap the
#'   columns for the other order.
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
#' nested within F1. The last component listed is the residual, which
#' combines the highest-order interaction and error.
#'
#' | `design` | Design | Components |
#' |---|---|---|
#' | `crossed` | P x F1 | P, F1, PxF1 |
#' | `nested` | F1:P | P, F1:P |
#' | `fully_crossed` | P x F1 x F2 | P, F1, F2, PxF1, PxF2, F1xF2, PxF1xF2 |
#' | `partially_nested` | (F2:F1) x P | P, F1, PxF1, F2:F1, F2:F1xP |
#' | `fully_nested` | F2:F1:P | P, F1:P, F2:F1:P |
#'
#' A facet declared with `fixed` keeps its observed conditions in every
#' D study. To fix a facet only in the D study, leave it random here and
#' set its sample size equal to its universe size in [d_study()].
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
#'   \item{\code{fit}}{The fitted `lmer` model.}
#'   \item{\code{fit_diagnostics}}{Singularity and convergence information,
#'     with messages and warnings from `lmer`.}
#'   \item{Design information}{\code{g_design}, \code{design_label},
#'     \code{formula}, \code{n_facets}, \code{facet_columns},
#'     \code{fixed_facet}, and \code{sample_sizes}.}
#' }
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [d_study()], [compare_d_studies()], and the example data
#'   [brennan_p_x_i], [brennan_i_in_p], [brennan_p_x_i_x_o], and
#'   [brennan_p_x_r_in_t].
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
#' @export
g_study <- function(data, score_column = "score", p, f1, f2 = NULL, design,
                    fixed = NULL) {
  design <- .build_g_study_design(score_column, p, f1, f2, design)
  fixed_facet <- .resolve_fixed_facet(fixed, design)
  dat_for_lmer <- .prepare_g_study_data(data, design)
  sample_sizes <- .g_study_sample_sizes(dat_for_lmer, design)
  effect_names <- c(design$random_effects, design$residual_effect)

  # A facet fixed in the G study enters the model as a fixed effect; its
  # interactions and the effects nested in it stay random.
  if (!is.null(fixed_facet)) {
    .check_fixed_facet_data(dat_for_lmer, sample_sizes, fixed)
    design$random_effects <- setdiff(design$random_effects, fixed_facet)
  }

  formula <- .lmer_formula("score", design$random_effects, fixed_facet)
  fit_result <- .fit_lmer_with_diagnostics(formula, dat_for_lmer)
  fit <- fit_result$fit
  fit_diagnostics <- fit_result$diagnostics

  components <- .extract_variance_components(fit, design)
  reported <- components
  if (!is.null(fixed_facet)) {
    reported <- .mixed_model_components(components, fixed_facet, sample_sizes[[fixed_facet]])
  }

  vc_table <- data.frame(
    effect = effect_names,
    effect_label = vapply(
      effect_names,
      .effect_label,
      character(1),
      facet_columns = design$facet_columns
    ),
    variance_estimate = unname(reported$vc[effect_names]),
    standard_error = unname(reported$se[effect_names]),
    row.names = NULL
  )

  # variance_component holds the fitted components, the input to d_study();
  # variance_components is the reported table.
  result <- list(
    fit = fit,
    fit_diagnostics = fit_diagnostics,
    variance_component = components$vc,
    variance_component_se = components$se,
    variance_component_vcov = components$vcov,
    design_label = design$label,
    formula = formula,
    g_design = design$g_design,
    n_facets = design$n_facets,
    facet_columns = design$facet_columns,
    fixed_facet = fixed_facet,
    sample_sizes = sample_sizes,
    variance_components = vc_table
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
  cat("Model formula:\n")
  cat(paste(deparse(x$formula, width.cutoff = 500), collapse = "\n"), "\n", sep = "")
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
  .print_lmer_fit_diagnostics(x$fit_diagnostics)
  return(invisible(x))
}

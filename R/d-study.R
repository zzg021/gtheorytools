# The printed facets of a D study: random, finite, or fixed, with their sizes.
# `shown_names` gives how each facet is shown, and d_sizes = NULL is the
# G-study pattern of unbalanced data.
.d_facet_labels <- function(facet_status, universe_sizes, d_sizes, fixed_facets,
                            shown_names) {
  result <- character(length(facet_status))

  for (k in seq_along(facet_status)) {
    facet <- names(facet_status)[k]
    status <- facet_status[[facet]]
    if (is.null(d_sizes)) {
      result[k] <- paste0(shown_names[[facet]], " random (G-study pattern)")
      next
    }
    sizes <- paste0("n' = ", d_sizes[[facet]])
    if (status == "finite") {
      sizes <- paste0(sizes, ", N' = ", universe_sizes[[facet]])
    }
    if (status == "fixed") {
      sizes <- paste0("n' = N' = ", d_sizes[[facet]])
      if (facet %in% fixed_facets) {
        status <- "fixed in the G study"
      }
    }
    result[k] <- paste0(shown_names[[facet]], " ", status, " (", sizes, ")")
  }

  return(result)
}

# The D-study steps that do not depend on d_sizes: checks on the G-study
# result, object of measurement, universe sizes, and design; recoding; and the
# single-observation components. compare_d_studies_general() runs them once
# for all sample sizes.
.prepare_d_study <- function(g_study_result, object_of_measurement, universe_sizes,
                             d_design = NULL) {
  .validate_g_study_result(g_study_result, "g_study_general")
  design <- .resolve_d_study_design(g_study_result, object_of_measurement, d_design)
  object_of_measurement <- design$object_of_measurement
  facet_names <- setdiff(names(g_study_result$parents), object_of_measurement)
  fixed_facets <- g_study_result$fixed_facets
  facet_columns <- g_study_result$facet_columns
  universe_sizes <- .fill_fixed_facet_sizes(
    universe_sizes, "universe_sizes", g_study_result$sample_sizes[fixed_facets],
    facet_columns
  )
  universe_values <- .check_size_vector(
    universe_sizes, "universe_sizes", facet_names, allow_inf = TRUE
  )[facet_names]

  se_available <- all(is.finite(g_study_result$variance_component_vcov))
  # Facets fixed in the G study are fixed in every D study, so they are
  # checked once here for all sizes.
  .check_g_fixed_facet_nesting(fixed_facets, design, facet_columns)
  recoded <- .recode_variance_components(
    g_study_result$variance_component, g_study_result$variance_component_vcov,
    design$d_parents, fixed_facets
  )
  variance_component <- recoded$variance_component
  variance_component_vcov <- recoded$variance_component_vcov
  effects <- names(variance_component)
  effect_labels <- vapply(effects, .effect_label, character(1), facet_columns = facet_columns)
  # The labels of the G-study components summed into each D-study component.
  source_effects <- vapply(
    recoded$component_sources,
    function(sources) {
      labels <- vapply(sources, .effect_label, character(1), facet_columns = facet_columns)
      return(paste(labels, collapse = " + "))
    },
    character(1)
  )

  weights <- .single_obs_weights(effects, universe_values, object_of_measurement)
  single_obs <- matrix(
    weights %*% variance_component, nrow = 1, dimnames = list(NULL, effects)
  )
  single_obs_cov <- NULL
  if (se_available) {
    single_obs_cov <- weights %*% variance_component_vcov %*% t(weights)
  }
  single_obs_results <- .component_table(
    single_obs, single_obs_cov, "single_obs_estimate", effect_labels, source_effects
  )

  return(list(
    g_study_result = g_study_result,
    facet_names = facet_names,
    universe_sizes = universe_values,
    design = design,
    se_available = se_available,
    component_mapping = recoded$component_mapping,
    variance_component = variance_component,
    variance_component_vcov = variance_component_vcov,
    effect_labels = effect_labels,
    source_effects = source_effects,
    object_of_measurement = object_of_measurement,
    single_obs = single_obs,
    single_obs_cov = single_obs_cov,
    single_obs_results = single_obs_results
  ))
}

# The D study that repeats the G-study pattern of unbalanced data: its
# mean-score multipliers and relative-error weights. Unequal numbers of
# conditions make the finite-universe factors (1 - n'/N') differ across
# levels, so every facet is random, and the pattern belongs to the G-study
# design.
.pattern_d_study <- function(prepared) {
  g_study_result <- prepared$g_study_result
  design <- prepared$design
  same_nesting <- vapply(
    names(design$g_parents),
    function(facet) setequal(design$d_parents[[facet]], design$g_parents[[facet]]),
    logical(1)
  )
  if (!all(same_nesting)) {
    stop(
      "Without d_sizes, the D study repeats the G-study pattern of these ",
      "unbalanced data, so d_design must be the G-study design (",
      design$g_design, "). Give d_sizes for a balanced D study in another design.",
      call. = FALSE
    )
  }
  if (any(is.finite(prepared$universe_sizes))) {
    stop(
      "Without d_sizes, the D study repeats the G-study pattern of these ",
      "unbalanced data, which requires every universe size to be Inf. Give ",
      "d_sizes for a balanced D study with finite universes.",
      call. = FALSE
    )
  }

  object_of_measurement <- prepared$object_of_measurement
  multipliers <- .pattern_multipliers(
    g_study_result$observed_design, names(prepared$variance_component),
    object_of_measurement
  )
  delta_weight <- multipliers$relative / multipliers$absolute
  delta_weight[object_of_measurement] <- 0
  facet_status <- rep("random", length(prepared$facet_names))
  names(facet_status) <- prepared$facet_names

  return(list(
    multiplier = multipliers$absolute,
    delta_weight = delta_weight,
    facet_status = facet_status
  ))
}

# The D-study steps for one set of sizes, given the result of .prepare_d_study().
# d_sizes = NULL repeats the G-study pattern of unbalanced data.
.evaluate_d_study <- function(prepared, d_sizes) {
  g_study_result <- prepared$g_study_result
  universe_values <- prepared$universe_sizes
  object_of_measurement <- prepared$object_of_measurement

  if (is.null(d_sizes)) {
    pattern <- .pattern_d_study(prepared)
    d_values <- NULL
    facet_status <- pattern$facet_status
    multiplier <- pattern$multiplier
    delta_weight <- pattern$delta_weight
  } else {
    fixed_facets <- g_study_result$fixed_facets
    d_sizes <- .fill_fixed_facet_sizes(
      d_sizes, "d_sizes", g_study_result$sample_sizes[fixed_facets],
      g_study_result$facet_columns
    )
    d_values <- .check_d_study_sizes(
      universe_values, d_sizes, prepared$facet_names, g_study_result$facet_columns
    )
    facet_status <- .d_facet_status(universe_values, d_values)
    .check_fixed_facets(facet_status, d_values, g_study_result, prepared$design)
    multiplier <- .d_mean_multiplier(
      prepared$single_obs, universe_values, d_values, prepared$design$d_parents,
      object_of_measurement
    )
    delta_weight <- NULL
  }
  shown_names <- vapply(
    prepared$facet_names, .facet_display, character(1),
    facet_columns = g_study_result$facet_columns
  )
  facet_labels <- .d_facet_labels(
    facet_status, universe_values, d_values, g_study_result$fixed_facets, shown_names
  )

  d_mean <- .d_mean_vc(prepared$single_obs, multiplier)
  d_mean_cov <- NULL
  if (!is.null(prepared$single_obs_cov)) {
    d_mean_cov <- .d_mean_vcov(prepared$single_obs, prepared$single_obs_cov, multiplier)
  }

  d_mean_results <- .component_table(
    d_mean, d_mean_cov, "mean_score_estimate", prepared$effect_labels,
    prepared$source_effects
  )
  if (is.null(d_sizes)) {
    d_mean_results$divisor <- unname(1 / multiplier[d_mean_results$effect])
    d_mean_results$delta_weight <- unname(delta_weight[d_mean_results$effect])
  }
  coefficient_results <- .g_coefficient(
    d_mean, d_mean_cov, object_of_measurement, delta_weight
  )

  result <- list(
    g_design = prepared$design$g_design,
    d_design = prepared$design$d_design,
    object_of_measurement = object_of_measurement,
    component_mapping = prepared$component_mapping,
    d_variance_component = prepared$variance_component,
    d_variance_component_vcov = prepared$variance_component_vcov,
    single_obs = prepared$single_obs,
    single_obs_vcov = prepared$single_obs_cov,
    single_obs_results = prepared$single_obs_results,
    d_mean = d_mean,
    d_mean_vcov = d_mean_cov,
    d_mean_results = d_mean_results,
    coefficient_results = coefficient_results,
    universe_sizes = universe_values,
    d_sizes = d_values,
    facet_status = facet_status,
    facet_labels = facet_labels,
    se_available = prepared$se_available
  )
  class(result) <- "d_study_general"
  return(result)
}

#' Decision (D) Study from `g_study_general()`
#'
#' Computes D-study variance components, error variances, and the
#' generalizability (G) and dependability (Phi) coefficients for one object
#' of measurement and one set of facet sample sizes. Standard errors use the
#' G-study covariance matrix and the delta method.
#'
#' @param g_study_result The object returned by [g_study_general()].
#' @param object_of_measurement The symbol of the facet whose levels are
#'   measured, usually `"p"` for persons: a nonnested facet that is random in
#'   the G study.
#' @param universe_sizes Universe sizes N' of the other facets, named by
#'   facet symbol: whole numbers or `Inf` for an infinite universe. Facets
#'   fixed in the G study may be omitted.
#' @param d_sizes D-study sample sizes n', whole numbers named as in
#'   `universe_sizes`, with `1 <= n' <= N'`. For a nested facet, n' is the
#'   number of conditions in each level of the facets it is nested within.
#'   `NULL`, the default, repeats the G study: its sample sizes for balanced
#'   data, and its pattern of scores for unbalanced data (see Unbalanced data).
#' @param d_design The D-study design, a design string with the facets of the
#'   G-study design, written as `design` in [g_study_general()], such as
#'   `"p x (i:o)"`. Defaults to the G-study design, and may be more nested:
#'   each facet is nested within at least the facets it is nested within in
#'   the G study.
#'
#' @section Coefficients:
#' The output labels `tau`, `delta`, and `Delta` denote universe-score
#' variance \eqn{\sigma^2(\tau)}, relative error variance
#' \eqn{\sigma^2(\delta)}, and absolute error variance
#' \eqn{\sigma^2(\Delta)}, respectively. The coefficients are
#' \deqn{G = \frac{\sigma^2(\tau)}{\sigma^2(\tau) + \sigma^2(\delta)},
#'   \qquad
#'   \Phi = \frac{\sigma^2(\tau)}{\sigma^2(\tau) + \sigma^2(\Delta)}}{
#'   G = tau / (tau + delta); Phi = tau / (tau + Delta).}
#' G is the generalizability coefficient for relative decisions, such as
#' ranking persons. Phi is the dependability coefficient for absolute
#' decisions, such as comparing scores with a cut score.
#'
#' @details
#' A facet is random when N' is infinite, finite when `n' < N' < Inf`, and
#' fixed in the D study when `n' = N'`. D-study effects are named as in
#' [g_study_general()] (`pr:t`). Recoding to a more nested design adds each
#' G-study component to the D-study effect that contains its facets and the
#' facets they are nested within; `component_mapping` lists the sums.
#'
#' @section Unbalanced data:
#' For unbalanced data, with every universe size `Inf` and `d_design` the
#' G-study design, `d_sizes = NULL` gives the D study with the G-study
#' pattern: the same numbers of conditions and the same missing scores, as
#' when a new form keeps a test's unequal numbers of items per passage
#' (Brennan, 2001, Sec. 7.2). The mean score of each level of the object of
#' measurement, such as a person, weights each cell of an effect by its share
#' of that level's scores, and each mean-score component is the variance
#' component divided by the `divisor` column: the reciprocal of the mean, over
#' the levels, of the sum of the squared weights. With persons as the object
#' of measurement, the divisor is the harmonic mean of the numbers of items
#' per person in i:p and in p x i with missing scores, and \eqn{n_{i+}} or
#' \eqn{n_{i+}^2 / \sum_h n_{i:h}^2} in p x (i:h) (Eqs. 7.18, 7.27-7.29, and
#' 7.41). `delta_weight` gives the share of each component in relative error,
#' for which \eqn{\sigma^2(\tau) + \sigma^2(\delta)} is the expected variance
#' of the mean scores. Effects without the object of measurement enter
#' relative error when the levels have different patterns, as with missing
#' scores (Eq. 7.45).
#'
#' @returns An object of class `d_study_general`, a list with the following
#'   elements:
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
#'     \code{object_of_measurement}, \code{component_mapping},
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
#' @seealso [compare_d_studies_general()] for several sets of sizes.
#'
#' @examples
#' # p x i with 12 items (Brennan, 2001, Table 2.4).
#' g_pi <- g_study_general(brennan_p_x_i, design = "p x i")
#' d_study_general(
#'   g_pi, object_of_measurement = "p", universe_sizes = c(i = Inf), d_sizes = c(i = 12)
#' )
#'
#' # p x i with missing scores (Table 7.7), and the D study that repeats the
#' # G-study pattern (Sec. 7.2.4).
#' g_missing <- g_study_general(
#'   brennan_p_x_i_missing, design = "p x i", method = "anova"
#' )
#' d_study_general(
#'   g_missing, object_of_measurement = "p", universe_sizes = c(i = Inf)
#' )
#'
#' # p x i x o with 4 items and 2 occasions (Table 4.4).
#' g_pio <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
#' d_study_general(
#'   g_pio, object_of_measurement = "p",
#'   universe_sizes = c(i = Inf, o = Inf), d_sizes = c(i = 4, o = 2)
#' )
#'
#' # The same sample sizes with different items on each occasion, p x (i:o)
#' # (Table 4.5).
#' d_study_general(
#'   g_pio, object_of_measurement = "p",
#'   universe_sizes = c(i = Inf, o = Inf), d_sizes = c(i = 4, o = 2),
#'   d_design = "p x (i:o)"
#' )
#'
#' # p x (r:t) with the 3 tasks fixed and 4 raters per task (Table 4.8).
#' g_prt <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
#' d_study_general(
#'   g_prt, object_of_measurement = "p",
#'   universe_sizes = c(r = Inf, t = 3), d_sizes = c(r = 4, t = 3)
#' )
#'
#' @export
d_study_general <- function(g_study_result, object_of_measurement, universe_sizes,
                            d_sizes = NULL, d_design = NULL) {
  prepared <- .prepare_d_study(g_study_result, object_of_measurement, universe_sizes, d_design)
  if (is.null(d_sizes) && g_study_result$balanced) {
    d_sizes <- g_study_result$sample_sizes[prepared$facet_names]
    exceeds_universe <- d_sizes > prepared$universe_sizes[names(d_sizes)]
    if (any(exceeds_universe)) {
      facet <- names(d_sizes)[exceeds_universe][1]
      stop(
        "Without d_sizes, the D study uses the G-study sample sizes, and ",
        g_study_result$facet_columns[[facet]], " has n = ", d_sizes[[facet]],
        " in the G study, more than its universe size, N' = ",
        prepared$universe_sizes[[facet]], ". Give d_sizes, or a universe size of ",
        "at least ", d_sizes[[facet]], ".",
        call. = FALSE
      )
    }
  }
  return(.evaluate_d_study(prepared, d_sizes))
}

#' @export
print.d_study_general <- function(x, ..., digits = 7) {
  cat("\n")
  cat("G-study design: ", x$g_design, "\n", sep = "")
  cat("D-study design: ", x$d_design, "\n", sep = "")
  cat("Object of measurement: ", x$object_of_measurement, "\n", sep = "")
  cat("Facets: ", paste(x$facet_labels, collapse = "; "), "\n", sep = "")
  .print_d_study_tables(x, digits)
  return(invisible(x))
}

# The printed D-study tables, after a note when standard errors are
# unavailable.
.print_d_study_tables <- function(x, digits) {
  if (!x$se_available) {
    cat("\n")
    cat("Note: G-study standard errors are unavailable (a REML variance component\n")
    cat("at the zero boundary, an unusable observed-information matrix, or an\n")
    cat("ANOVA sampling variance below zero), so all D-study standard errors are\n")
    cat("reported as NA.\n")
  }
  cat("\n")
  cat("D-study single-observation variance components:\n")
  print(x$single_obs_results, digits = digits, row.names = FALSE)

  cat("\n")
  cat("D-study mean-score variance components:\n")
  print(x$d_mean_results, digits = digits, row.names = FALSE)

  cat("\n")
  cat("G coefficient and dependability:\n")
  print(x$coefficient_results, digits = digits, row.names = FALSE)
  return(invisible(NULL))
}

.d_facet_labels <- function(facet_status, universe_sizes, d_sizes, g_study_result) {
  universe_values <- unlist(universe_sizes)
  d_values <- unlist(d_sizes)
  result <- character(length(facet_status))

  for (k in seq_along(facet_status)) {
    facet <- names(facet_status)[k]
    status <- facet_status[[facet]]
    sizes <- paste0("n' = ", d_values[[facet]])
    if (status == "finite") {
      sizes <- paste0(sizes, ", N' = ", universe_values[[facet]])
    }
    if (status == "fixed") {
      sizes <- paste0("n' = N' = ", d_values[[facet]])
      if (identical(facet, g_study_result$fixed_facet)) {
        status <- "fixed in the G study"
      }
    }
    result[k] <- paste0(g_study_result$facet_columns[[facet]], " ", status, " (", sizes, ")")
  }

  return(result)
}

.effect_label_lookup <- function(g_study_result) {
  result <- g_study_result$variance_components$effect_label
  names(result) <- g_study_result$variance_components$effect
  return(result)
}

.target_effect_label_lookup <- function(effect_names, facet_columns) {
  result <- vapply(
    effect_names,
    .effect_label,
    character(1),
    facet_columns = facet_columns
  )
  names(result) <- effect_names
  return(result)
}

.source_effect_label_lookup <- function(component_sources, source_effect_labels) {
  result <- vapply(
    component_sources,
    function(source_names) {
      paste(unname(source_effect_labels[source_names]), collapse = " + ")
    },
    character(1)
  )
  names(result) <- names(component_sources)
  return(result)
}

.add_effect_metadata <- function(result_table, effect_labels, source_effects) {
  result_table$effect_label <- unname(effect_labels[result_table$effect])
  result_table$source_effects <- unname(source_effects[result_table$effect])
  column_order <- c(
    "effect",
    "effect_label",
    "source_effects",
    setdiff(names(result_table), c("effect", "effect_label", "source_effects"))
  )
  return(result_table[, column_order])
}

# The D-study steps that do not depend on d_sizes: checks on the G-study
# result, universe sizes, and design; recoding; and the single-observation
# components. compare_d_studies() runs them once for all sample sizes.
.prepare_d_study <- function(g_study_result, universe_sizes,
                             d_design = NULL, d_nesting = NULL) {
  .validate_g_study_result(g_study_result)
  n_facets <- g_study_result$n_facets
  facet_names <- .d_study_facet_names(n_facets)
  fixed_facet <- g_study_result$fixed_facet
  if (!is.null(fixed_facet)) {
    universe_sizes <- .fill_fixed_facet_sizes(
      universe_sizes, "universe_sizes", fixed_facet,
      g_study_result$sample_sizes[[fixed_facet]],
      g_study_result$facet_columns[[fixed_facet]]
    )
  }
  universe_values <- .check_size_vector(
    universe_sizes, "universe_sizes", facet_names, allow_inf = TRUE
  )[facet_names]

  variance_component <- g_study_result$variance_component
  variance_component_vcov <- g_study_result$variance_component_vcov
  se_available <- all(is.finite(variance_component_vcov))
  design <- .resolve_d_study_design(
    g_design = g_study_result$g_design,
    d_design = d_design,
    d_nesting = d_nesting,
    n_facets = n_facets
  )
  # A facet fixed in the G study is fixed in every D study, so it is checked
  # once here for all sizes.
  if (!is.null(fixed_facet)) {
    .check_fixed_nested_within_p(fixed_facet, design, g_study_result$facet_columns)
  }
  recoded_components <- .recode_variance_components(
    variance_component = variance_component,
    variance_component_vcov = variance_component_vcov,
    g_design = design$g_design,
    d_design = design$d_design,
    d_nesting = design$d_nesting,
    fixed_facet = fixed_facet,
    facet_columns = g_study_result$facet_columns
  )
  variance_component <- recoded_components$variance_component
  variance_component_vcov <- recoded_components$variance_component_vcov
  effect_labels <- .target_effect_label_lookup(
    names(variance_component),
    g_study_result$facet_columns
  )
  source_effects <- .source_effect_label_lookup(
    recoded_components$component_sources,
    .effect_label_lookup(g_study_result)
  )
  object_of_measurement <- "P"

  single_obs <- .single_obs_vc(variance_component, universe_values, object_of_measurement)
  single_obs_cov <- NULL

  if (se_available) {
    single_obs_cov <- .single_obs_vcov(
      variance_component,
      universe_values,
      variance_component_vcov,
      object_of_measurement
    )
  }
  single_obs_results <- .single_obs_vc_table(single_obs, single_obs_cov)
  single_obs_results <- .add_effect_metadata(
    single_obs_results, effect_labels, source_effects
  )

  return(list(
    g_study_result = g_study_result,
    facet_names = facet_names,
    universe_sizes = universe_values,
    design = design,
    se_available = se_available,
    recoded_components = recoded_components,
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

# The D-study steps for one set of sizes, given the result of .prepare_d_study().
.evaluate_d_study <- function(prepared, d_sizes) {
  g_study_result <- prepared$g_study_result
  universe_values <- prepared$universe_sizes
  object_of_measurement <- prepared$object_of_measurement
  fixed_facet <- g_study_result$fixed_facet
  if (!is.null(fixed_facet)) {
    d_sizes <- .fill_fixed_facet_sizes(
      d_sizes, "d_sizes", fixed_facet,
      g_study_result$sample_sizes[[fixed_facet]],
      g_study_result$facet_columns[[fixed_facet]]
    )
  }
  d_values <- .check_d_study_sizes(universe_values, d_sizes, prepared$facet_names)
  facet_status <- .d_facet_status(universe_values, d_values)
  .check_fixed_facets(facet_status, d_values, g_study_result, prepared$design)

  d_mean <- .d_mean_vc(
    prepared$variance_component, prepared$single_obs, universe_values, d_values,
    object_of_measurement
  )
  d_mean_cov <- NULL
  if (!is.null(prepared$single_obs_cov)) {
    d_mean_cov <- .d_mean_vcov(
      prepared$variance_component,
      prepared$single_obs,
      prepared$single_obs_cov,
      universe_values,
      d_values,
      object_of_measurement
    )
  }

  d_mean_results <- .d_mean_vc_table(d_mean, d_mean_cov)
  d_mean_results <- .add_effect_metadata(
    d_mean_results, prepared$effect_labels, prepared$source_effects
  )
  coefficient_results <- .g_coefficient(d_mean, d_mean_cov, object_of_measurement)

  result <- list(
    g_design = prepared$design$g_design,
    d_design = prepared$design$d_design,
    d_nesting = prepared$design$d_nesting,
    n_facets = g_study_result$n_facets,
    component_mapping = prepared$recoded_components$component_mapping,
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
    facet_labels = .d_facet_labels(facet_status, universe_values, d_values, g_study_result),
    se_available = prepared$se_available
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
#' @param g_study_result The object returned by [g_study()].
#' @param universe_sizes Universe sizes N', named `F1`, or `F1` and `F2`:
#'   whole numbers or `Inf` for an infinite universe. A facet fixed in the
#'   G study may be omitted.
#' @param d_sizes D-study sample sizes n', whole numbers named as in
#'   `universe_sizes`, with `1 <= n' <= N'`. For a nested facet, n' is the
#'   number of conditions in each level of its parent.
#' @param d_design The D-study design, one of the `design` values of
#'   [g_study()]. Defaults to the G-study design; it may be more nested but
#'   not less.
#' @param d_nesting `"F2:F1"` or `"F1:F2"`, required when a fully crossed
#'   G study is recoded to a nested design. A nested G study fixes it as
#'   `"F2:F1"`.
#'
#' @details
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
#'     components and standard errors.}
#'   \item{\code{d_variance_component}}{Variance components after recoding
#'     to the D-study design.}
#'   \item{\code{single_obs}, \code{d_mean}}{The single-observation and
#'     mean-score component estimates as one-row matrices.}
#'   \item{\code{d_variance_component_vcov}, \code{single_obs_vcov},
#'     \code{d_mean_vcov}}{Covariance matrices for the corresponding component
#'     estimates. `single_obs_vcov` and `d_mean_vcov` are `NULL` when
#'     standard errors are unavailable.}
#'   \item{Design information}{\code{g_design}, \code{d_design},
#'     \code{d_nesting}, \code{n_facets}, \code{component_mapping},
#'     \code{universe_sizes}, \code{d_sizes}, \code{facet_status}, and
#'     \code{facet_labels}.}
#'   \item{\code{se_available}}{Whether the G-study covariance matrix can
#'     be used to calculate standard errors.}
#' }
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' @seealso [compare_d_studies()] for several sets of sizes.
#'
#' @examples
#' # p x I with 12 items (Brennan, 2001, Table 2.4).
#' g_pi <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed")
#' d_study(g_pi, universe_sizes = c(F1 = Inf), d_sizes = c(F1 = 12))
#'
#' # p x I x O with 4 items and 2 occasions (Table 4.4).
#' g_pio <- g_study(
#'   brennan_p_x_i_x_o, p = "p", f1 = "i", f2 = "o", design = "fully_crossed"
#' )
#' d_study(g_pio, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 4, F2 = 2))
#'
#' # The same sample sizes with different items on each occasion, p x (I:O)
#' # (Table 4.5).
#' d_study(
#'   g_pio,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = c(F1 = 4, F2 = 2),
#'   d_design = "partially_nested",
#'   d_nesting = "F1:F2"
#' )
#'
#' # p x (R:T) with the 3 tasks fixed and 4 raters per task (Table 4.8).
#' g_prt <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
#' )
#' d_study(g_prt, universe_sizes = c(F1 = 3, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
#'
#' @export
d_study <- function(g_study_result, universe_sizes, d_sizes,
                    d_design = NULL, d_nesting = NULL) {
  prepared <- .prepare_d_study(g_study_result, universe_sizes, d_design, d_nesting)
  return(.evaluate_d_study(prepared, d_sizes))
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
  if (!x$se_available) {
    cat("\n")
    cat("Note: G-study standard errors are unavailable (variance component at\n")
    cat("the zero boundary, or an unusable observed-information matrix), so\n")
    cat("all D-study standard errors are reported as NA.\n")
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
  return(invisible(x))
}

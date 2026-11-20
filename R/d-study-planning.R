# A named list gives a Cartesian grid, the first of `facets` varying fastest;
# a data frame gives explicit rows. `fixed_sizes` holds the G-study n of each
# facet fixed in the G study. Size values are checked by .evaluate_d_study()
# for each row, so an invalid size is reported in its own row.
.d_study_candidates <- function(d_sizes, facets, fixed_sizes) {
  fixed_facets <- names(fixed_sizes)
  if (is.numeric(d_sizes) && is.null(dim(d_sizes))) {
    d_sizes <- as.list(d_sizes)
  }
  if (!is.list(d_sizes) || is.null(names(d_sizes)) ||
      anyDuplicated(names(d_sizes)) ||
      !all(names(d_sizes) %in% facets) ||
      !all(setdiff(facets, fixed_facets) %in% names(d_sizes))) {
    stop("d_sizes must name each facet (", paste(facets, collapse = ", "),
         ") exactly once; facets fixed in the G study may be omitted.",
         .case_note(names(d_sizes), facets),
         call. = FALSE)
  }
  for (facet in names(d_sizes)) {
    values <- d_sizes[[facet]]
    if (!is.numeric(values) || length(values) == 0L) {
      stop("Each element of d_sizes must be a nonempty numeric vector.", call. = FALSE)
    }
  }
  for (fixed_facet in setdiff(fixed_facets, names(d_sizes))) {
    d_sizes[[fixed_facet]] <- fixed_sizes[[fixed_facet]]
  }
  d_sizes <- d_sizes[facets]
  if (is.data.frame(d_sizes)) {
    candidates <- as.data.frame(d_sizes)
  } else {
    candidates <- expand.grid(d_sizes, KEEP.OUT.ATTRS = FALSE)
  }
  rownames(candidates) <- NULL
  return(candidates)
}

#' Compare D Studies from `g_study_general()` at Different Sample Sizes
#'
#' Evaluates combinations of facet sample sizes from one G study, with one
#' object of measurement, one D-study design, and one set of universe sizes.
#'
#' @inheritParams d_study_general
#' @param d_sizes D-study sample sizes, named by facet symbol: a named list
#'   of numeric vectors for all combinations (the first facet of the design
#'   varies fastest), a data frame with one combination per row, or a named
#'   numeric vector for one D study.
#'
#' @returns An object of class `d_study_comparison_general`, a list with the
#'   following elements:
#' \describe{
#'   \item{\code{scenarios}}{A data frame with one row per D study. Columns
#'     contain the row number (`scenario`), facet sample sizes, estimates of
#'     `tau`, `delta`, `Delta`, `G`, and `Phi`, their `se_` columns,
#'     `se_available`, and `status` (`"ok"`, `"warning"`, or `"error"`) with
#'     the `warning` and `error` messages. See [d_study_general()] for the
#'     definitions of the variances and coefficients. From a G study with
#'     `method = "bayes"`, the estimates are posterior medians, the `se_`
#'     columns posterior standard deviations, and `lower_G`, `upper_G`,
#'     `lower_Phi`, and `upper_Phi` the bounds of the credible intervals.}
#'   \item{\code{coefficient_draws}, \code{level}}{From a G study with
#'     `method = "bayes"`, a list with one matrix of posterior draws of
#'     `tau`, `delta`, `Delta`, `G`, and `Phi` per D study, `NULL` for a D
#'     study that failed, and the level of the credible intervals; otherwise
#'     `NULL`.}
#'   \item{Design information}{\code{g_design}, \code{d_design},
#'     \code{object_of_measurement}, \code{universe_sizes},
#'     \code{facet_columns}, and \code{fixed_facets}, which apply to all rows.}
#' }
#' @seealso [find_d_designs()] to select D studies that meet a target.
#'
#' @examples
#' # p x i x o: 4, 8, or 16 items on 1 or 2 occasions (Brennan, 2001,
#' # Table 4.4).
#' g_pio <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
#' compare_d_studies_general(
#'   g_pio, object_of_measurement = "p",
#'   universe_sizes = c(i = Inf, o = Inf),
#'   d_sizes = list(i = c(4, 8, 16), o = c(1, 2))
#' )
#'
#' # p x (r:t): the listed numbers of tasks and raters per task, with at most
#' # 12 ratings per person (Table 4.6).
#' g_prt <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
#' compare_d_studies_general(
#'   g_prt, object_of_measurement = "p",
#'   universe_sizes = c(r = Inf, t = Inf),
#'   d_sizes = data.frame(r = c(12, 6, 4, 3, 2, 2), t = 1:6)
#' )
#'
#' # Tasks fixed in the G study: the task sizes can be omitted (Table 4.8).
#' g_prt_fixed <- g_study_general(
#'   brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t"
#' )
#' compare_d_studies_general(
#'   g_prt_fixed, object_of_measurement = "p",
#'   universe_sizes = c(r = Inf), d_sizes = list(r = 1:4)
#' )
#'
#' @export
compare_d_studies_general <- function(g_study_result, object_of_measurement,
                                      universe_sizes, d_sizes, d_design = NULL) {
  prepared <- .prepare_d_study(g_study_result, object_of_measurement, universe_sizes, d_design)
  design <- prepared$design
  fixed_facets <- g_study_result$fixed_facets
  candidates <- .d_study_candidates(
    d_sizes, prepared$facet_names, g_study_result$sample_sizes[fixed_facets]
  )

  n_candidates <- nrow(candidates)
  quantities <- c("tau", "delta", "Delta", "G", "Phi")
  scenarios <- data.frame(scenario = seq_len(n_candidates), candidates)
  for (quantity in c(quantities, paste0("se_", quantities))) {
    scenarios[[quantity]] <- rep(NA_real_, n_candidates)
  }
  scenarios$se_available <- rep(FALSE, n_candidates)
  has_draws <- !is.null(prepared$single_obs_draws)
  coefficient_draws <- NULL
  if (has_draws) {
    for (column in c("lower_G", "upper_G", "lower_Phi", "upper_Phi")) {
      scenarios[[column]] <- rep(NA_real_, n_candidates)
    }
    coefficient_draws <- vector("list", n_candidates)
  }
  scenarios$status <- rep("ok", n_candidates)
  scenarios$warning <- rep("", n_candidates)
  scenarios$error <- rep("", n_candidates)

  warning_messages <- character(0)
  record_warning <- function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  for (i in seq_len(n_candidates)) {
    candidate_sizes <- unlist(candidates[i, , drop = FALSE])
    warning_messages <- character(0)
    result <- tryCatch(
      withCallingHandlers(
        .evaluate_d_study(prepared, candidate_sizes),
        warning = record_warning
      ),
      error = function(e) e
    )
    scenarios$warning[i] <- paste(unique(warning_messages), collapse = "\n")
    if (inherits(result, "error")) {
      scenarios$status[i] <- "error"
      scenarios$error[i] <- conditionMessage(result)
      next
    }
    if (length(warning_messages) > 0L) {
      scenarios$status[i] <- "warning"
    }
    estimates <- result$coefficient_results
    for (quantity in quantities) {
      k <- match(quantity, estimates$quantity)
      scenarios[[quantity]][i] <- estimates$estimate[k]
      scenarios[[paste0("se_", quantity)]][i] <- estimates$standard_error[k]
    }
    scenarios$se_available[i] <- all(is.finite(estimates$standard_error))
    if (has_draws) {
      draw_summary <- .draw_summary(result$coefficient_draws, prepared$level)
      for (quantity in c("G", "Phi")) {
        scenarios[[paste0("lower_", quantity)]][i] <- draw_summary$lower[[quantity]]
        scenarios[[paste0("upper_", quantity)]][i] <- draw_summary$upper[[quantity]]
      }
      coefficient_draws[[i]] <- result$coefficient_draws
    }
  }

  n_errors <- sum(scenarios$status == "error")
  n_warnings <- sum(nzchar(scenarios$warning))
  if (n_errors > 0L || n_warnings > 0L) {
    warning(n_errors, " of ", n_candidates, " D-study scenarios failed; ",
            n_warnings, " produced warnings. See $scenarios$error and $scenarios$warning.",
            call. = FALSE)
  }
  return(structure(
    list(
      scenarios = scenarios,
      g_design = design$g_design,
      d_design = design$d_design,
      object_of_measurement = design$object_of_measurement,
      universe_sizes = prepared$universe_sizes,
      facet_columns = g_study_result$facet_columns,
      fixed_facets = g_study_result$fixed_facets,
      coefficient_draws = coefficient_draws,
      level = prepared$level
    ),
    class = "d_study_comparison_general"
  ))
}

#' @export
print.d_study_comparison_general <- function(x, ..., digits = 7, row.names = FALSE) {
  design_line <- x$d_design
  if (x$d_design != x$g_design) {
    design_line <- paste0(x$g_design, " -> ", x$d_design)
  }
  cat("D-study comparison: ", design_line, "\n", sep = "")
  cat("Object of measurement: ", x$object_of_measurement, "\n", sep = "")
  if (grepl(":", x$d_design, fixed = TRUE)) {
    cat("Sizes of nested facets are per level of the facets they are nested within.\n")
  }
  facets <- names(x$universe_sizes)
  for (facet in facets) {
    cat(.facet_display(facet, x$facet_columns), ": N' = ",
        x$universe_sizes[[facet]], "\n", sep = "")
  }
  .print_scenarios(x, facets, digits, row.names, ...)
  return(invisible(x))
}

# The printed scenarios of a comparison: sizes, G, Phi, and status when a
# scenario has a warning or error.
.print_scenarios <- function(x, facets, digits, row.names, ...) {
  columns <- c("scenario", facets, "G", "Phi")
  if (!is.null(x$coefficient_draws)) {
    columns <- c("scenario", facets, "G", "lower_G", "upper_G", "Phi", "lower_Phi", "upper_Phi")
  }
  if (any(x$scenarios$status != "ok")) {
    columns <- c(columns, "status")
  }
  print(x$scenarios[, columns, drop = FALSE], digits = digits, row.names = row.names, ...)
  if (!is.null(x$coefficient_draws)) {
    cat("G and Phi are posterior medians, with ", format(100 * x$level),
        "% equal-tailed credible intervals.\n", sep = "")
  }
  n_missing_se <- sum(x$scenarios$status != "error" & !x$scenarios$se_available)
  if (n_missing_se > 0L) {
    cat("SEs are unavailable for some or all quantities in ", n_missing_se,
        " scenarios.\n", sep = "")
  }
  return(invisible(NULL))
}

#' Select D Studies That Meet a Coefficient Target
#'
#' Returns the D studies in `comparison` whose estimated G or Phi is at
#' least `target`. From a G study with `method = "bayes"`, the estimate is the
#' posterior median, and `probability` instead selects the D studies whose
#' posterior probability of meeting the target is at least `probability`.
#'
#' @param comparison The object returned by [compare_d_studies()] or
#'   [compare_d_studies_general()].
#' @param coefficient `"G"` or `"Phi"` (see [d_study()]).
#' @param target A number between 0 and 1. An estimate equal to the target up
#'   to rounding error, `sqrt(.Machine$double.eps)`, qualifies.
#' @param probability Optional number between 0 and 1: the minimum posterior
#'   probability that the coefficient is at least `target`. Requires a
#'   comparison from a G study with `method = "bayes"`. A draw equal to the
#'   target up to the rounding error above counts as meeting it. The default,
#'   `NULL`, selects by the estimate.
#'
#' @returns An object of class `d_study_selection`, a list with the following
#'   elements:
#' \describe{
#'   \item{\code{designs}}{The rows of `comparison$scenarios` that meet the
#'     target, in their original order.}
#'   \item{\code{coefficient}}{The selected coefficient, `"G"` or `"Phi"`.}
#'   \item{\code{target}}{The minimum coefficient value.}
#'   \item{\code{meets_target}}{Whether each D study meets the target; `NA`
#'     when the coefficient has no finite estimate.}
#'   \item{\code{comparison}}{The input comparison object.}
#'   \item{\code{probability}}{The minimum posterior probability, or `NULL`.}
#'   \item{\code{posterior_probability}}{From a G study with
#'     `method = "bayes"`, the posterior probability that each D study meets
#'     the target, also a column of `designs`; otherwise `NULL`.}
#' }
#'
#' @examples
#' # p x (r:t) designs with at most 12 ratings per person (Brennan, 2001,
#' # Table 4.6) whose G coefficient is at least .55.
#' g_prt <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
#' )
#' comparison <- compare_d_studies(
#'   g_prt,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = data.frame(F1 = 1:6, F2 = c(12, 6, 4, 3, 2, 2))
#' )
#' find_d_designs(comparison, coefficient = "G", target = 0.55)
#'
#' # From Bayesian estimates, which need brms and a C++ toolchain, the designs
#' # whose G coefficient is at least .55 with posterior probability .50 or
#' # more.
#' \dontrun{
#' g_bayes <- g_study(
#'   brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested",
#'   method = "bayes", bayes_options = list(seed = 1)
#' )
#' comparison <- compare_d_studies(
#'   g_bayes,
#'   universe_sizes = c(F1 = Inf, F2 = Inf),
#'   d_sizes = data.frame(F1 = 1:6, F2 = c(12, 6, 4, 3, 2, 2))
#' )
#' find_d_designs(comparison, coefficient = "G", target = 0.55, probability = 0.5)
#' }
#'
#' @export
find_d_designs <- function(comparison, coefficient, target, probability = NULL) {
  if (!inherits(comparison, c("d_study_comparison", "d_study_comparison_general"))) {
    stop(
      "comparison must be the result of compare_d_studies() or ",
      "compare_d_studies_general().",
      call. = FALSE
    )
  }
  if (length(coefficient) != 1L || !(coefficient %in% c("G", "Phi"))) {
    stop("coefficient must be \"G\" or \"Phi\".", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L ||
      !is.finite(target) || target < 0 || target > 1) {
    stop("target must be a single finite number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(probability)) {
    if (!is.numeric(probability) || length(probability) != 1L ||
        !is.finite(probability) || probability < 0 || probability > 1) {
      stop("probability must be a single finite number between 0 and 1.", call. = FALSE)
    }
    if (is.null(comparison$coefficient_draws)) {
      stop(
        "probability requires a comparison from a G study with method = \"bayes\".",
        call. = FALSE
      )
    }
  }
  scenarios <- comparison$scenarios
  values <- scenarios[[coefficient]]
  eligible <- is.finite(values)
  # NA means that the D study failed or its selected coefficient is undefined.
  meets_target <- rep(NA, nrow(scenarios))
  # An estimate equal to the target up to rounding error qualifies; the
  # tolerance is the default of all.equal().
  rounding_tolerance <- sqrt(.Machine$double.eps)
  posterior_probability <- NULL
  if (!is.null(comparison$coefficient_draws)) {
    posterior_probability <- vapply(
      comparison$coefficient_draws,
      function(draws) {
        if (is.null(draws)) {
          return(NA_real_)
        }
        return(mean(draws[, coefficient] >= target - rounding_tolerance))
      },
      numeric(1)
    )
    scenarios$posterior_probability <- posterior_probability
  }
  if (is.null(probability)) {
    meets_target[eligible] <- values[eligible] >= target - rounding_tolerance
  } else {
    eligible <- eligible & !is.na(posterior_probability)
    meets_target[eligible] <- posterior_probability[eligible] >= probability
  }
  designs <- scenarios[which(meets_target), , drop = FALSE]
  rownames(designs) <- NULL
  return(structure(
    list(
      designs = designs,
      coefficient = coefficient,
      target = target,
      meets_target = meets_target,
      comparison = comparison,
      probability = probability,
      posterior_probability = posterior_probability
    ),
    class = "d_study_selection"
  ))
}

#' @export
print.d_study_selection <- function(x, ..., digits = 7, row.names = FALSE) {
  if (is.null(x$probability) && is.null(x$posterior_probability)) {
    cat("Target: estimated ", x$coefficient, " >= ", x$target, "\n", sep = "")
  } else if (is.null(x$probability)) {
    cat("Target: posterior median of ", x$coefficient, " >= ", x$target, "\n", sep = "")
  } else {
    cat("Target: posterior probability that ", x$coefficient, " >= ", x$target,
        " is at least ", x$probability, "\n", sep = "")
  }
  if (nrow(x$designs) == 0L) {
    cat("No D study meets the target.\n")
  } else {
    cat(nrow(x$designs), " of ", nrow(x$comparison$scenarios),
        " designs meet the target.\n", sep = "")
    columns <- c("scenario", names(x$comparison$universe_sizes), "G", "Phi")
    if (!is.null(x$posterior_probability)) {
      columns <- c(
        "scenario", names(x$comparison$universe_sizes),
        "G", "lower_G", "upper_G", "Phi", "lower_Phi", "upper_Phi", "posterior_probability"
      )
    }
    if (any(x$designs$status != "ok")) {
      columns <- c(columns, "status")
    }
    print(x$designs[, columns, drop = FALSE], digits = digits, row.names = row.names, ...)
  }
  n_unassessed <- sum(is.na(x$meets_target))
  if (n_unassessed > 0L) {
    cat("D studies with no finite estimate of the selected coefficient: ",
        n_unassessed, "\n", sep = "")
  }
  return(invisible(x))
}

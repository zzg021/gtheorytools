.check_single_string <- function(value, argument_name) {
  if (!is.character(value) || length(value) != 1 || is.na(value) || value == "") {
    stop(argument_name, " must be a single non-empty string.", call. = FALSE)
  }
  return(invisible(NULL))
}

.check_size_vector <- function(value, argument_name, facet_names, allow_inf) {
  values <- unlist(value)

  if (!is.numeric(values) || length(values) != length(facet_names) ||
      is.null(names(values)) || !setequal(names(values), facet_names)) {
    stop(
      argument_name, " must be a named numeric vector (or list) containing ",
      "exactly: ", paste(facet_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (any(is.na(values))) {
    stop(argument_name, " must not contain missing values.", call. = FALSE)
  }

  if (!allow_inf && any(!is.finite(values))) {
    stop(argument_name, " values must be finite.", call. = FALSE)
  }

  if (any(values < 1)) {
    stop(argument_name, " values must be at least 1.", call. = FALSE)
  }

  finite_values <- values[is.finite(values)]
  if (any(finite_values != round(finite_values))) {
    stop(
      argument_name, " values must be whole numbers",
      if (allow_inf) " (or Inf)" else "", ".",
      call. = FALSE
    )
  }

  return(values)
}

.validate_g_study_result <- function(g_study_result) {
  required_names <- c(
    "variance_component",
    "variance_component_vcov",
    "variance_components",
    "g_design",
    "n_facets",
    "facet_columns"
  )
  missing_names <- setdiff(required_names, names(g_study_result))

  if (length(missing_names) > 0) {
    stop(
      "g_study_result must be the list returned by g_study(). Missing: ",
      paste(missing_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

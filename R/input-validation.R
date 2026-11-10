.check_single_string <- function(value, argument_name) {
  if (!is.character(value) || length(value) != 1 || is.na(value) || value == "") {
    stop(argument_name, " must be a single non-empty string.", call. = FALSE)
  }
  return(invisible(NULL))
}

# A sentence for an error message, naming the names in `given` that match
# names in `valid` only when case is ignored: "Names are case-sensitive: I is
# not i." It is empty when there are none.
.case_note <- function(given, valid) {
  differs_in_case <- !(given %in% valid) & tolower(given) %in% tolower(valid)
  if (!any(differs_in_case)) {
    return("")
  }
  given <- unique(given[differs_in_case])
  matched <- valid[match(tolower(given), tolower(valid))]
  return(paste0(
    " Names are case-sensitive: ",
    paste0(given, " is not ", matched, collapse = ", "), "."
  ))
}

.check_size_vector <- function(value, argument_name, facet_names, allow_inf) {
  values <- unlist(value)

  if (!is.numeric(values) || length(values) != length(facet_names) ||
      is.null(names(values)) || !setequal(names(values), facet_names)) {
    stop(
      argument_name, " must be a named numeric vector (or list) containing ",
      "exactly: ", paste(facet_names, collapse = ", "), ".",
      .case_note(names(values), facet_names),
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

# `g_function` is "g_study" or "g_study_general", the function whose result
# the D-study function takes; it is also the class of that result.
.validate_g_study_result <- function(g_study_result, g_function) {
  if (inherits(g_study_result, g_function)) {
    return(invisible(NULL))
  }

  other_functions <- if (g_function == "g_study") {
    "g_study_general() goes to d_study_general() and compare_d_studies_general()"
  } else {
    "g_study() goes to d_study() and compare_d_studies()"
  }
  stop(
    "g_study_result must be the result of ", g_function, "(). The result of ",
    other_functions, ".",
    call. = FALSE
  )
}

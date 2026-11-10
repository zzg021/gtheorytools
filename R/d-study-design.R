# A facet fixed in the G study keeps its observed conditions in the D study,
# n' = N' = n, and may be omitted from universe_sizes and d_sizes.
# `fixed_sizes` holds the G-study n of each such facet, and `facet_columns`
# the data column of each facet.
.fill_fixed_facet_sizes <- function(sizes, argument_name, fixed_sizes, facet_columns) {
  values <- unlist(sizes)

  for (fixed_facet in names(fixed_sizes)) {
    n_fixed <- fixed_sizes[[fixed_facet]]
    if (!(fixed_facet %in% names(values))) {
      values[fixed_facet] <- n_fixed
      next
    }

    if (!isTRUE(values[[fixed_facet]] == n_fixed)) {
      column <- facet_columns[[fixed_facet]]
      stop(
        column, " is fixed in the G study, so the D study uses its ", n_fixed,
        " observed conditions, and its size in ", argument_name, " must be ",
        n_fixed, ". To sample from the fixed set, fit the G study with ",
        column, " random and give N' in universe_sizes.",
        call. = FALSE
      )
    }
  }

  return(values)
}

# Check one set of D-study sizes; n' cannot exceed a finite N'.
.check_d_study_sizes <- function(universe_values, d_sizes, facet_names, facet_columns) {
  d_values <- .check_size_vector(
    d_sizes, "d_sizes", facet_names, allow_inf = FALSE
  )

  for (facet in facet_names) {
    if (d_values[[facet]] > universe_values[[facet]]) {
      stop(
        "The D-study sample size of ", facet_columns[[facet]], ", n' = ",
        d_values[[facet]], ", exceeds its universe size, N' = ",
        universe_values[[facet]], ".",
        call. = FALSE
      )
    }
  }

  return(d_values[facet_names])
}

# The object of measurement and the D-study design. The object of measurement
# is a facet that is random in the G study and nonnested in the D study
# (Brennan 2001, App. F). A D-study design is admissible when each facet's
# D-study parents include its G-study parents: nested components can then be
# summed into the more nested design.
.resolve_d_study_design <- function(g_study_result, object_of_measurement, d_design) {
  g_parents <- g_study_result$parents
  facets <- names(g_parents)

  .check_single_string(object_of_measurement, "object_of_measurement")
  if (!(object_of_measurement %in% facets)) {
    stop(
      "object_of_measurement must be one of the facet symbols of the G-study ",
      "design: ", paste(facets, collapse = ", "), ".",
      .case_note(object_of_measurement, facets),
      call. = FALSE
    )
  }
  if (object_of_measurement %in% g_study_result$fixed_facets) {
    stop(
      object_of_measurement, " is fixed in the G study, so it cannot be the ",
      "object of measurement.",
      call. = FALSE
    )
  }

  d_string <- g_study_result$g_design
  d_parents <- g_parents
  if (!is.null(d_design)) {
    parsed <- .parse_design(d_design, "d_design")
    if (!setequal(parsed$facets, facets)) {
      stop(
        "d_design must have the facets of the G-study design: ",
        paste(facets, collapse = ", "), ".",
        .case_note(parsed$facets, facets),
        call. = FALSE
      )
    }
    d_string <- parsed$design
    d_parents <- parsed$parents

    for (facet in facets) {
      missing_parents <- setdiff(g_parents[[facet]], d_parents[[facet]])
      if (length(missing_parents) > 0) {
        stop(
          "d_design cannot be less nested than the G-study design: ", facet,
          " is nested within ", paste(missing_parents, collapse = " and "),
          " in the G study but not in d_design. Variance components can be summed ",
          "into a more nested D-study design, but nested components cannot be ",
          "divided into a less nested design.",
          call. = FALSE
        )
      }
    }
  }

  if (length(d_parents[[object_of_measurement]]) > 0) {
    stop(
      "The object of measurement must be a nonnested facet (Brennan 2001, ",
      "App. F), but ", object_of_measurement, " is nested within ",
      paste(d_parents[[object_of_measurement]], collapse = " and "),
      " in the D-study design.",
      call. = FALSE
    )
  }

  return(list(
    g_design = g_study_result$g_design,
    d_design = d_string,
    object_of_measurement = object_of_measurement,
    g_parents = g_parents,
    d_parents = d_parents
  ))
}

# Each G-study component is added to the smallest D-study effect that
# contains its facets: those facets and every facet they are nested within
# in the D-study design. From a crossed G study this is Brennan's (2001,
# p. 63) rule for confounded effects. A D-study effect whose facets are all
# fixed in the G study is a fixed effect and has no component.
.target_component_sources <- function(source_effects, d_parents, fixed_facets) {
  d_effects <- .design_effects(d_parents)
  result <- list()

  for (target_effect in names(d_effects)) {
    target_facets <- d_effects[[target_effect]]
    if (all(target_facets %in% fixed_facets)) {
      next
    }

    is_source <- vapply(
      source_effects,
      function(source_effect) {
        source_facets <- .parse_effect(source_effect)
        closure <- union(source_facets, unlist(d_parents[source_facets]))
        return(setequal(closure, target_facets))
      },
      logical(1)
    )
    result[[target_effect]] <- unname(source_effects[is_source])
  }

  return(result)
}

# A facet fixed in the G study is one shared set of conditions, so the D-study
# design keeps its G-study parents.
.check_g_fixed_facet_nesting <- function(fixed_facets, design, facet_columns) {
  for (fixed_facet in fixed_facets) {
    new_parents <- setdiff(design$d_parents[[fixed_facet]], design$g_parents[[fixed_facet]])
    if (length(new_parents) > 0) {
      stop(
        facet_columns[[fixed_facet]], " is fixed in the G study, so it is one shared ",
        "set of conditions and cannot be nested within ",
        paste(facet_columns[new_parents], collapse = " and "), " in the D study.",
        call. = FALSE
      )
    }
  }
  return(invisible(NULL))
}

# Each facet is random (N' = Inf), finite (n' < N' < Inf), or fixed (n' = N').
.d_facet_status <- function(universe_sizes, d_sizes) {
  result <- vapply(
    names(d_sizes),
    function(facet) {
      if (is.infinite(universe_sizes[[facet]])) {
        return("random")
      }
      if (d_sizes[[facet]] == universe_sizes[[facet]]) {
        return("fixed")
      }
      return("finite")
    },
    FUN.VALUE = character(1)
  )

  return(result)
}

# A facet with n' = N' is one set of conditions common to all objects of
# measurement, so it is crossed with the object of measurement.
.check_fixed_nested_within_object <- function(fixed_facets, design, facet_columns) {
  object <- design$object_of_measurement
  object_column <- facet_columns[[object]]
  d_parents <- design$d_parents
  nested_within_object <- names(d_parents)[vapply(
    d_parents,
    function(parents) object %in% parents,
    logical(1)
  )]
  for (facet in intersect(nested_within_object, fixed_facets)) {
    stop(
      facet_columns[[facet]], " is fixed (n' = N') and nested within ", object_column,
      " in the D-study design. A fixed facet is one set of conditions common to ",
      "all levels of ", object_column, ", so it is crossed with ", object_column, ".",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Check restrictions and approximations for facets fixed in the D study.
.check_fixed_facets <- function(facet_status, d_sizes, g_study_result, design) {
  fixed_facets <- names(facet_status)[facet_status == "fixed"]

  if (length(fixed_facets) == length(facet_status)) {
    stop(
      "Every facet has n' = N', so there is no error of measurement. ",
      "At least one facet needs n' < N'.",
      call. = FALSE
    )
  }

  # Facets fixed in the G study were checked by .check_g_fixed_facet_nesting().
  facet_columns <- g_study_result$facet_columns
  d_only_fixed <- setdiff(fixed_facets, g_study_result$fixed_facets)
  .check_fixed_nested_within_object(d_only_fixed, design, facet_columns)

  # Facets fixed only in the D study, whose conditions should be the G-study ones.
  for (facet in d_only_fixed) {
    if (!isTRUE(d_sizes[[facet]] == g_study_result$sample_sizes[[facet]])) {
      warning(
        facet_columns[[facet]], " is fixed at n' = N' = ", d_sizes[[facet]],
        ", and the G study ",
        "observed ", g_study_result$sample_sizes[[facet]], " conditions, so the ",
        "results are approximations (Brennan 2001, Sec. 4.3).",
        call. = FALSE
      )
    }
  }

  # Facets fixed only in the D study and nested there within a facet they are
  # crossed with in the G study. Nesting within the object of measurement was
  # rejected above.
  for (facet in d_only_fixed) {
    new_parents <- setdiff(design$d_parents[[facet]], design$g_parents[[facet]])
    if (length(new_parents) > 0) {
      parent_label <- paste(facet_columns[new_parents], collapse = " and ")
      warning(
        facet_columns[[facet]], " is fixed and nested within ", parent_label,
        " in the D study but crossed with it in the G study, so each level of ",
        parent_label, " needs its own conditions and the results are ",
        "approximations (Brennan 2001, Sec. 4.3).",
        call. = FALSE
      )
    }
  }
  return(invisible(NULL))
}

# The D-study components are sums of G-study components, W sigma^2 with W
# marking the sources of each, and their covariance is W V W'.
.recode_variance_components <- function(variance_component, variance_component_vcov,
                                        d_parents, fixed_facets) {
  g_effects <- names(variance_component)
  component_sources <- .target_component_sources(g_effects, d_parents, fixed_facets)
  d_effects <- names(component_sources)
  weights <- matrix(
    0,
    nrow = length(d_effects),
    ncol = length(g_effects),
    dimnames = list(d_effects, g_effects)
  )
  for (effect in d_effects) {
    weights[effect, component_sources[[effect]]] <- 1
  }

  d_variance_component <- as.numeric(weights %*% variance_component)
  names(d_variance_component) <- d_effects
  d_variance_component_vcov <- weights %*%
    variance_component_vcov[g_effects, g_effects, drop = FALSE] %*%
    t(weights)

  component_mapping <- data.frame(
    effect = d_effects,
    source_effects = vapply(component_sources, paste, character(1), collapse = " + "),
    row.names = NULL
  )

  return(list(
    variance_component = d_variance_component,
    variance_component_vcov = d_variance_component_vcov,
    component_sources = component_sources,
    component_mapping = component_mapping
  ))
}

.d_study_facet_names <- function(n_facets) {
  if (n_facets == 1) {
    return("F1")
  }
  return(c("F1", "F2"))
}

# A facet fixed in the G study keeps its observed conditions in the D study,
# n' = N' = n, and may be omitted from universe_sizes and d_sizes.
.fill_fixed_facet_sizes <- function(sizes, argument_name, fixed_facet, n_fixed, fixed_label) {
  values <- unlist(sizes)

  if (!(fixed_facet %in% names(values))) {
    values[fixed_facet] <- n_fixed
    return(values)
  }

  if (!isTRUE(values[[fixed_facet]] == n_fixed)) {
    stop(
      fixed_label, " is fixed in the G study, so the D study uses its ", n_fixed,
      " observed conditions: ", argument_name, "[\"", fixed_facet, "\"] must be ",
      n_fixed, ". To sample from the fixed set, fit the G study without `fixed` ",
      "and give N' in universe_sizes.",
      call. = FALSE
    )
  }

  return(values)
}

# Check one set of D-study sizes; n' cannot exceed a finite N'.
.check_d_study_sizes <- function(universe_values, d_sizes, facet_names) {
  d_values <- .check_size_vector(
    d_sizes, "d_sizes", facet_names, allow_inf = FALSE
  )

  for (facet in facet_names) {
    if (d_values[[facet]] > universe_values[[facet]]) {
      stop(
        "d_sizes[\"", facet, "\"] exceeds universe_sizes[\"", facet, "\"]. ",
        "The D-study sample size cannot exceed a finite universe size.",
        call. = FALSE
      )
    }
  }

  return(d_values[facet_names])
}

.resolve_d_study_design <- function(g_design, d_design, d_nesting, n_facets) {
  if (length(.d_study_facet_names(n_facets)) == 1L) {
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
      "Variance components can be aggregated into a more nested D-study design, ",
      "but already nested components cannot be decomposed into a less nested design.",
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

  # A nested two-facet G-study already fixes the nesting order as F2:F1.
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

  return(list(
    g_design = g_design,
    d_design = d_design,
    d_nesting = d_nesting
  ))
}

.same_design_sources <- function(effect_names) {
  result <- as.list(effect_names)
  names(result) <- effect_names
  return(result)
}

.target_component_sources <- function(g_design, d_design, d_nesting, source_effects) {
  if (g_design == d_design) {
    return(.same_design_sources(source_effects))
  }

  if (g_design == "crossed" && d_design == "nested") {
    result <- list(
      "P",
      c("F1", "PxF1")
    )
    names(result) <- c("P", "F1:P")
    return(result)
  }

  if (g_design == "fully_crossed" && d_design == "partially_nested") {
    if (d_nesting == "F2:F1") {
      result <- list(
        "P",
        "F1",
        "PxF1",
        c("F2", "F1xF2"),
        c("PxF2", "PxF1xF2")
      )
      names(result) <- c("P", "F1", "PxF1", "F2:F1", "F2:F1xP")
      return(result)
    }

    result <- list(
      "P",
      "F2",
      "PxF2",
      c("F1", "F1xF2"),
      c("PxF1", "PxF1xF2")
    )
    names(result) <- c("P", "F2", "PxF2", "F1:F2", "F1:F2xP")
    return(result)
  }

  if (g_design == "fully_crossed" && d_design == "fully_nested") {
    if (d_nesting == "F2:F1") {
      result <- list(
        "P",
        c("F1", "PxF1"),
        c("F2", "PxF2", "F1xF2", "PxF1xF2")
      )
      names(result) <- c("P", "F1:P", "F2:F1:P")
      return(result)
    }

    result <- list(
      "P",
      c("F2", "PxF2"),
      c("F1", "PxF1", "F1xF2", "PxF1xF2")
    )
    names(result) <- c("P", "F2:P", "F1:F2:P")
    return(result)
  }

  if (g_design == "partially_nested" && d_design == "fully_nested") {
    result <- list(
      "P",
      c("F1", "PxF1"),
      c("F2:F1", "F2:F1xP")
    )
    names(result) <- c("P", "F1:P", "F2:F1:P")
    return(result)
  }
}

# A facet fixed in the G study has no fitted main-effect variance. A D-study
# component formed from that main effect alone is its quadratic form and is
# dropped. A component that absorbs the main effect together with other
# effects nests the fixed facet within another facet. That is an error: the
# fixed facet is one shared set of conditions, and nesting it would need a
# separate set in every level of the parent.
.drop_fixed_main_effect <- function(target_sources, fixed_facet, facet_columns) {
  fixed_label <- facet_columns[[fixed_facet]]
  result <- list()

  for (target_effect in names(target_sources)) {
    sources <- target_sources[[target_effect]]

    if (!(fixed_facet %in% sources)) {
      result[[target_effect]] <- sources
      next
    }
    if (length(sources) == 1) {
      next
    }

    # Nesting within P is rejected earlier by .check_fixed_nested_within_p().
    parent <- setdiff(.parse_effect(target_effect), c(fixed_facet, "P"))
    stop(
      fixed_label, " is fixed in the G study, so it is one shared set of conditions ",
      "and cannot be nested within ", facet_columns[[parent]], " in the D study.",
      call. = FALSE
    )
  }

  return(result)
}

# Each facet is random (N' = Inf), finite (n' < N' < Inf), or fixed (n' = N').
.d_facet_status <- function(universe_sizes, d_sizes) {
  universe_values <- unlist(universe_sizes)
  d_values <- unlist(d_sizes)
  facets <- sort(names(d_values))
  result <- vapply(
    facets,
    function(facet) {
      if (is.infinite(universe_values[[facet]])) {
        return("random")
      }
      if (d_values[[facet]] == universe_values[[facet]]) {
        return("fixed")
      }
      return("finite")
    },
    FUN.VALUE = character(1)
  )

  return(result)
}

# A fixed set of conditions common to all persons cannot be nested within P.
.check_fixed_nested_within_p <- function(fixed_facets, design, facet_columns) {
  nested_within_p <- switch(
    design$d_design,
    nested = "F1",
    fully_nested = c("F1", "F2"),
    character(0)
  )
  for (facet in intersect(nested_within_p, fixed_facets)) {
    stop(
      facet_columns[[facet]], " is fixed (n' = N') but nested within ",
      facet_columns[["P"]], " in the D-study design: a fixed set of conditions ",
      "common to all persons cannot be nested within persons.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Check restrictions and approximations for facets fixed in the D study.
.check_fixed_facets <- function(facet_status, d_sizes, g_study_result, design) {
  facet_columns <- g_study_result$facet_columns
  fixed_facets <- names(facet_status)[facet_status == "fixed"]

  if (length(fixed_facets) == length(facet_status)) {
    stop(
      "Every facet has n' = N', so there is no error of measurement. ",
      "At least one facet needs n' < N'.",
      call. = FALSE
    )
  }

  # A facet fixed in the G study was checked by .prepare_d_study().
  .check_fixed_nested_within_p(
    setdiff(fixed_facets, g_study_result$fixed_facet), design, facet_columns
  )

  # Facets fixed only in the D study, whose conditions should be the G-study ones.
  d_only_fixed <- setdiff(fixed_facets, g_study_result$fixed_facet)
  d_values <- unlist(d_sizes)
  for (facet in d_only_fixed) {
    if (!isTRUE(d_values[[facet]] == g_study_result$sample_sizes[[facet]])) {
      warning(
        facet_columns[[facet]], " is fixed at n' = N' = ", d_values[[facet]],
        ", a different number of conditions than the G study observed, so the ",
        "D-study conditions are not the G-study conditions and the results are ",
        "approximations (Brennan 2001, Sec. 4.3).",
        call. = FALSE
      )
    }
  }

  if (design$g_design == "fully_crossed" && design$d_design == "partially_nested") {
    nested_facet <- if (design$d_nesting == "F2:F1") "F2" else "F1"
    parent <- setdiff(c("F1", "F2"), nested_facet)
    if (nested_facet %in% d_only_fixed) {
      warning(
        facet_columns[[nested_facet]], " is fixed and nested within ",
        facet_columns[[parent]], " in the D study but crossed with it in the G study, ",
        "so each level of ", facet_columns[[parent]], " needs its own conditions ",
        "and the results are approximations (Brennan 2001, Sec. 4.3).",
        call. = FALSE
      )
    }
  }
  return(invisible(NULL))
}

.component_transform_matrix <- function(variance_component, target_sources) {
  source_effects <- names(variance_component)
  target_effects <- names(target_sources)

  transform_matrix <- matrix(
    0,
    nrow = length(target_effects),
    ncol = length(source_effects),
    dimnames = list(target_effects, source_effects)
  )

  for (target_effect in target_effects) {
    transform_matrix[target_effect, target_sources[[target_effect]]] <- 1
  }

  return(transform_matrix)
}

.recode_variance_components <- function(variance_component, variance_component_vcov,
                                        g_design, d_design, d_nesting,
                                        fixed_facet = NULL, facet_columns = NULL) {
  target_sources <- .target_component_sources(
    g_design = g_design,
    d_design = d_design,
    d_nesting = d_nesting,
    source_effects = names(variance_component)
  )
  if (!is.null(fixed_facet)) {
    target_sources <- .drop_fixed_main_effect(target_sources, fixed_facet, facet_columns)
  }
  transform_matrix <- .component_transform_matrix(variance_component, target_sources)
  g_effects <- colnames(transform_matrix)

  d_variance_component <- as.numeric(transform_matrix %*% variance_component[g_effects])
  names(d_variance_component) <- rownames(transform_matrix)

  d_variance_component_vcov <- transform_matrix %*%
    variance_component_vcov[g_effects, g_effects, drop = FALSE] %*%
    t(transform_matrix)

  component_mapping <- data.frame(
    effect = rownames(transform_matrix),
    source_effects = vapply(target_sources, paste, character(1), collapse = " + "),
    row.names = NULL
  )

  return(list(
    variance_component = d_variance_component,
    variance_component_vcov = d_variance_component_vcov,
    component_sources = target_sources,
    component_mapping = component_mapping,
    transform_matrix = transform_matrix
  ))
}

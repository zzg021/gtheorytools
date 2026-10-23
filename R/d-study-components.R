# For each component (target), collect the components (sources) that contain
# all of its facets. A target without the object of measurement collects only
# sources without it.
.components_by_target <- function(variance_components, fn,
                                  object_of_measurement = "P") {
  result <- list()

  for (target_name in names(variance_components)) {
    target <- fn(target_name)
    target_has_object <- object_of_measurement %in% target

    component_names <- names(variance_components)[vapply(
      names(variance_components),
      function(source_name) {
        source <- fn(source_name)

        all(target %in% source) &&
          (target_has_object || !(object_of_measurement %in% source))
      },
      logical(1)
    )]

    result[[target_name]] <- variance_components[component_names]
  }

  return(result)
}

# The facets of each component other than the object of measurement.
.child_components_without_object <- function(variance_components,
                                             object_of_measurement = "P",
                                             split_pattern = "[:x]") {
  result <- lapply(names(variance_components), function(component_name) {
    facets <- unlist(strsplit(component_name, split_pattern, perl = TRUE))
    return(facets[facets != "" & facets != object_of_measurement])
  })
  if (length(result) > 0) {
    names(result) <- names(variance_components)
  }

  return(result)
}

# Single-observation components in the universe of generalization (Brennan
# Eq. 5.4): add each contributing higher-order component divided by the
# product of N' over its added facets.
.single_obs_transform_matrix <- function(variance_component, universe_sizes,
                                         object_of_measurement = "P") {
  component_list <- .components_by_target(
    variance_component,
    fn = .parse_effect,
    object_of_measurement
  )
  component_names <- names(component_list)
  transform_matrix <- matrix(
    0,
    nrow = length(component_names),
    ncol = length(variance_component),
    dimnames = list(component_names, names(variance_component))
  )

  for (target_name in component_names) {
    target <- .parse_effect(target_name)
    for (source_name in names(component_list[[target_name]])) {
      denominator <- 1
      product_facets <- setdiff(.parse_effect(source_name), target)

      for (facet in product_facets) {
        denominator <- denominator * universe_sizes[[facet]]
      }

      transform_matrix[target_name, source_name] <- 1 / denominator
    }
  }

  return(transform_matrix)
}

.single_obs_vc <- function(variance_component, universe_sizes, object_of_measurement = "P") {
  transform_matrix <- .single_obs_transform_matrix(
    variance_component,
    universe_sizes,
    object_of_measurement
  )
  result_single_obs <- as.numeric(
    transform_matrix %*% variance_component[colnames(transform_matrix)]
  )
  return(matrix(
    result_single_obs, 1, dimnames = list(NULL, rownames(transform_matrix))
  ))
}

.single_obs_vcov <- function(variance_component, universe_sizes, variance_component_vcov,
                             object_of_measurement = "P") {
  transform_matrix <- .single_obs_transform_matrix(
    variance_component,
    universe_sizes,
    object_of_measurement
  )
  component_names <- rownames(transform_matrix)

  result <- transform_matrix %*%
    variance_component_vcov[names(variance_component), names(variance_component), drop = FALSE] %*%
    t(transform_matrix)

  dimnames(result) <- list(component_names, component_names)
  return(result)
}

# Multiplier C_A / h_A from the single-observation to the mean-score component
# (Brennan Eqs. 5.2-5.4). h_A is the product of n' over the facets of A other
# than P. C_A is the product of (1 - n'/N') over the facets of A, other than P,
# that precede any ":": every facet of a crossed component, and only the
# innermost facet of a nested one.
.d_mean_multiplier <- function(variance_component, single_obs, universe_sizes,
                               d_sizes, object_of_measurement = "P") {
  component_children <- .child_components_without_object(
    variance_component,
    object_of_measurement
  )
  result_multiplier <- numeric(length(single_obs))
  names(result_multiplier) <- colnames(single_obs)

  for (i in seq_along(single_obs)) {
    multiplier <- 1
    component_name <- colnames(single_obs)[i]
    nested_portion <- strsplit(component_name, "x", fixed = TRUE)[[1]][1]
    if (grepl(":", nested_portion, fixed = TRUE)) {
      nested_facets <- unlist(strsplit(nested_portion, ":", fixed = TRUE))
      nested_facets <- nested_facets[
        nested_facets != "" & nested_facets != object_of_measurement
      ]
      innermost_facet <- nested_facets[1]
      outer_sizes <- vapply(
        nested_facets[-1],
        function(facet) d_sizes[[facet]],
        numeric(1)
      )
      multiplier <- (1 - d_sizes[[innermost_facet]] / universe_sizes[[innermost_facet]]) /
        d_sizes[[innermost_facet]] / prod(outer_sizes)
    } else {
      for (facet in component_children[[component_name]]) {
        multiplier <- multiplier *
          (1 - d_sizes[[facet]] / universe_sizes[[facet]]) / d_sizes[[facet]]
      }
    }
    result_multiplier[i] <- multiplier
  }

  return(result_multiplier)
}

.d_mean_vc <- function(variance_component, single_obs, universe_sizes, d_sizes,
                       object_of_measurement = "P") {
  multiplier <- .d_mean_multiplier(
    variance_component, single_obs, universe_sizes, d_sizes, object_of_measurement
  )
  result_mean_score <- as.numeric(single_obs) * multiplier
  return(matrix(
    result_mean_score, 1, dimnames = list(NULL, colnames(single_obs))
  ))
}

.d_mean_vcov <- function(variance_component, single_obs, single_obs_vcov,
                         universe_sizes, d_sizes, object_of_measurement = "P") {
  multiplier <- .d_mean_multiplier(
    variance_component, single_obs, universe_sizes, d_sizes, object_of_measurement
  )
  component_names <- colnames(single_obs)

  transform_matrix <- diag(multiplier)
  dimnames(transform_matrix) <- list(component_names, component_names)

  result <- transform_matrix %*%
    single_obs_vcov[component_names, component_names, drop = FALSE] %*%
    t(transform_matrix)

  dimnames(result) <- list(component_names, component_names)
  return(result)
}

.single_obs_vc_table <- function(single_obs, single_obs_vcov = NULL) {
  single_obs_values <- as.numeric(single_obs)
  names(single_obs_values) <- colnames(single_obs)
  standard_error <- rep(NA_real_, length(single_obs_values))

  if (!is.null(single_obs_vcov)) {
    standard_error <- sqrt(diag(single_obs_vcov)[names(single_obs_values)])
  }

  return(data.frame(
    effect = names(single_obs_values),
    single_obs_estimate = unname(single_obs_values),
    standard_error = unname(standard_error),
    row.names = NULL
  ))
}

.d_mean_vc_table <- function(d_mean, d_mean_vcov = NULL) {
  d_mean_values <- as.numeric(d_mean)
  names(d_mean_values) <- colnames(d_mean)
  standard_error <- rep(NA_real_, length(d_mean_values))

  if (!is.null(d_mean_vcov)) {
    standard_error <- sqrt(diag(d_mean_vcov)[names(d_mean_values)])
  }

  return(data.frame(
    effect = names(d_mean_values),
    mean_score_estimate = unname(d_mean_values),
    standard_error = unname(standard_error),
    row.names = NULL
  ))
}

.g_coefficient <- function(d_mean, d_mean_vcov = NULL, object_of_measurement = "P") {
  d_mean_values <- as.numeric(d_mean)
  names(d_mean_values) <- colnames(d_mean)

  has_object <- vapply(
    names(d_mean_values),
    function(component_name) object_of_measurement %in% .parse_effect(component_name),
    logical(1)
  )
  names(has_object) <- names(d_mean_values)

  tau_weight <- numeric(length(d_mean_values))
  delta_weight <- numeric(length(d_mean_values))
  Delta_weight <- numeric(length(d_mean_values))
  names(tau_weight) <- names(d_mean_values)
  names(delta_weight) <- names(d_mean_values)
  names(Delta_weight) <- names(d_mean_values)

  tau_weight[object_of_measurement] <- 1
  delta_weight[has_object] <- 1
  delta_weight[object_of_measurement] <- 0
  Delta_weight <- delta_weight
  Delta_weight[!has_object] <- 1

  tau <- sum(tau_weight * d_mean_values)
  delta <- sum(delta_weight * d_mean_values)
  Delta <- sum(Delta_weight * d_mean_values)
  G <- tau / (tau + delta)
  Phi <- tau / (tau + Delta)

  standard_error <- rep(NA_real_, 5)

  if (!is.null(d_mean_vcov)) {
    component_names <- names(d_mean_values)
    d_mean_vcov <- d_mean_vcov[component_names, component_names, drop = FALSE]

    G_weight <- (delta / (tau + delta)^2) * tau_weight -
      (tau / (tau + delta)^2) * delta_weight
    Phi_weight <- (Delta / (tau + Delta)^2) * tau_weight -
      (tau / (tau + Delta)^2) * Delta_weight

    se_weight <- rbind(
      tau = tau_weight,
      delta = delta_weight,
      Delta = Delta_weight,
      G = G_weight,
      Phi = Phi_weight
    )

    standard_error <- sqrt(diag(se_weight %*% d_mean_vcov %*% t(se_weight)))
  }

  return(data.frame(
    quantity = c("tau", "delta", "Delta", "G", "Phi"),
    estimate = c(tau, delta, Delta, G, Phi),
    standard_error = unname(standard_error),
    row.names = NULL
  ))
}

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

# Multipliers for a D study with the G-study pattern of unbalanced data
# (Brennan 2001, Sec. 7.2), every facet random. Person p's mean score gives the
# alpha cell c the weight a_pc = n_pc / n_p+, its share of the person's scores.
#
# The absolute-error multiplier of sigma^2(alpha) is the mean over persons of
# sum_c a_pc^2: 1 over the harmonic mean of the n_i:p in i:p (Eqs. 7.18-7.19),
# 1 / n~_h for h and ph and 1 / n_i+ for i:h and pi:h in p x (I:H) (Eqs. 7.27
# and 7.28), and 1 over the harmonic mean of the n_p in p x I with missing
# scores (Eqs. 7.41-7.42).
#
# The relative-error multiplier makes sigma^2(delta) the expected variance of
# the person means less sigma^2(p), the denominator of E rho^2 in Eq. 7.21:
# sum_p sum_c (a_pc - abar_c)^2 / (n_p - 1), with abar_c the mean of a_pc over
# all persons. Persons share no cell of an effect that contains p, so for such
# an effect it equals the absolute multiplier. For an effect without p it is
# the absolute multiplier less the mean over pairs of persons of
# sum_c a_pc a_p'c (Eq. 7.45), and it is zero when every person has the same
# number of scores in each cell, as in Eq. 7.29.
.pattern_multipliers <- function(observed_design, effect_names) {
  n_persons <- nlevels(observed_design$P)
  person <- as.integer(observed_design$P)
  scores_per_person <- tabulate(person, nbins = n_persons)
  absolute <- numeric(length(effect_names))
  relative <- numeric(length(effect_names))
  names(absolute) <- effect_names
  names(relative) <- effect_names

  for (effect in effect_names) {
    facets <- .parse_effect(effect)
    person_cells <- .cell_index(observed_design[union("P", facets)])
    n_person_cells <- max(person_cells)
    counts <- tabulate(person_cells, nbins = n_person_cells)
    first_score <- match(seq_len(n_person_cells), person_cells)
    share <- counts / scores_per_person[person[first_score]]
    absolute[[effect]] <- sum(share^2) / n_persons

    if ("P" %in% facets) {
      relative[[effect]] <- absolute[[effect]]
      next
    }
    cells <- .cell_index(observed_design[facets])
    cell_of_share <- cells[first_score]
    persons_in_cell <- tabulate(cell_of_share, nbins = max(cells))
    first_in_cell <- match(cell_of_share, cell_of_share)
    if (all(persons_in_cell == n_persons) && all(counts == counts[first_in_cell])) {
      relative[[effect]] <- 0
      next
    }
    mean_share <- as.vector(rowsum(share, cell_of_share)) / n_persons
    deviation <- share - mean_share[cell_of_share]
    relative[[effect]] <- (sum(deviation^2) +
      sum((n_persons - persons_in_cell) * mean_share^2)) / (n_persons - 1)
  }

  return(list(absolute = absolute, relative = relative))
}

.d_mean_vc <- function(single_obs, multiplier) {
  result_mean_score <- as.numeric(single_obs) * multiplier[colnames(single_obs)]
  return(matrix(
    result_mean_score, 1, dimnames = list(NULL, colnames(single_obs))
  ))
}

.d_mean_vcov <- function(single_obs, single_obs_vcov, multiplier) {
  component_names <- colnames(single_obs)

  transform_matrix <- diag(multiplier[component_names], nrow = length(component_names))
  dimnames(transform_matrix) <- list(component_names, component_names)

  result <- transform_matrix %*%
    single_obs_vcov[component_names, component_names, drop = FALSE] %*%
    t(transform_matrix)

  dimnames(result) <- list(component_names, component_names)
  return(result)
}

# Standard errors from estimated sampling variances. With unbalanced ANOVA
# estimates, the estimated covariance of the G-study components can have a
# negative eigenvalue, so a linear combination of them can have an estimated
# variance below zero; its standard error is NA, as in the G study.
.standard_errors <- function(variances) {
  negative <- !is.na(variances) & variances < 0
  if (any(negative)) {
    warning(
      "The estimated sampling variance of ", paste(names(variances)[negative], collapse = ", "),
      " is below zero; its standard error is NA.",
      call. = FALSE
    )
  }
  result <- sqrt(ifelse(negative, NA_real_, variances))
  names(result) <- names(variances)
  return(result)
}

.single_obs_vc_table <- function(single_obs, single_obs_vcov = NULL) {
  single_obs_values <- as.numeric(single_obs)
  names(single_obs_values) <- colnames(single_obs)
  standard_error <- rep(NA_real_, length(single_obs_values))

  if (!is.null(single_obs_vcov)) {
    standard_error <- .standard_errors(diag(single_obs_vcov)[names(single_obs_values)])
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
    standard_error <- .standard_errors(diag(d_mean_vcov)[names(d_mean_values)])
  }

  return(data.frame(
    effect = names(d_mean_values),
    mean_score_estimate = unname(d_mean_values),
    standard_error = unname(standard_error),
    row.names = NULL
  ))
}

# tau, delta, Delta, G, and Phi from the mean-score components. delta_weight
# gives the share of each component in relative error; by default it is 1 for
# the components with the object of measurement, other than its own, and 0
# for the others.
.g_coefficient <- function(d_mean, d_mean_vcov = NULL, object_of_measurement = "P",
                           delta_weight = NULL) {
  d_mean_values <- as.numeric(d_mean)
  names(d_mean_values) <- colnames(d_mean)

  has_object <- vapply(
    names(d_mean_values),
    function(component_name) object_of_measurement %in% .parse_effect(component_name),
    logical(1)
  )
  names(has_object) <- names(d_mean_values)

  tau_weight <- numeric(length(d_mean_values))
  Delta_weight <- rep(1, length(d_mean_values))
  names(tau_weight) <- names(d_mean_values)
  names(Delta_weight) <- names(d_mean_values)
  tau_weight[object_of_measurement] <- 1
  Delta_weight[object_of_measurement] <- 0

  if (is.null(delta_weight)) {
    delta_weight <- as.numeric(has_object)
    names(delta_weight) <- names(d_mean_values)
    delta_weight[object_of_measurement] <- 0
  }
  delta_weight <- delta_weight[names(d_mean_values)]

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

    standard_error <- .standard_errors(diag(se_weight %*% d_mean_vcov %*% t(se_weight)))
  }

  return(data.frame(
    quantity = c("tau", "delta", "Delta", "G", "Phi"),
    estimate = c(tau, delta, Delta, G, Phi),
    standard_error = unname(standard_error),
    row.names = NULL
  ))
}

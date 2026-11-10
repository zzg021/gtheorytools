# Single-observation components in the universe of generalization (Brennan
# Eq. 5.4) as weights on the D-study components: each component adds the
# components that contain all of its facets, each divided by the product of N'
# over the facets it adds. A component without the object of measurement adds
# only components without it.
.single_obs_weights <- function(effects, universe_sizes, object_of_measurement) {
  result <- matrix(
    0,
    nrow = length(effects),
    ncol = length(effects),
    dimnames = list(effects, effects)
  )

  for (target in effects) {
    target_facets <- .parse_effect(target)
    target_has_object <- object_of_measurement %in% target_facets
    for (source in effects) {
      source_facets <- .parse_effect(source)
      if (!all(target_facets %in% source_facets)) {
        next
      }
      if (!target_has_object && object_of_measurement %in% source_facets) {
        next
      }
      added_facets <- setdiff(source_facets, target_facets)
      result[target, source] <- 1 / prod(universe_sizes[added_facets])
    }
  }

  return(result)
}

# Multiplier C_A / h_A from the single-observation to the mean-score component
# (Brennan Eqs. 5.2-5.4). h_A is the product of n' over the facets of A other
# than the object of measurement, and C_A is the product of (1 - n'/N') over
# the primary indices of A other than the object of measurement.
.d_mean_multiplier <- function(single_obs, universe_sizes, d_sizes, d_parents,
                               object_of_measurement) {
  component_names <- colnames(single_obs)
  result_multiplier <- numeric(length(component_names))
  names(result_multiplier) <- component_names

  for (component_name in component_names) {
    facets <- .parse_effect(component_name)
    primary_facets <- .primary_indices(facets, d_parents)
    multiplier <- 1

    for (facet in setdiff(primary_facets, object_of_measurement)) {
      multiplier <- multiplier * (1 - d_sizes[[facet]] / universe_sizes[[facet]])
    }
    for (facet in setdiff(facets, object_of_measurement)) {
      multiplier <- multiplier / d_sizes[[facet]]
    }
    result_multiplier[[component_name]] <- multiplier
  }

  return(result_multiplier)
}

# Multipliers for a D study with the G-study pattern of unbalanced data
# (Brennan 2001, Sec. 7.2), every facet random. With o the object of
# measurement, the mean score of level o gives the alpha cell c the weight
# a_oc = n_oc / n_o+, its share of the scores of o.
#
# The absolute-error multiplier of sigma^2(alpha) is the mean over the levels
# of o of sum_c a_oc^2: with persons as o, 1 over the harmonic mean of the
# n_i:p in i:p (Eqs. 7.18-7.19), 1 / n~_h for h and ph and 1 / n_i+ for i:h
# and pi:h in p x (i:h) (Eqs. 7.27 and 7.28), and 1 over the harmonic mean of
# the n_p in p x i with missing scores (Eqs. 7.41-7.42).
#
# The relative-error multiplier makes sigma^2(delta) the expected variance of
# the means of the levels of o less sigma^2(o), the denominator of E rho^2 in
# Eq. 7.21: sum_o sum_c (a_oc - abar_c)^2 / (n_o - 1), with abar_c the mean of
# a_oc over all levels of o. Levels of o share no cell of an effect that
# contains o, so for such an effect it equals the absolute multiplier. For an
# effect without o it is the absolute multiplier less the mean over pairs of
# levels of sum_c a_oc a_o'c (Eq. 7.45), and it is zero when every level has
# the same number of scores in each cell, as in Eq. 7.29.
.pattern_multipliers <- function(observed_design, effect_names, object_of_measurement) {
  object_levels <- observed_design[[object_of_measurement]]
  n_objects <- nlevels(object_levels)
  object <- as.integer(object_levels)
  scores_per_object <- tabulate(object, nbins = n_objects)
  absolute <- numeric(length(effect_names))
  relative <- numeric(length(effect_names))
  names(absolute) <- effect_names
  names(relative) <- effect_names

  for (effect in effect_names) {
    facets <- .parse_effect(effect)
    object_cells <- .cell_index(observed_design[union(object_of_measurement, facets)])
    n_object_cells <- max(object_cells)
    counts <- tabulate(object_cells, nbins = n_object_cells)
    first_score <- match(seq_len(n_object_cells), object_cells)
    share <- counts / scores_per_object[object[first_score]]
    absolute[[effect]] <- sum(share^2) / n_objects

    if (object_of_measurement %in% facets) {
      relative[[effect]] <- absolute[[effect]]
      next
    }
    cells <- .cell_index(observed_design[facets])
    cell_of_share <- cells[first_score]
    objects_in_cell <- tabulate(cell_of_share, nbins = max(cells))
    first_in_cell <- match(cell_of_share, cell_of_share)
    if (all(objects_in_cell == n_objects) && all(counts == counts[first_in_cell])) {
      relative[[effect]] <- 0
      next
    }
    mean_share <- as.vector(rowsum(share, cell_of_share)) / n_objects
    deviation <- share - mean_share[cell_of_share]
    relative[[effect]] <- (sum(deviation^2) +
      sum((n_objects - objects_in_cell) * mean_share^2)) / (n_objects - 1)
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
  weights <- diag(multiplier[component_names], nrow = length(component_names))
  dimnames(weights) <- list(component_names, component_names)

  result <- weights %*%
    single_obs_vcov[component_names, component_names, drop = FALSE] %*%
    weights
  return(result)
}

# Standard errors from estimated sampling variances. With unbalanced ANOVA
# estimates, the estimated covariance of the G-study components can have a
# negative eigenvalue, so a linear combination of them can have an estimated
# variance below zero; its standard error is NA, as in the G study.
.standard_errors <- function(variances, labels = names(variances)) {
  negative <- !is.na(variances) & variances < 0
  if (any(negative)) {
    warning(
      "The estimated sampling variance of ", paste(labels[negative], collapse = ", "),
      " is below zero; its standard error is NA.",
      call. = FALSE
    )
  }
  result <- sqrt(ifelse(negative, NA_real_, variances))
  names(result) <- names(variances)
  return(result)
}

# A table of D-study components from their estimates, a one-row matrix, and
# their covariance matrix, which is NULL when standard errors are unavailable.
# Each effect has its label and the labels of the G-study components summed
# into it.
.component_table <- function(estimates, vcov, estimate_column, effect_labels,
                             source_effects) {
  effects <- colnames(estimates)
  standard_error <- rep(NA_real_, length(effects))
  if (!is.null(vcov)) {
    standard_error <- .standard_errors(diag(vcov)[effects], effect_labels[effects])
  }

  result <- data.frame(
    effect = effects,
    effect_label = unname(effect_labels[effects]),
    source_effects = unname(source_effects[effects]),
    row.names = NULL
  )
  result[[estimate_column]] <- as.numeric(estimates)
  result$standard_error <- unname(standard_error)
  return(result)
}

# tau, delta, Delta, G, and Phi from the mean-score components. delta_weight
# gives the share of each component in relative error; by default it is 1 for
# the components with the object of measurement, other than its own, and 0
# for the others.
.g_coefficient <- function(d_mean, d_mean_vcov, object_of_measurement,
                           delta_weight = NULL) {
  component_names <- colnames(d_mean)
  d_mean_values <- as.numeric(d_mean)
  tau_weight <- as.numeric(component_names == object_of_measurement)
  Delta_weight <- 1 - tau_weight
  if (is.null(delta_weight)) {
    has_object <- vapply(
      component_names,
      function(component_name) object_of_measurement %in% .parse_effect(component_name),
      logical(1)
    )
    delta_weight <- as.numeric(has_object) - tau_weight
  } else {
    delta_weight <- unname(delta_weight[component_names])
  }

  tau <- sum(tau_weight * d_mean_values)
  delta <- sum(delta_weight * d_mean_values)
  Delta <- sum(Delta_weight * d_mean_values)
  G <- tau / (tau + delta)
  Phi <- tau / (tau + Delta)

  standard_error <- rep(NA_real_, 5)

  if (!is.null(d_mean_vcov)) {
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

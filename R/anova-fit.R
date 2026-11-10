# For each effect, the other effects whose facets are all among its facets.
.lower_effects <- function(facet_sets) {
  result <- lapply(names(facet_sets), function(effect) {
    is_lower <- vapply(
      names(facet_sets),
      function(other) other != effect && all(facet_sets[[other]] %in% facet_sets[[effect]]),
      FUN.VALUE = logical(1)
    )
    return(names(facet_sets)[is_lower])
  })
  names(result) <- names(facet_sets)
  return(result)
}

# Degrees of freedom, T terms, sums of squares, and mean squares of each effect
# (Brennan 2001, Sec. 7.1.1). T(alpha) sums, over the alpha cells, the squared
# total of the scores in the cell divided by their number (Eq. 7.1). SS(alpha)
# is T(alpha) less T(mu) and the sums of squares of the effects whose facets
# are among those of alpha, and df(alpha) is the number of alpha cells less one
# and the df of those effects. For balanced data these are the usual sums of
# squares (Eq. 3.11); otherwise they are the analogous sums of squares of
# Henderson's Method 1. The T terms are computed from scores centered at the
# grand mean, which lowers each by n+ times the squared grand mean and leaves
# the sums of squares unchanged.
.anova_table <- function(dat, effect_names) {
  facet_sets <- lapply(effect_names, .parse_effect)
  names(facet_sets) <- effect_names
  lower_effects <- .lower_effects(facet_sets)
  grand_mean <- mean(dat$score)
  centered_score <- dat$score - grand_mean
  t_mean <- nrow(dat) * grand_mean^2

  df <- numeric(0)
  t_term <- numeric(0)
  sum_of_squares <- numeric(0)
  for (effect in effect_names[order(lengths(facet_sets))]) {
    cells <- .cell_index(dat[facet_sets[[effect]]])
    n_cells <- max(cells)
    cell_totals <- as.vector(rowsum(centered_score, cells))
    cell_counts <- tabulate(cells, nbins = n_cells)
    centered_t <- sum(cell_totals^2 / cell_counts)
    lower <- lower_effects[[effect]]
    df[[effect]] <- n_cells - 1 - sum(df[lower])
    t_term[[effect]] <- centered_t + t_mean
    sum_of_squares[[effect]] <- centered_t - sum(sum_of_squares[lower])
  }

  result <- data.frame(
    effect = effect_names,
    df = unname(df[effect_names]),
    t_term = unname(t_term[effect_names]),
    sum_of_squares = unname(sum_of_squares[effect_names]),
    mean_square = unname(sum_of_squares[effect_names] / df[effect_names]),
    row.names = NULL
  )
  return(result)
}

# Every effect needs positive degrees of freedom. With missing scores an
# interaction can have none, as when persons share no conditions.
.check_anova_df <- function(anova_table, effect_labels) {
  no_df <- anova_table$df < 1
  if (any(no_df)) {
    stop(
      "method = \"anova\" requires positive degrees of freedom for every effect. ",
      "In these data, ", paste(effect_labels[no_df], collapse = ", "),
      if (sum(no_df) == 1) " has " else " have ", "none.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Coefficient of sigma^2(beta) in the expected T term of alpha (Brennan 2001,
# Eq. 7.3): for each alpha cell, the squared numbers of scores in the cells of
# the facets of alpha and beta together, summed and divided by the number of
# scores in the alpha cell, then summed over the alpha cells. alpha = mu has a
# single cell holding all scores.
.expected_t_coefficient <- function(dat, alpha_facets, beta_facets) {
  joint_cells <- .cell_index(dat[union(alpha_facets, beta_facets)])
  n_joint_cells <- max(joint_cells)
  joint_counts <- tabulate(joint_cells, nbins = n_joint_cells)
  if (length(alpha_facets) == 0) {
    return(sum(joint_counts^2) / nrow(dat))
  }

  alpha_cells <- .cell_index(dat[alpha_facets])
  alpha_counts <- tabulate(alpha_cells, nbins = max(alpha_cells))
  first_score <- match(seq_len(n_joint_cells), joint_cells)
  squares_by_alpha_cell <- as.vector(rowsum(joint_counts^2, alpha_cells[first_score]))
  return(sum(squares_by_alpha_cell / alpha_counts))
}

# Expected mean squares: row alpha and column beta give the coefficient of
# sigma^2(beta) in E[MS(alpha)]. The expected sums of squares follow from the
# expected T terms as the sums of squares follow from the T terms, and the
# coefficient of mu^2, n+ in every T term, cancels. For balanced data the
# coefficient is the number of scores in each beta cell when beta contains the
# facets of alpha, and zero otherwise (Eq. 3.24).
.ems_coefficients <- function(dat, effect_names, df) {
  facet_sets <- lapply(effect_names, .parse_effect)
  names(facet_sets) <- effect_names
  lower_effects <- .lower_effects(facet_sets)

  expected_t <- function(alpha_facets) {
    result <- vapply(
      effect_names,
      function(beta) .expected_t_coefficient(dat, alpha_facets, facet_sets[[beta]]),
      FUN.VALUE = numeric(1)
    )
    return(result)
  }
  expected_t_mean <- expected_t(character(0))

  expected_ss <- matrix(
    0,
    nrow = length(effect_names),
    ncol = length(effect_names),
    dimnames = list(effect_names, effect_names)
  )
  for (effect in effect_names[order(lengths(facet_sets))]) {
    lower <- lower_effects[[effect]]
    expected_ss[effect, ] <- expected_t(facet_sets[[effect]]) - expected_t_mean -
      colSums(expected_ss[lower, , drop = FALSE])
  }

  result <- expected_ss / df
  return(result)
}

# Random-effects variance components solve E[MS] = MS (Brennan 2001, Eq. 3.25),
# so each estimate is a linear combination of mean squares. For balanced data
# the mean squares are independent under normality, and the covariance of the
# estimates uses 2 MS^2 / (df + 2) as the variance of each mean square
# (Eq. 6.1). For unbalanced data the analogous mean squares are correlated, and
# the covariance is estimated from all of them together. An estimated
# sampling variance below zero gives that component an NA standard error.
.anova_components <- function(anova_table, ems_coefficients, dat, balanced, effect_labels) {
  effects <- anova_table$effect
  mean_squares <- anova_table$mean_square
  ms_weights <- solve(ems_coefficients)
  vc <- as.vector(ms_weights %*% mean_squares)
  names(vc) <- effects

  if (balanced) {
    ms_variance <- 2 * mean_squares^2 / (anova_table$df + 2)
    vcov <- ms_weights %*% diag(ms_variance, nrow = length(ms_variance)) %*% t(ms_weights)
  } else {
    vcov <- .anova_covariance_unbalanced(
      dat, effects, anova_table$df, ems_coefficients, vc
    )
    negative <- diag(vcov) < 0
    if (any(negative)) {
      warning(
        "The estimated sampling variance of ",
        paste(effect_labels[negative], collapse = ", "),
        " is below zero; its standard error is NA.",
        call. = FALSE
      )
      vcov[negative, ] <- NA_real_
      vcov[, negative] <- NA_real_
    }
  }
  dimnames(vcov) <- list(effects, effects)
  se <- sqrt(diag(vcov))
  names(se) <- effects

  return(list(vc = vc, se = se, vcov = vcov))
}

# The expected mean squares written out, as in GENOVA's output, with the
# coefficients rounded to three decimals.
.print_expected_mean_squares <- function(ems_coefficients, effect_labels) {
  effects <- rownames(ems_coefficients)
  left_sides <- format(paste0("E[MS(", effect_labels[effects], ")]"))
  names(left_sides) <- effects

  cat("Expected mean squares:\n")
  for (effect in effects) {
    coefficients <- round(ems_coefficients[effect, ], 3)
    sources <- rev(effects[coefficients != 0])
    terms <- character(length(sources))
    for (k in seq_along(sources)) {
      coefficient <- coefficients[[sources[[k]]]]
      sign <- if (coefficient < 0) "- " else if (k == 1) "" else "+ "
      size <- if (abs(coefficient) == 1) "" else paste0(format(abs(coefficient)), " ")
      terms[[k]] <- paste0(sign, size, "sigma^2(", effect_labels[[sources[[k]]]], ")")
    }
    cat("  ", left_sides[[effect]], " = ", paste(terms, collapse = " "), "\n", sep = "")
  }
  return(invisible(NULL))
}

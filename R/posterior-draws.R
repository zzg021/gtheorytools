# `summary`, such as median(), of each column of a matrix of posterior draws,
# named by the columns.
.column_summary <- function(draws, summary) {
  result <- vapply(colnames(draws), function(column) summary(draws[, column]), numeric(1))
  return(result)
}

# Posterior median, standard deviation, and equal-tailed interval of each
# column of a matrix of posterior draws. A quantity that is undefined in a
# draw, as G is when tau and delta are both zero, is NA in all four.
.draw_summary <- function(draws, level) {
  tail_probability <- (1 - level) / 2
  is_defined <- colSums(is.na(draws)) == 0
  defined_draws <- draws[, is_defined, drop = FALSE]
  lower <- rep(NA_real_, ncol(draws))
  upper <- rep(NA_real_, ncol(draws))
  names(lower) <- colnames(draws)
  names(upper) <- colnames(draws)
  lower[is_defined] <- .column_summary(defined_draws, function(x) {
    return(quantile(x, tail_probability, names = FALSE))
  })
  upper[is_defined] <- .column_summary(defined_draws, function(x) {
    return(quantile(x, 1 - tail_probability, names = FALSE))
  })

  return(list(
    estimate = .column_summary(draws, median),
    se = .column_summary(draws, sd),
    lower = lower,
    upper = upper
  ))
}

# Add the bounds of the credible intervals to a result table. The column
# names give the level in percent: lower_95 and upper_95 for a level of 0.95.
.add_interval_columns <- function(table, lower, upper, level) {
  percent <- format(100 * level)
  table[[paste0("lower_", percent)]] <- unname(lower)
  table[[paste0("upper_", percent)]] <- unname(upper)
  return(table)
}

# Variance components from their posterior draws, as the REML and ANOVA
# estimates are returned: posterior medians, posterior standard deviations,
# and the covariance matrix of the draws.
.posterior_components <- function(draws) {
  vc <- .column_summary(draws, median)
  se <- .column_summary(draws, sd)
  return(list(vc = vc, se = se, vcov = cov(draws)))
}

# The posterior medians as a one-row matrix, the form of the single-observation
# and mean-score components in a D study.
.posterior_median_row <- function(draws) {
  medians <- .column_summary(draws, median)
  return(matrix(medians, nrow = 1, dimnames = list(NULL, colnames(draws))))
}

# The note under a table of Bayesian estimates.
.print_posterior_note <- function(level) {
  cat("Estimates are posterior medians and standard errors posterior standard\n")
  cat("deviations, with ", format(100 * level), "% equal-tailed credible intervals.\n", sep = "")
  return(invisible(NULL))
}

# The note that names the random effects with fewer than 5 levels, by the
# labels of the variance-component table. Their data say little about their
# variance, so the prior shows in the estimate.
.print_few_levels_note <- function(n_levels, variance_components) {
  has_few <- n_levels < 5
  if (!any(has_few)) {
    return(invisible(NULL))
  }
  effect_rows <- match(names(n_levels)[has_few], variance_components$effect)
  effect_labels <- variance_components$effect_label[effect_rows]
  few_labels <- paste0(effect_labels, " (", n_levels[has_few], ")")
  cat(
    "Effects with fewer than 5 levels, whose estimates depend on the priors: ",
    paste(few_labels, collapse = ", "), "\n",
    sep = ""
  )
  return(invisible(NULL))
}

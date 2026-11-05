# ANOVA estimates use the expected mean squares of a balanced design, which
# hold only with one non-missing score in every cell and the same number of
# nested conditions in each level of their parent.
.check_anova_data <- function(dat_for_lmer, sample_sizes) {
  if (!.is_balanced(dat_for_lmer, sample_sizes)) {
    stop(
      "method = \"anova\" requires complete, balanced data: one non-missing ",
      "score for every person in every condition, with the same number of ",
      "nested conditions in each level of their parent. Use method = \"reml\" ",
      "for other data.",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# Degrees of freedom, sums of squares, and mean squares of each effect. The
# observed score effect of alpha is the mean of each alpha cell less the grand
# mean and the score effects of the components whose facets are among those of
# alpha (Brennan 2001, Sec. 3.3.1), and SS(alpha) sums its square over all
# scores (Eq. 3.11). df(alpha) is the number of alpha cells less one and the df
# of those components.
.anova_table <- function(dat, effect_names) {
  facet_sets <- lapply(effect_names, .parse_effect)
  names(facet_sets) <- effect_names
  grand_mean <- mean(dat$score)

  score_effects <- list()
  n_cells <- numeric(0)
  df <- numeric(0)
  sum_of_squares <- numeric(0)
  for (effect in effect_names[order(lengths(facet_sets))]) {
    cells <- interaction(dat[facet_sets[[effect]]], drop = TRUE)
    lower_effects <- character(0)
    for (lower in names(score_effects)) {
      if (all(facet_sets[[lower]] %in% facet_sets[[effect]])) {
        lower_effects <- c(lower_effects, lower)
      }
    }

    lower_fit <- grand_mean + Reduce(`+`, score_effects[lower_effects], 0)
    score_effects[[effect]] <- ave(dat$score, cells) - lower_fit
    n_cells[[effect]] <- nlevels(cells)
    df[[effect]] <- nlevels(cells) - 1 - sum(df[lower_effects])
    sum_of_squares[[effect]] <- sum(score_effects[[effect]]^2)
  }

  table <- data.frame(
    effect = effect_names,
    df = unname(df[effect_names]),
    sum_of_squares = unname(sum_of_squares[effect_names]),
    mean_square = unname(sum_of_squares[effect_names] / df[effect_names]),
    row.names = NULL
  )
  return(list(table = table, n_cells = n_cells[effect_names]))
}

# Random-effects variance components from the mean squares. sigma^2(beta)
# enters E[MS(alpha)] when beta contains the facets of alpha, with the number
# of scores in each beta cell as its coefficient (Brennan 2001, Eq. 3.24), and
# the estimates solve E[MS] = MS (Eq. 3.25). Each estimate is a linear
# combination of mean squares, which are independent under normality, so the
# covariance of the estimates uses 2 MS^2 / (df + 2) as the variance of each
# mean square (Eq. 6.1).
.anova_components <- function(anova, n_scores) {
  effects <- anova$table$effect
  facet_sets <- lapply(effects, .parse_effect)
  ems_coefficients <- matrix(
    0,
    nrow = length(effects),
    ncol = length(effects),
    dimnames = list(effects, effects)
  )
  for (a in seq_along(effects)) {
    for (b in seq_along(effects)) {
      if (all(facet_sets[[a]] %in% facet_sets[[b]])) {
        ems_coefficients[a, b] <- n_scores / anova$n_cells[[b]]
      }
    }
  }

  mean_squares <- anova$table$mean_square
  ms_weights <- solve(ems_coefficients)
  vc <- as.vector(ms_weights %*% mean_squares)
  names(vc) <- effects

  ms_variance <- 2 * mean_squares^2 / (anova$table$df + 2)
  vcov <- ms_weights %*% diag(ms_variance, nrow = length(ms_variance)) %*% t(ms_weights)
  dimnames(vcov) <- list(effects, effects)
  se <- sqrt(diag(vcov))

  return(list(vc = vc, se = se, vcov = vcov, ems_coefficients = ems_coefficients))
}

# The expected mean squares written out, as in GENOVA's output.
.print_expected_mean_squares <- function(ems_coefficients, effect_labels) {
  effects <- rownames(ems_coefficients)
  left_sides <- format(paste0("E[MS(", effect_labels[effects], ")]"))
  names(left_sides) <- effects

  cat("Expected mean squares:\n")
  for (effect in effects) {
    sources <- rev(effects[ems_coefficients[effect, ] != 0])
    terms <- vapply(
      sources,
      function(source) {
        coefficient <- ems_coefficients[effect, source]
        prefix <- if (coefficient == 1) "" else paste0(format(coefficient), " ")
        return(paste0(prefix, "sigma^2(", effect_labels[[source]], ")"))
      },
      FUN.VALUE = character(1)
    )
    cat("  ", left_sides[[effect]], " = ", paste(terms, collapse = " + "), "\n", sep = "")
  }
  return(invisible(NULL))
}

# D studies from posterior draws against the plug-in D study of each draw, the
# checks and printed notes of method = "bayes", and Bayesian estimates against
# the posterior on a grid. Each sampled quantile is within 8% of the grid: from
# 20000 draws, its largest error was below 4% with each of 12 seeds, and a
# prior with twice the default scale moves a quantile by up to 11%.

# A G study with posterior draws of its variance components, as a Bayesian G
# study has them.
with_draws <- function(g, draws, level = 0.95) {
  g$variance_component_draws <- draws
  g$bayes_settings <- list(level = level)
  return(g)
}

# A G study whose posterior draws all equal its estimates. The posterior is
# then a point mass, so every D-study estimate from the draws must equal the
# plug-in estimate, with a standard error of zero.
with_point_mass_draws <- function(g, n_draws = 50) {
  estimates <- g$variance_component
  draws <- matrix(
    estimates, nrow = n_draws, ncol = length(estimates), byrow = TRUE,
    dimnames = list(NULL, names(estimates))
  )
  return(with_draws(g, draws))
}

# The D study from point-mass draws equals the D study without draws.
expect_plug_in_d_study <- function(g, universe_sizes, d_sizes = NULL, d_design = NULL) {
  plug_in <- d_study_general(
    g, object_of_measurement = "p", universe_sizes = universe_sizes,
    d_sizes = d_sizes, d_design = d_design
  )
  from_draws <- d_study_general(
    with_point_mass_draws(g), object_of_measurement = "p",
    universe_sizes = universe_sizes, d_sizes = d_sizes, d_design = d_design
  )
  expect_equal(
    from_draws$coefficient_results$estimate, plug_in$coefficient_results$estimate
  )
  expect_equal(from_draws$coefficient_results$standard_error, rep(0, 5))
  expect_equal(
    from_draws$d_mean_results$mean_score_estimate,
    plug_in$d_mean_results$mean_score_estimate
  )
  expect_equal(
    from_draws$d_mean_results$standard_error, rep(0, nrow(plug_in$d_mean_results))
  )
  expect_equal(from_draws$d_mean_results$divisor, plug_in$d_mean_results$divisor)
}

expect_options_error <- function(bayes_options, message) {
  expect_error(.resolve_bayes_options(bayes_options), message, fixed = TRUE)
}

expect_scale_error <- function(scores, scale) {
  dat <- brennan_p_x_i
  dat$score <- scores
  expect_error(
    g_study_general(dat, design = "p x i", method = "bayes"),
    paste0("must be finite and above zero. It is ", scale, "."),
    fixed = TRUE
  )
}

# The posterior quantiles of the variance components of a balanced p x i
# design on a grid of standard deviations, without sampling. The sums of
# squares are independent, and the likelihood of SS(alpha) is that of
# EMS(alpha) times a chi-square with df(alpha) degrees of freedom (Brennan
# 2001, Eq. 3.24) when the prior of the grand mean is flat; the Student-t
# prior that brms gives the grand mean is flat over the range of the data.
# Each standard deviation has the default half Student-t prior.
grid_posterior_p_x_i <- function(dat, probabilities) {
  n_p <- nlevels(factor(dat$p))
  n_i <- nlevels(factor(dat$i))
  grand_mean <- mean(dat$score)
  ss_p <- n_i * sum((tapply(dat$score, dat$p, mean) - grand_mean)^2)
  ss_i <- n_p * sum((tapply(dat$score, dat$i, mean) - grand_mean)^2)
  ss_pi <- sum((dat$score - grand_mean)^2) - ss_p - ss_i
  df_p <- n_p - 1
  df_i <- n_i - 1
  prior_scale <- sd(dat$score)
  log_prior <- function(sd_values) {
    return(dt(sd_values / prior_scale, df = 3, log = TRUE))
  }
  log_likelihood <- function(ss, df, ems) {
    return(-(df * log(ems) + ss / ems) / 2)
  }

  sd_p <- seq(0.0025, 2.5, length.out = 500)
  sd_i <- sd_p
  sd_pi <- seq(0.25, 0.55, length.out = 121)
  weight_p <- numeric(length(sd_p))
  weight_i <- numeric(length(sd_i))
  weight_pi <- numeric(length(sd_pi))
  # Adding 60 to each log density keeps exp() from underflowing; the constant
  # cancels when the weights are normalized.
  for (k in seq_along(sd_pi)) {
    residual <- sd_pi[[k]]^2
    density_p <- exp(
      log_likelihood(ss_p, df_p, residual + n_i * sd_p^2) + log_prior(sd_p) + 60
    )
    density_i <- exp(
      log_likelihood(ss_i, df_i, residual + n_p * sd_i^2) + log_prior(sd_i) + 60
    )
    density_pi <- exp(
      log_likelihood(ss_pi, df_p * df_i, residual) + log_prior(sd_pi[[k]]) + 60
    )
    weight_p <- weight_p + density_p * sum(density_i) * density_pi
    weight_i <- weight_i + density_i * sum(density_p) * density_pi
    weight_pi[[k]] <- sum(density_p) * sum(density_i) * density_pi
  }

  grid_quantiles <- function(sd_values, weights) {
    cumulative <- cumsum(weights) / sum(weights)
    return(vapply(
      probabilities,
      function(probability) sd_values[[which(cumulative >= probability)[[1]]]]^2,
      numeric(1)
    ))
  }
  return(rbind(
    p = grid_quantiles(sd_p, weight_p),
    i = grid_quantiles(sd_i, weight_i),
    pi = grid_quantiles(sd_pi, weight_pi)
  ))
}

test_that("D studies from point-mass draws equal the plug-in D studies", {
  g_prt <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  expect_plug_in_d_study(g_prt, c(r = Inf, t = Inf), c(r = 4, t = 3))
  expect_plug_in_d_study(g_prt, c(r = 8, t = 3), c(r = 4, t = 3))

  g_fixed <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t")
  expect_plug_in_d_study(g_fixed, c(r = Inf), c(r = 4))

  g_pio <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
  expect_plug_in_d_study(g_pio, c(i = Inf, o = Inf), c(i = 4, o = 2))
  expect_plug_in_d_study(g_pio, c(i = Inf, o = Inf), c(i = 4, o = 2), "p x (i:o)")
})

test_that("The G-study pattern from point-mass draws equals the plug-in D study", {
  g <- g_study_general(brennan_i_in_p_unbalanced, design = "i:p")
  expect_plug_in_d_study(g, c(i = Inf))

  g <- g_study_general(brennan_p_x_i_in_h_unbalanced, design = "p x (i:h)")
  expect_plug_in_d_study(g, c(i = Inf, h = Inf))

  g <- g_study_general(brennan_p_x_i_missing, design = "p x i")
  expect_plug_in_d_study(g, c(i = Inf))
})

test_that("Each draw of a D study equals the plug-in D study of that draw", {
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  draws <- rbind(
    c(p = 0.2, t = 0.1, "r:t" = 0.9, pt = 0.3, "pr:t" = 2.0),
    c(p = 0.9, t = 0.6, "r:t" = 0.2, pt = 0.8, "pr:t" = 3.1),
    c(p = 0.5, t = 0.0, "r:t" = 0.5, pt = 0.1, "pr:t" = 2.4),
    c(p = 0.1, t = 0.4, "r:t" = 0.7, pt = 0.6, "pr:t" = 1.6)
  )
  from_draws <- d_study_general(
    with_draws(g, draws, level = 0.9), object_of_measurement = "p",
    universe_sizes = c(r = Inf, t = 8), d_sizes = c(r = 4, t = 3)
  )

  plug_in_coefficients <- matrix(NA_real_, nrow = nrow(draws), ncol = 5)
  plug_in_d_mean <- matrix(NA_real_, nrow = nrow(draws), ncol = ncol(draws))
  for (k in seq_len(nrow(draws))) {
    g$variance_component <- draws[k, ]
    plug_in <- d_study_general(
      g, object_of_measurement = "p",
      universe_sizes = c(r = Inf, t = 8), d_sizes = c(r = 4, t = 3)
    )
    plug_in_coefficients[k, ] <- plug_in$coefficient_results$estimate
    plug_in_d_mean[k, ] <- plug_in$d_mean_results$mean_score_estimate
  }
  expect_equal(unname(from_draws$coefficient_draws), plug_in_coefficients)

  # The tables summarize the draws: medians, standard deviations, and the 5%
  # and 95% quantiles for a level of 0.9.
  coefficients <- from_draws$coefficient_results
  expect_identical(
    names(coefficients),
    c("quantity", "estimate", "standard_error", "lower_90", "upper_90")
  )
  expect_equal(coefficients$estimate, apply(plug_in_coefficients, 2, median))
  expect_equal(coefficients$standard_error, apply(plug_in_coefficients, 2, sd))
  expect_equal(coefficients$lower_90, apply(plug_in_coefficients, 2, quantile, probs = 0.05))
  expect_equal(coefficients$upper_90, apply(plug_in_coefficients, 2, quantile, probs = 0.95))
  d_mean <- from_draws$d_mean_results
  expect_equal(d_mean$mean_score_estimate, apply(plug_in_d_mean, 2, median))
  expect_equal(d_mean$lower_90, apply(plug_in_d_mean, 2, quantile, probs = 0.05))
  expect_identical(
    names(from_draws$single_obs_results)[6:7], c("lower_90", "upper_90")
  )
  expect_identical(from_draws$level, 0.9)
})

test_that("A coefficient that is undefined in a draw is NA", {
  g <- g_study_general(brennan_p_x_i, design = "p x i")
  draws <- rbind(c(p = 0.1, i = 0.2, pi = 0.3), c(p = 0, i = 0, pi = 0))
  d <- d_study_general(
    with_draws(g, draws), object_of_measurement = "p",
    universe_sizes = c(i = Inf), d_sizes = c(i = 12)
  )
  is_ratio <- d$coefficient_results$quantity %in% c("G", "Phi")
  expect_true(all(is.na(d$coefficient_results$estimate[is_ratio])))
  expect_true(all(is.na(d$coefficient_results$lower_95[is_ratio])))
  expect_false(anyNA(d$coefficient_results$estimate[!is_ratio]))
})

test_that("d_study() and compare_d_studies() carry the draws of g_study()", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
  )
  g$general <- with_point_mass_draws(g$general)
  g$general$variance_component_draws[, "p"] <- seq(0.1, 1, length.out = 50)
  general <- d_study_general(
    g$general, object_of_measurement = "p",
    universe_sizes = c(a = Inf, b = Inf), d_sizes = c(a = 3, b = 4)
  )

  d <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
  expect_equal(d$coefficient_draws, general$coefficient_draws)
  expect_equal(d$coefficient_results, general$coefficient_results)
  # d_study() reports PxF1 before F2:F1; d_study_general() has them reversed.
  expect_identical(d$d_mean_results$effect, c("P", "F1", "PxF1", "F2:F1", "F2:F1xP"))
  expect_equal(d$d_mean_results$lower_95, general$d_mean_results$lower_95[c(1, 2, 4, 3, 5)])
  expect_identical(d$level, 0.95)

  comparison <- compare_d_studies(
    g, universe_sizes = c(F1 = Inf, F2 = Inf),
    d_sizes = data.frame(F1 = c(1, 3), F2 = c(12, 4))
  )
  expect_equal(comparison$coefficient_draws[[2]], general$coefficient_draws)
  expect_identical(comparison$level, 0.95)
  expect_equal(
    comparison$scenarios$lower_G[[2]],
    quantile(general$coefficient_draws[, "G"], probs = 0.025, names = FALSE)
  )
  expect_equal(
    comparison$scenarios$upper_Phi[[2]],
    quantile(general$coefficient_draws[, "Phi"], probs = 0.975, names = FALSE)
  )
})

test_that("find_d_designs() selects by the posterior probability of the target", {
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  d_sizes <- data.frame(r = c(12, 4, 2), t = c(1, 3, 6))
  # Only the person variance differs across the draws, from 0.1 to 1, so
  # G = p / (p + delta) is at least 0.55 where p is at least
  # delta * 0.55 / 0.45, with delta the relative error variance of the
  # plug-in D study.
  person_variance <- seq(0.1, 1, length.out = 50)
  delta <- vapply(
    seq_len(nrow(d_sizes)),
    function(k) {
      d <- d_study_general(
        g, object_of_measurement = "p",
        universe_sizes = c(r = Inf, t = Inf), d_sizes = unlist(d_sizes[k, ])
      )
      return(d$coefficient_results$estimate[[2]])
    },
    numeric(1)
  )
  probability <- vapply(
    delta, function(value) mean(person_variance >= value * 0.55 / 0.45), numeric(1)
  )

  g_draws <- with_point_mass_draws(g)
  g_draws$variance_component_draws[, "p"] <- person_variance
  comparison <- compare_d_studies_general(
    g_draws, object_of_measurement = "p",
    universe_sizes = c(r = Inf, t = Inf), d_sizes = d_sizes
  )
  selection <- find_d_designs(comparison, coefficient = "G", target = 0.55, probability = 0.5)
  expect_equal(selection$posterior_probability, probability)
  expect_identical(selection$meets_target, probability >= 0.5)
  expect_identical(selection$designs$scenario, which(probability >= 0.5))
  expect_equal(selection$designs$posterior_probability, probability[probability >= 0.5])

  # Without probability, the posterior medians are compared with the target.
  selection <- find_d_designs(comparison, coefficient = "G", target = 0.55)
  expect_identical(selection$meets_target, comparison$scenarios$G >= 0.55)
  expect_equal(selection$posterior_probability, probability)
  expect_output(print(selection), "Target: posterior median of G >= 0.55", fixed = TRUE)
  expect_output(print(selection), "lower_G", fixed = TRUE)

  # Phi is below G in every draw, so its probability is no larger.
  selection <- find_d_designs(comparison, coefficient = "Phi", target = 0.55)
  expect_true(all(selection$posterior_probability <= probability))
})

test_that("A failed D study has no draws and no posterior probability", {
  g <- with_point_mass_draws(g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)"))
  expect_warning(
    comparison <- compare_d_studies_general(
      g, object_of_measurement = "p",
      universe_sizes = c(r = Inf, t = 5), d_sizes = data.frame(r = c(4, 4), t = c(3, 6))
    ),
    "1 of 2 D-study scenarios failed"
  )
  expect_identical(comparison$scenarios$status, c("ok", "error"))
  expect_length(comparison$coefficient_draws, 2)
  expect_null(comparison$coefficient_draws[[2]])
  expect_true(is.na(comparison$scenarios$lower_G[[2]]))

  selection <- find_d_designs(comparison, coefficient = "G", target = 0.5, probability = 0.5)
  expect_identical(selection$posterior_probability, c(1, NA))
  expect_identical(selection$meets_target, c(TRUE, NA))
})

test_that("The notes and diagnostics of a Bayesian G study are printed", {
  # A G study with the elements that method = "bayes" adds, without a fit.
  g <- g_study_general(brennan_p_x_r_in_t, design = "p x (r:t)")
  g$method <- "bayes"
  g$bayes_settings <- list(
    prior_lines = "residual standard deviation: half-student_t(3, 0, 2.043)",
    prior_scale = 2.043,
    level = 0.9,
    n_levels = c(p = 10, t = 3, "r:t" = 12, pt = 30)
  )
  g$fit_diagnostics <- list(
    max_rhat = 1.0234, min_ess = 231.4, n_divergent = 2,
    chains = 4, iter = 4000, warmup = 2000, n_draws = 8000,
    issues = c(
      "The largest R-hat of the standard deviations, 1.023, is above 1.01.",
      "Divergent transitions: 2."
    ),
    brms_messages = character(0),
    brms_warnings = "There were 2 divergent transitions after warmup. See\nhttps://mc-stan.org",
    has_issue = TRUE
  )
  expected_lines <- c(
    "\n\nPriors:\n  residual standard deviation: half-student_t(3, 0, 2.043)",
    "Sampling: 4 chains, 4000 iterations (2000 warmup), 8000 draws",
    "Standard deviations: largest R-hat 1.023, smallest effective sample size 231",
    "deviations, with 90% equal-tailed credible intervals.",
    "Effects with fewer than 5 levels, whose estimates depend on the priors: t (3)",
    "G-study fit diagnostics:\n  The largest R-hat of the standard deviations, 1.023,",
    "    - There were 2 divergent transitions after warmup. See\n      https://mc-stan.org"
  )
  for (expected_line in expected_lines) {
    expect_output(print(g), expected_line, fixed = TRUE)
  }

  g$fit_diagnostics$chains <- 1
  g$fit_diagnostics$issues <- character(0)
  g$fit_diagnostics$brms_warnings <- character(0)
  g$fit_diagnostics$has_issue <- FALSE
  g$bayes_settings$n_levels <- c(p = 10, t = 6)
  printed <- paste(capture.output(print(g)), collapse = "\n")
  expect_match(printed, "Sampling: 1 chain, 4000 iterations", fixed = TRUE)
  expect_no_match(printed, "G-study fit diagnostics", fixed = TRUE)
  expect_no_match(printed, "fewer than 5 levels", fixed = TRUE)

  # g_study() names the effects with few levels by their columns.
  g_keyword <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
  )
  g_keyword$method <- "bayes"
  g_keyword$bayes_settings <- g$bayes_settings
  g_keyword$bayes_settings$n_levels <- c(P = 10, F1 = 3)
  g_keyword$fit_diagnostics <- g$fit_diagnostics
  expect_output(print(g_keyword), "depend on the priors: t (3)", fixed = TRUE)
})

test_that("The elements that method = \"bayes\" adds come last in every result", {
  bayes_elements <- c("variance_component_draws", "bayes_settings")
  g <- g_study_general(brennan_p_x_i, design = "p x i")
  expect_identical(tail(names(g), 2), bayes_elements)
  expect_null(g$variance_component_draws)
  g <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed", method = "anova")
  expect_identical(tail(names(g), 2), bayes_elements)

  comparison <- compare_d_studies(
    g, universe_sizes = c(F1 = Inf), d_sizes = data.frame(F1 = c(6, 12))
  )
  expect_identical(tail(names(comparison), 2), c("coefficient_draws", "level"))
  selection <- find_d_designs(comparison, coefficient = "G", target = 0.8)
  expect_identical(tail(names(selection), 2), c("probability", "posterior_probability"))
})

test_that("D studies keep the level of the G study", {
  g <- g_study(
    brennan_p_x_r_in_t, p = "p", f1 = "t", f2 = "r", design = "partially_nested"
  )
  g$general <- with_draws(
    g$general, with_point_mass_draws(g$general, 20)$variance_component_draws, level = 0.8
  )
  g$general$variance_component_draws[, "p"] <- seq(0.1, 1, length.out = 20)
  d <- d_study(g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = c(F1 = 3, F2 = 4))
  expect_identical(d$level, 0.8)
  expect_identical(names(d$coefficient_results)[4:5], c("lower_80", "upper_80"))
  expect_equal(
    d$coefficient_results$lower_80[[4]],
    quantile(d$coefficient_draws[, "G"], probs = 0.1, names = FALSE)
  )
  expect_output(print(d), "80% equal-tailed credible intervals", fixed = TRUE)

  comparison <- compare_d_studies(
    g, universe_sizes = c(F1 = Inf, F2 = Inf), d_sizes = data.frame(F1 = 3, F2 = 4)
  )
  expect_identical(comparison$level, 0.8)
  expect_equal(
    comparison$scenarios$upper_G[[1]],
    quantile(d$coefficient_draws[, "G"], probs = 0.9, names = FALSE)
  )
})

test_that("The settings in bayes_options are checked before the model is fitted", {
  expect_error(
    g_study_general(
      brennan_p_x_i, design = "p x i", method = "bayes", bayes_options = list(Iter = 2000)
    ),
    "Names are case-sensitive: Iter is not iter.",
    fixed = TRUE
  )
  expect_options_error(
    list(data = brennan_p_x_i), "bayes_options must be a named list that sets any of chains"
  )
  expect_options_error(list(2000), "bayes_options must be a named list")
  expect_options_error(data.frame(chains = 2), "bayes_options must be a named list")
  expect_options_error(
    list(Iter = 2000), "such as list(iter = 2000). Names are case-sensitive: Iter is not iter."
  )
  expect_options_error(
    list(iter = 3000, iter = 5000), "These settings are repeated in bayes_options: iter."
  )
  expect_options_error(list(iter = 1), "bayes_options$iter must be at least 2")
  expect_options_error(
    list(iter = 2000, warmup = 2000),
    "bayes_options$iter must be larger than bayes_options$warmup."
  )
  expect_options_error(
    list(chains = 0), "bayes_options$chains must be a single whole number from 1 to 2147483647."
  )
  expect_options_error(
    list(chains = 3e9), "bayes_options$chains must be a single whole number from 1 to"
  )
  expect_options_error(
    list(cores = "2"), "bayes_options$cores must be a single whole number from 1 to"
  )
  expect_options_error(
    list(max_treedepth = -1), "bayes_options$max_treedepth must be a single whole number"
  )
  expect_options_error(
    list(seed = "abc"), "bayes_options$seed must be a single whole number or NA."
  )
  expect_options_error(
    list(seed = 1e10), "bayes_options$seed must be a single whole number or NA."
  )
  expect_options_error(
    list(seed = NaN), "bayes_options$seed must be a single whole number or NA."
  )
  expect_options_error(
    list(adapt_delta = 1),
    "bayes_options$adapt_delta must be a single number between 0 and 1, exclusive."
  )
  expect_options_error(
    list(level = 0), "bayes_options$level must be a single number between 0 and 1, exclusive."
  )
  expect_options_error(
    list(prior_scale = 0),
    "bayes_options$prior_scale must be a single finite number above zero."
  )
  expect_options_error(
    list(prior = "normal(0, 1)"), "bayes_options$prior must be made with brms::set_prior()"
  )
  expect_error(
    g_study_general(brennan_p_x_i, design = "p x i", bayes_options = list(seed = 1)),
    "bayes_options is only used with method = \"bayes\".",
    fixed = TRUE
  )
  expect_error(
    g_study(
      brennan_p_x_i, p = "p", f1 = "i", design = "crossed", method = "anova",
      bayes_options = list(seed = 1)
    ),
    "bayes_options is only used with method = \"bayes\".",
    fixed = TRUE
  )
})

test_that("The scores are checked before the model is fitted", {
  scores <- brennan_p_x_i$score
  # Equal scores, a score that is not finite, and scores whose variance
  # overflows.
  expect_scale_error(rep(3, length(scores)), "0")
  expect_scale_error(replace(scores, 1, Inf), "NaN")
  expect_scale_error(scores * 1e307, "Inf")

  dat <- brennan_p_x_i
  dat$score <- as.character(dat$score)
  expect_error(
    g_study_general(dat, design = "p x i", method = "bayes"),
    "but the score column, score, is of class character.",
    fixed = TRUE
  )
})

test_that("Priors of standard deviations centered at zero are printed as half distributions", {
  expect_identical(
    .is_centered_at_zero(c(
      "normal(0, 1)", "normal(0.0, 1)", "student_t(3, 0, 2.043)", "cauchy(0,2)",
      "normal(1, 2)", "normal(0.5, 1)", "student_t(3, 10, 2)", "exponential(1)"
    )),
    c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )
})

test_that("A prior for one effect is added to the defaults and names a model term", {
  skip_if_not_installed("brms", minimum_version = "2.19.0")
  default_prior <- "student_t(3, 0, gt_prior_scale)"
  prior <- .brms_prior(brms::set_prior("normal(0, 1)", class = "sd", group = "p"))
  expect_identical(prior$prior, c(default_prior, default_prior, "normal(0, 1)"))
  expect_identical(prior$class, c("sd", "sigma", "sd"))
  expect_identical(prior$group, c("", "", "p"))
  # A prior for a whole class replaces the default of that class.
  prior <- .brms_prior(brms::set_prior("exponential(1)", class = "sd"))
  expect_identical(prior$prior, c(default_prior, "exponential(1)"))
  expect_identical(prior$class, c("sigma", "sd"))

  # g_study() writes the model formula with P and F1.
  expect_error(
    g_study(
      brennan_p_x_i, p = "p", f1 = "i", design = "crossed", method = "bayes",
      bayes_options = list(prior = brms::set_prior("normal(0, 1)", class = "sd", group = "p"))
    ),
    "the model formula: P, F1. Not a random effect: p. Names are case-sensitive: p is not P.",
    fixed = TRUE
  )
})

test_that("Selecting by probability requires posterior draws and a probability", {
  g <- g_study_general(brennan_p_x_i, design = "p x i")
  comparison <- compare_d_studies_general(
    g, object_of_measurement = "p", universe_sizes = c(i = Inf), d_sizes = list(i = 12)
  )
  expect_error(
    find_d_designs(comparison, coefficient = "G", target = 0.8, probability = 0.9),
    "probability requires a comparison from a G study with method = \"bayes\".",
    fixed = TRUE
  )
  comparison <- compare_d_studies_general(
    with_point_mass_draws(g), object_of_measurement = "p",
    universe_sizes = c(i = Inf), d_sizes = list(i = 12)
  )
  expect_error(
    find_d_designs(comparison, coefficient = "G", target = 0.8, probability = 1.5),
    "probability must be a single finite number between 0 and 1.",
    fixed = TRUE
  )
})

test_that("Bayesian estimates agree with the grid posterior for Brennan Table 2.2 (p x i)", {
  skip_on_cran()
  skip_if_not_installed("brms", minimum_version = "2.19.0")
  dat <- brennan_p_x_i
  g <- g_study(
    dat, p = "p", f1 = "i", design = "crossed", method = "bayes",
    bayes_options = list(iter = 10000, seed = 1)
  )
  # g_study() renames the effects of the draws, which g$general has as p, a,
  # and pa.
  draws <- g$variance_component_draws
  expect_identical(colnames(draws), c("P", "F1", "PxF1"))
  expect_identical(nrow(draws), 20000L)
  expect_identical(colnames(g$general$variance_component_draws), c("p", "a", "pa"))
  expect_equal(unname(g$general$variance_component_draws), unname(draws))
  expect_identical(g$bayes_settings$n_levels, c(P = 10, F1 = 12))

  probabilities <- c(0.025, 0.5, 0.975)
  grid <- grid_posterior_p_x_i(dat, probabilities)
  sampled <- t(apply(draws, 2, quantile, probs = probabilities))
  expect_lt(max(abs(unname(sampled) / unname(grid) - 1)), 0.08)
  expect_equal(g$variance_components$variance_estimate, unname(sampled[, 2]))
  expect_equal(g$variance_components$lower_95, unname(sampled[, 1]))
  expect_equal(g$variance_components$upper_95, unname(sampled[, 3]))

  expect_output(print(g), "half-student_t(3, 0, 0.4987)", fixed = TRUE)
  d <- d_study(g, universe_sizes = c(F1 = Inf), d_sizes = c(F1 = 12))
  expect_identical(dim(d$coefficient_draws), c(20000L, 5L))
  expect_equal(
    d$coefficient_draws[, "G"], draws[, "P"] / (draws[, "P"] + draws[, "PxF1"] / 12)
  )
})

test_that("bayes_options reach brms, and the reported table follows the draws", {
  skip_on_cran()
  skip_if_not_installed("brms", minimum_version = "2.19.0")
  # A short run with the sampler settings, the seed, the level, and the prior
  # scale changed; a prior for one effect and one for the residual; and tasks
  # fixed.
  g <- g_study_general(
    brennan_p_x_r_in_t, design = "p x (r:t)", fixed = "t", method = "bayes",
    bayes_options = list(
      chains = 2, iter = 200, warmup = 100, seed = 7, adapt_delta = 0.95,
      max_treedepth = 10, level = 0.8, prior_scale = 2,
      prior = c(
        brms::set_prior("normal(0, 0.5)", class = "sd", group = "r:t"),
        brms::set_prior("normal(1, 2)", class = "sigma")
      )
    )
  )
  draws <- g$variance_component_draws
  expect_identical(nrow(draws), 200L)
  expect_identical(g$bayes_settings$prior_scale, 2)
  expect_identical(g$fit$stanvars$gt_prior_scale$sdata, 2)
  if (identical(g$fit$backend, "rstan")) {
    stan_args <- g$fit$fit@stan_args[[1]]
    expect_identical(as.numeric(stan_args$seed), 7)
    expect_identical(stan_args$control$adapt_delta, 0.95)
    expect_identical(stan_args$control$max_treedepth, 10)
  }
  # A prior centered at zero is printed as a half distribution.
  expected_lines <- c(
    "standard deviations of the random effects: half-student_t(3, 0, 2)",
    "standard deviation of r:t: half-normal(0, 0.5)",
    "residual standard deviation: normal(1, 2)",
    "Sampling: 2 chains, 200 iterations (100 warmup), 200 draws"
  )
  for (expected_line in expected_lines) {
    expect_output(print(g), expected_line, fixed = TRUE)
  }

  # The draws are the squares of the sampled standard deviations.
  expect_equal(
    unname(draws[, "r:t"]), unname(as.matrix(g$fit, variable = "sd_r:t__Intercept")[, 1]^2)
  )
  expect_equal(unname(draws[, "pr:t"]), unname(as.matrix(g$fit, variable = "sigma")[, 1]^2))

  # The reported p row applies Eq. 3.34 to each draw, sigma^2(p) +
  # sigma^2(pt) / 3, and its interval has the level of bayes_options.
  reported_p <- draws[, "p"] + draws[, "pt"] / 3
  vc_table <- g$variance_components
  expect_identical(names(vc_table)[5:6], c("lower_80", "upper_80"))
  expect_equal(vc_table$variance_estimate[[1]], median(reported_p))
  expect_equal(vc_table$lower_80[[1]], quantile(reported_p, probs = 0.1, names = FALSE))
  expect_equal(vc_table$upper_80[[1]], quantile(reported_p, probs = 0.9, names = FALSE))

  # 200 draws cannot reach an effective sample size of 400.
  expect_lt(g$fit_diagnostics$min_ess, 400)
  expect_true(any(grepl("is below 400", g$fit_diagnostics$issues, fixed = TRUE)))
  expect_output(print(g), "G-study fit diagnostics", fixed = TRUE)
})

# The settings of method = "bayes" that `bayes_options` can change, with
# their defaults. A `warmup` of NULL is half of `iter`, and a `prior_scale`
# of NULL is the standard deviation of the scores.
.bayes_option_defaults <- list(
  chains = 4,
  iter = 4000,
  warmup = NULL,
  seed = NA,
  cores = 1,
  adapt_delta = 0.99,
  max_treedepth = 12,
  prior_scale = NULL,
  prior = NULL,
  level = 0.95
)

# The oldest brms whose R-hat and effective sample size are those of the
# posterior package.
.brms_minimum_version <- "2.19.0"

# Check `bayes_options` and return the settings: the defaults, with those the
# user gives. An empty list gives the defaults.
.resolve_bayes_options <- function(bayes_options) {
  if (is.null(bayes_options)) {
    bayes_options <- list()
  }
  option_names <- names(bayes_options)
  valid_names <- names(.bayes_option_defaults)
  if (!is.list(bayes_options) || is.data.frame(bayes_options) ||
      (length(bayes_options) > 0 && is.null(option_names)) ||
      !all(option_names %in% valid_names)) {
    stop(
      "bayes_options must be a named list that sets any of ",
      paste(valid_names, collapse = ", "), ", such as list(iter = 2000).",
      .case_note(option_names, valid_names),
      call. = FALSE
    )
  }
  if (anyDuplicated(option_names) > 0) {
    stop(
      "These settings are repeated in bayes_options: ",
      paste(unique(option_names[duplicated(option_names)]), collapse = ", "), ".",
      call. = FALSE
    )
  }
  result <- .bayes_option_defaults
  result[option_names] <- bayes_options

  for (name in c("chains", "iter", "cores", "max_treedepth")) {
    .check_whole_number(result[[name]], paste0("bayes_options$", name))
  }
  if (result$iter < 2) {
    stop(
      "bayes_options$iter must be at least 2, for at least one warmup and one ",
      "sampling iteration.",
      call. = FALSE
    )
  }
  if (is.null(result$warmup)) {
    result$warmup <- result$iter %/% 2
  }
  .check_whole_number(result$warmup, "bayes_options$warmup")
  if (result$iter <= result$warmup) {
    stop("bayes_options$iter must be larger than bayes_options$warmup.", call. = FALSE)
  }

  # A seed of NULL or NA leaves the seed to brms; NaN is an error.
  seed <- result$seed
  if (is.null(seed) ||
      (is.atomic(seed) && length(seed) == 1 && is.na(seed) && !is.nan(seed))) {
    result$seed <- NA
  } else if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) ||
             seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("bayes_options$seed must be a single whole number or NA.", call. = FALSE)
  }

  for (name in c("adapt_delta", "level")) {
    value <- result[[name]]
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value) ||
        value <= 0 || value >= 1) {
      stop(
        "bayes_options$", name, " must be a single number between 0 and 1, exclusive.",
        call. = FALSE
      )
    }
  }
  if (!is.null(result$prior_scale)) {
    .check_positive_number(result$prior_scale, "bayes_options$prior_scale")
  }
  if (!is.null(result$prior) && !inherits(result$prior, "brmsprior")) {
    stop(
      "bayes_options$prior must be made with brms::set_prior() or brms::prior().",
      call. = FALSE
    )
  }
  return(result)
}

.check_brms_installed <- function() {
  version_check <- list(op = ">=", version = .brms_minimum_version)
  if (!requireNamespace("brms", versionCheck = version_check, quietly = TRUE)) {
    stop(
      "method = \"bayes\" requires the brms package, version ", .brms_minimum_version,
      " or later. Install it with install.packages(\"brms\").",
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# The scale of the default priors: `prior_scale`, or the standard deviation
# of the observed scores when it is NULL.
.resolve_prior_scale <- function(scores, prior_scale) {
  if (!is.null(prior_scale)) {
    return(prior_scale)
  }
  result <- sd(scores)
  if (!is.finite(result) || result <= 0) {
    stop(
      "The default prior scale is the standard deviation of the scores, which ",
      "must be finite and above zero. It is ", format(result), ". Check that ",
      "the scores are finite and not all equal, or set bayes_options$prior_scale.",
      call. = FALSE
    )
  }
  return(result)
}

# Check that the prior for each random effect names a term of the model
# formula, `group_names`, so that brms does not stop with a message about
# Stan parameters.
.check_prior_groups <- function(user_prior, group_names) {
  if (is.null(user_prior)) {
    return(invisible(NULL))
  }
  groups <- unique(user_prior$group[nzchar(user_prior$group)])
  unknown_groups <- setdiff(groups, group_names)
  if (length(unknown_groups) > 0) {
    stop(
      "The group of a prior in bayes_options$prior must be a random effect of the ",
      "model formula: ", paste(group_names, collapse = ", "), ". Not a random ",
      "effect: ", paste(unknown_groups, collapse = ", "), ".",
      .case_note(unknown_groups, group_names),
      call. = FALSE
    )
  }
  return(invisible(NULL))
}

# The priors of the standard deviations: half Student-t with 3 degrees of
# freedom and scale `gt_prior_scale`, which is passed to Stan as data. A user
# prior for a whole class ("sd" or "sigma") replaces the default of that
# class; a user prior for one random effect is added, and the other effects
# keep the default.
.brms_prior <- function(user_prior) {
  result <- c(
    brms::set_prior("student_t(3, 0, gt_prior_scale)", class = "sd"),
    brms::set_prior("student_t(3, 0, gt_prior_scale)", class = "sigma")
  )
  if (is.null(user_prior)) {
    return(result)
  }

  is_whole_class <- !nzchar(user_prior$group) & !nzchar(user_prior$coef)
  is_replaced <- result$class %in% user_prior$class[is_whole_class]
  return(c(result[!is_replaced, ], user_prior))
}

# Whether a normal, Cauchy, or Student-t prior is centered at zero. brms
# truncates the priors of standard deviations at zero, so these are half
# distributions.
.is_centered_at_zero <- function(distribution) {
  compact <- gsub(" ", "", distribution, fixed = TRUE)
  result <- grepl("^((normal|cauchy)\\(|student_t\\([^,]+,)-?0(\\.0*)?,", compact)
  return(result)
}

# The printed name of a row of brms::prior_summary().
.brms_prior_label <- function(class, coef, group) {
  if (class == "Intercept") {
    return("intercept")
  }
  if (class == "b") {
    if (nzchar(coef)) {
      return(paste0("fixed effect ", coef))
    }
    return("fixed effects")
  }
  if (class == "sigma") {
    return("residual standard deviation")
  }
  if (nzchar(group)) {
    return(paste0("standard deviation of ", group))
  }
  return("standard deviations of the random effects")
}

# The priors of a fit as printed lines, read from the fit so that they are
# the priors brms used: the intercept, the fixed effects of facets fixed in
# the G study (flat unless the user gives a prior), the standard deviations,
# and the residual.
.brms_prior_lines <- function(fit, prior_scale) {
  used <- as.data.frame(brms::prior_summary(fit))
  is_all_fixed_effects <- used$class == "b" & !nzchar(used$coef)
  used <- used[nzchar(used$prior) | is_all_fixed_effects, ]
  # brms versions differ in the order of the rows.
  class_position <- match(used$class, c("Intercept", "b", "sd", "sigma"))
  is_one_effect <- nzchar(used$coef) | nzchar(used$group)
  used <- used[order(class_position, is_one_effect), ]

  result <- character(nrow(used))
  for (k in seq_len(nrow(used))) {
    class <- used$class[[k]]
    distribution <- gsub(
      "gt_prior_scale", format(signif(prior_scale, 4)), used$prior[[k]], fixed = TRUE
    )
    label <- .brms_prior_label(class, used$coef[[k]], used$group[[k]])
    if (class == "b" && !nzchar(distribution)) {
      distribution <- "flat"
    }
    if (class %in% c("sd", "sigma") && .is_centered_at_zero(distribution)) {
      distribution <- paste0("half-", distribution)
    }
    result[[k]] <- paste0(label, ": ", distribution)
  }
  return(result)
}

# R-hat and effective sample size of the standard deviations, the number of
# divergent transitions, and the sampler settings. `issues` has a sentence
# for each of them that calls for a longer or more careful run: an R-hat
# above 1.01, an effective sample size below 400, a value that brms cannot
# compute, or a divergent transition. brms messages are kept but are not an
# issue, except those that begin with "Warning", which is how the cmdstanr
# backend reports sampler problems.
.build_brms_fit_diagnostics <- function(fit, bayes_options, captured_messages,
                                        captured_warnings) {
  parameters <- grep("^sd_|^sigma$", brms::variables(fit), value = TRUE)
  n_draws <- brms::ndraws(fit)
  max_rhat <- max(brms::rhat(fit, pars = parameters))
  min_ess <- min(brms::neff_ratio(fit, pars = parameters)) * n_draws
  sampler <- brms::nuts_params(fit)
  n_divergent <- sum(sampler$Value[sampler$Parameter == "divergent__"])
  captured_messages <- unique(trimws(captured_messages))
  captured_warnings <- unique(trimws(captured_warnings))
  is_sampler_warning <- grepl("^Warning", captured_messages)
  captured_warnings <- c(captured_warnings, captured_messages[is_sampler_warning])
  captured_messages <- captured_messages[!is_sampler_warning]

  issues <- character()
  if (is.na(max_rhat) || is.na(min_ess)) {
    issues <- c(
      issues,
      "R-hat or the effective sample size of a standard deviation is missing."
    )
  }
  if (!is.na(max_rhat) && max_rhat > 1.01) {
    issues <- c(issues, paste0(
      "The largest R-hat of the standard deviations, ",
      format(round(max_rhat, 3), nsmall = 3), ", is above 1.01."
    ))
  }
  if (!is.na(min_ess) && min_ess < 400) {
    issues <- c(issues, paste0(
      "The smallest effective sample size of the standard deviations, ",
      round(min_ess), ", is below 400."
    ))
  }
  if (n_divergent > 0) {
    issues <- c(issues, paste0("Divergent transitions: ", n_divergent, "."))
  }
  has_issue <- length(issues) > 0 || length(captured_warnings) > 0

  return(list(
    max_rhat = max_rhat,
    min_ess = min_ess,
    n_divergent = n_divergent,
    chains = bayes_options$chains,
    iter = bayes_options$iter,
    warmup = bayes_options$warmup,
    n_draws = n_draws,
    issues = issues,
    brms_messages = captured_messages,
    brms_warnings = captured_warnings,
    has_issue = has_issue
  ))
}

# The model of method = "reml", fitted by Markov chain Monte Carlo with the
# priors of .brms_prior().
.fit_brms_with_diagnostics <- function(formula, dat_for_brms, bayes_options, prior_scale) {
  dat_for_brms <- dat_for_brms[!is.na(dat_for_brms$score), ]
  captured_messages <- character()
  captured_warnings <- character()

  fit <- withCallingHandlers(
    brms::brm(
      formula,
      data = dat_for_brms,
      family = gaussian(),
      prior = .brms_prior(bayes_options$prior),
      stanvars = brms::stanvar(prior_scale, name = "gt_prior_scale"),
      chains = bayes_options$chains,
      iter = bayes_options$iter,
      warmup = bayes_options$warmup,
      seed = bayes_options$seed,
      cores = bayes_options$cores,
      control = list(
        adapt_delta = bayes_options$adapt_delta,
        max_treedepth = bayes_options$max_treedepth
      ),
      refresh = 0,
      silent = 2
    ),
    message = function(message_condition) {
      captured_messages <<- c(captured_messages, conditionMessage(message_condition))
      invokeRestart("muffleMessage")
    },
    warning = function(warning_condition) {
      captured_warnings <<- c(captured_warnings, conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }
  )

  return(list(
    fit = fit,
    diagnostics = .build_brms_fit_diagnostics(
      fit,
      bayes_options = bayes_options,
      captured_messages = captured_messages,
      captured_warnings = captured_warnings
    )
  ))
}

# The printed priors and sampler settings of a Bayesian fit.
.print_brms_settings <- function(bayes_settings, diagnostics) {
  cat("Priors:\n")
  cat(paste0("  ", bayes_settings$prior_lines), sep = "\n")
  chain_word <- if (diagnostics$chains == 1) " chain, " else " chains, "
  cat(
    "Sampling: ", diagnostics$chains, chain_word, diagnostics$iter, " iterations (",
    diagnostics$warmup, " warmup), ", diagnostics$n_draws, " draws\n",
    sep = ""
  )
  cat(
    "Standard deviations: largest R-hat ", format(round(diagnostics$max_rhat, 3), nsmall = 3),
    ", smallest effective sample size ", round(diagnostics$min_ess), "\n",
    sep = ""
  )
  return(invisible(NULL))
}

.print_brms_fit_diagnostics <- function(diagnostics) {
  if (!diagnostics$has_issue) {
    return(invisible(NULL))
  }

  cat("\n")
  cat("G-study fit diagnostics:\n")
  for (issue in diagnostics$issues) {
    cat("  ", issue, "\n", sep = "")
  }

  if (length(diagnostics$brms_messages) > 0) {
    cat("  brms messages:\n")
    for (message in diagnostics$brms_messages) {
      cat("    - ", gsub("\n", "\n      ", message, fixed = TRUE), "\n", sep = "")
    }
  }

  if (length(diagnostics$brms_warnings) > 0) {
    cat("  brms warnings:\n")
    for (warning in diagnostics$brms_warnings) {
      cat("    - ", gsub("\n", "\n      ", warning, fixed = TRUE), "\n", sep = "")
    }
  }

  return(invisible(NULL))
}

# Posterior draws of the variance components, one column per effect. brms
# samples standard deviations, and their squares are the variances. The
# brms names are renamed to G-theory effect names, for example
# sd_p:r:t__Intercept to pr:t.
.extract_variance_component_draws <- function(fit, design) {
  effect_names <- c(design$random_effects, design$residual_effect)
  lmer_names <- unname(design$lmer_terms[design$random_effects])
  brms_names <- c(paste0("sd_", lmer_names, "__Intercept"), "sigma")

  sd_draws <- as.matrix(fit, variable = brms_names)[, brms_names, drop = FALSE]
  result <- sd_draws^2
  dimnames(result) <- list(NULL, effect_names)
  return(result)
}

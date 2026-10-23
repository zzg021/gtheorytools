# The lmer formula: a random intercept for each random effect, and the main
# effect of a facet fixed in the G study as a fixed effect.
.lmer_formula <- function(score_column, random_effects, fixed_effects = NULL) {
  random_terms <- paste0("(1 | ", .lmer_effect_name(random_effects), ")")
  right_side <- paste(random_terms, collapse = " + ")
  if (!is.null(fixed_effects)) {
    right_side <- paste(paste(fixed_effects, collapse = " + "), "+", right_side)
  }
  return(as.formula(paste(score_column, "~", right_side)))
}

.build_lmer_fit_diagnostics <- function(fit, captured_messages, captured_warnings) {
  is_singular <- lme4::isSingular(fit)
  convergence_code <- fit@optinfo$conv$opt
  captured_messages <- unique(trimws(captured_messages))
  captured_warnings <- unique(trimws(captured_warnings))

  has_issue <- is_singular ||
    convergence_code != 0 ||
    length(captured_messages) > 0 ||
    length(captured_warnings) > 0

  return(list(
    is_singular = is_singular,
    optimizer = fit@optinfo$optimizer,
    convergence_code = convergence_code,
    lmer_messages = captured_messages,
    lmer_warnings = captured_warnings,
    has_issue = has_issue
  ))
}

.fit_lmer_with_diagnostics <- function(formula, dat_for_lmer) {
  captured_messages <- character()
  captured_warnings <- character()

  fit <- withCallingHandlers(
    lme4::lmer(
      formula,
      data = dat_for_lmer,
      REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
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
    diagnostics = .build_lmer_fit_diagnostics(
      fit,
      captured_messages = captured_messages,
      captured_warnings = captured_warnings
    )
  ))
}

.print_lmer_fit_diagnostics <- function(diagnostics) {
  if (!diagnostics$has_issue) {
    return(invisible(NULL))
  }

  cat("\n")
  cat("G-study fit diagnostics:\n")
  cat("  Optimizer: ", diagnostics$optimizer, "\n", sep = "")
  cat("  Convergence code: ", diagnostics$convergence_code, "\n", sep = "")
  cat("  Singular fit: ", diagnostics$is_singular, "\n", sep = "")

  if (length(diagnostics$lmer_messages) > 0) {
    cat("  lmer messages:\n")
    for (message in diagnostics$lmer_messages) {
      cat("    - ", message, "\n", sep = "")
    }
  }

  if (length(diagnostics$lmer_warnings) > 0) {
    cat("  lmer warnings:\n")
    for (warning in diagnostics$lmer_warnings) {
      cat("    - ", warning, "\n", sep = "")
    }
  }

  return(invisible(NULL))
}

# Variance components and their covariance from an lmer fit. The covariance is
# the inverse observed information of the REML deviance on the SD scale, mapped
# to the variance scale by the delta method.
.extract_vc_and_se <- function(fit) {
  vc_raw <- as.data.frame(lme4::VarCorr(fit))
  n_re <- nrow(vc_raw) - 1L

  vc <- c(vc_raw$vcov[seq_len(n_re)], sigma(fit)^2)
  names(vc) <- c(vc_raw$grp[seq_len(n_re)], "residual")
  sd_est <- c(vc_raw$sdcor[seq_len(n_re)], sigma(fit))

  # The profiled REML deviance function from lmer(devFunOnly = TRUE). Evaluating
  # it at theta updates the pp and resp objects in its environment, from which
  # the REML deviance with sigma as a free parameter is rebuilt (as in
  # lmerTest's devfun_vp()). Building a new deviance function leaves the fit
  # unchanged.
  profiled_devfun <- lme4::lmer(
    formula(fit),
    data = model.frame(fit),
    REML = TRUE,
    devFunOnly = TRUE
  )
  devfun_env <- environment(profiled_devfun)
  pp <- devfun_env$pp
  resp <- devfun_env$resp
  ld_w <- sum(log(resp$weights))
  n <- nrow(pp$X)
  p <- ncol(pp$X)

  # Every random term is a scalar intercept, so theta_k = sd_k / sigma, in the
  # order of the VarCorr rows.
  reml_sd_devfun <- function(pars) {
    sigma_hat <- pars[length(pars)]
    theta <- pars[seq_len(n_re)] / sigma_hat
    profiled_devfun(theta)

    sigma2_hat <- sigma_hat^2
    pwrss <- resp$wrss() + pp$sqrL(1)

    reml_deviance <- pp$ldL2() + pp$ldRX2() - ld_w +
      pwrss / sigma2_hat +
      (n - p) * log(2 * pi * sigma2_hat)
    return(reml_deviance)
  }

  vcov_var <- matrix(
    NA_real_,
    nrow = length(sd_est),
    ncol = length(sd_est),
    dimnames = list(names(vc), names(vc))
  )

  hessian_result <- try(
    {
      hessian <- optimHess(sd_est, reml_sd_devfun)
      hessian <- (hessian + t(hessian)) / 2
      eigenvalues <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigenvalues) <= 0) {
        stop("the observed information matrix is not positive definite", call. = FALSE)
      }
      if (min(eigenvalues) < max(eigenvalues) * 1e-12) {
        stop("the observed information matrix is numerically singular", call. = FALSE)
      }
      sd_vcov <- 2 * solve(hessian)
      jacobian <- diag(2 * sd_est)
      jacobian %*% sd_vcov %*% jacobian
    },
    silent = TRUE
  )

  if (inherits(hessian_result, "try-error")) {
    warning(
      "Could not obtain a valid observed-information covariance (",
      conditionMessage(attr(hessian_result, "condition")),
      "). Returning NA standard errors.",
      call. = FALSE
    )
  } else {
    vcov_var <- hessian_result
    dimnames(vcov_var) <- list(names(vc), names(vc))
  }

  # At a zero-boundary SD the delta-method Jacobian row is zero, which makes the
  # propagated SE 0; these SEs are set to NA.
  boundary_tol <- 1e-4
  is_boundary <- c(vc_raw$sdcor[seq_len(n_re)] < boundary_tol * sigma(fit), FALSE)
  names(is_boundary) <- names(vc)
  if (any(is_boundary)) {
    vcov_var[is_boundary, ] <- NA_real_
    vcov_var[, is_boundary] <- NA_real_
  }

  se_var <- sqrt(diag(vcov_var))
  names(se_var) <- names(vc)

  return(list(vc = vc, se = se_var, vcov = vcov_var, is_boundary = is_boundary))
}

# Rename lmer components to G-theory effect names, for example P:F1 to PxF1.
.extract_variance_components <- function(fit, design) {
  raw <- .extract_vc_and_se(fit)
  effect_names <- c(design$random_effects, design$residual_effect)
  lmer_names <- c(.lmer_effect_name(design$random_effects), "residual")

  boundary_effects <- effect_names[raw$is_boundary[lmer_names]]
  if (length(boundary_effects) > 0) {
    warning(
      "Variance component(s) estimated at the zero boundary: ",
      paste(boundary_effects, collapse = ", "),
      ". Delta-method standard errors are invalid at the boundary; returning NA.",
      call. = FALSE
    )
  }

  vc <- raw$vc[lmer_names]
  names(vc) <- effect_names

  vc_se <- raw$se[lmer_names]
  names(vc_se) <- effect_names

  vc_vcov <- raw$vcov[lmer_names, lmer_names, drop = FALSE]
  dimnames(vc_vcov) <- list(effect_names, effect_names)

  return(list(vc = vc, se = vc_se, vcov = vc_vcov))
}

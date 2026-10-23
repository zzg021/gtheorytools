#' gtheorytools: Tools for Generalizability Theory
#'
#' @description
#' Estimates variance components for one- and two-facet G studies with crossed
#' or nested designs, and computes D studies with random, finite, or fixed
#' facets. Reports standard errors, compares D studies across sample sizes,
#' and selects designs that meet a target coefficient.
#'
#' - [g_study()] fits a one- or two-facet G study.
#' - [d_study()] computes a D study for one set of sample sizes.
#' - [compare_d_studies()] compares D studies across combinations of facet sample sizes.
#' - [find_d_designs()] selects D studies that meet a G or Phi target.
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software*,
#' *67*(1), 1--48. \doi{10.18637/jss.v067.i01}
#'
#' @importFrom stats as.formula formula model.frame optimHess sigma
#' @keywords internal
"_PACKAGE"

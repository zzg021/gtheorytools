#' gtheorytools: Tools for Generalizability Theory
#'
#' @details
#' - [g_study()] estimates a one- or two-facet G study.
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
#' @seealso <https://zzg021.github.io/gtheorytools/>
#'
#' @importFrom methods new
#' @importFrom stats as.formula formula model.frame optimHess sigma
#' @keywords internal
"_PACKAGE"

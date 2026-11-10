#' gtheorytools: Tools for Generalizability Theory
#'
#' @details
#' - [g_study()] fits a one- or two-facet G study named by keyword, such as
#'   `design = "partially_nested"`, with persons as the object of measurement.
#' - [d_study()] computes a D study for one set of sample sizes.
#' - [compare_d_studies()] compares D studies across combinations of facet
#'   sample sizes.
#' - [g_study_general()], [d_study_general()], and
#'   [compare_d_studies_general()] do the same for any design, written as a
#'   design string such as `"p x (r:t) x o"`, and any nonnested random facet
#'   as the object of measurement.
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
#' @importFrom utils combn
#' @keywords internal
"_PACKAGE"

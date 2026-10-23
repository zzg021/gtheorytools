#' Brennan's Synthetic Data Set No. 1: Persons x Items
#'
#' Dichotomous scores of 10 persons on 12 items, a p x i design (Brennan,
#' 2001, Table 2.2). Brennan's G study is in Table 2.3 and his D studies in
#' Table 2.4.
#'
#' @format A data frame with 120 rows and 3 variables:
#' \describe{
#'   \item{p}{Person, 1 to 10.}
#'   \item{i}{Item, 1 to 12.}
#'   \item{score}{Item score, 0 or 1.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 2.2.
#'
#' @seealso [g_study()] with `design = "crossed"`.
#'
#' @examples
#' head(brennan_p_x_i)
"brennan_p_x_i"

#' Brennan's Synthetic Data Set No. 2: Items Nested Within Persons
#'
#' Scores from 0 to 9 of 10 persons on 8 items each, an i:p design in which
#' every person takes a different set of items (Brennan, 2001, Table 2.8).
#' Brennan's G and D studies are in Table 2.9. The scores are those of
#' [brennan_p_x_i_x_o], with items 1-4 on occasion 1 and items 5-8 on
#' occasion 2.
#'
#' @format A data frame with 80 rows and 3 variables:
#' \describe{
#'   \item{p}{Person, 1 to 10.}
#'   \item{i}{Item within person, 1 to 8; the same number denotes a different
#'     item for each person.}
#'   \item{score}{Score, 0 to 9.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 2.8.
#'
#' @seealso [g_study()] with `design = "nested"`.
#'
#' @examples
#' head(brennan_i_in_p)
"brennan_i_in_p"

#' Brennan's Synthetic Data Set No. 3: Persons x Items x Occasions
#'
#' Scores from 0 to 9 of 10 persons on 4 items on each of 2 occasions, a
#' p x i x o design (Brennan, 2001, Table 3.1). Brennan's G study is in
#' Table 3.3 and his D studies in Tables 4.4 (p x I x O) and 4.5 (p x (I:O)).
#'
#' @format A data frame with 80 rows and 4 variables:
#' \describe{
#'   \item{p}{Person, 1 to 10.}
#'   \item{i}{Item, 1 to 4.}
#'   \item{o}{Occasion, 1 or 2.}
#'   \item{score}{Score, 0 to 9.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 3.1.
#'
#' @seealso [g_study()] with `design = "fully_crossed"`.
#'
#' @examples
#' head(brennan_p_x_i_x_o)
"brennan_p_x_i_x_o"

#' Brennan's Synthetic Data Set No. 4: Persons x (Raters Nested Within Tasks)
#'
#' Ratings of 10 persons on 3 tasks, each task scored by its own 4 raters, a
#' p x (r:t) design (Brennan, 2001, Table 3.2). Brennan's G study is in
#' Table 3.4 and his D studies in Tables 4.6 (random model) and 4.8 (random
#' model and tasks fixed). Appendix F gives the GENOVA run with tasks fixed.
#'
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{p}{Person, 1 to 10.}
#'   \item{t}{Task, 1 to 3.}
#'   \item{r}{Rater, 1 to 12: raters 1-4 score task 1, 5-8 task 2, and 9-12
#'     task 3.}
#'   \item{score}{Rating, 1 to 9.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 3.2.
#'
#' @seealso [g_study()] with `design = "partially_nested"`.
#'
#' @examples
#' head(brennan_p_x_r_in_t)
"brennan_p_x_r_in_t"

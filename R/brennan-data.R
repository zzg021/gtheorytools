#' Brennan's Synthetic Data Set No. 1: Persons x Items
#'
#' Dichotomous scores of 10 persons on 12 items, a p x i design. Brennan's
#' (2001) G study is in Table 2.3 and his D studies in Table 2.4.
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
#' @seealso [g_study()]
#'
#' @examples
#' head(brennan_p_x_i)
"brennan_p_x_i"

#' Brennan's Synthetic Data Set No. 2: Items Nested Within Persons
#'
#' Scores from 0 to 9 of 10 persons on 8 items each, an i:p design.
#' Brennan's (2001) G and D studies are in Table 2.9. The scores are those of
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
#' @seealso [g_study()]
#'
#' @examples
#' head(brennan_i_in_p)
"brennan_i_in_p"

#' Brennan's Synthetic Data Set No. 3: Persons x Items x Occasions
#'
#' Scores from 0 to 9 of 10 persons on 4 items on each of 2 occasions, a
#' p x i x o design. Brennan's (2001) G study is in Table 3.3 and his
#' D studies in Tables 4.4 (p x I x O) and 4.5 (p x (I:O)).
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
#' @seealso [g_study()]
#'
#' @examples
#' head(brennan_p_x_i_x_o)
"brennan_p_x_i_x_o"

#' Brennan's Synthetic Data Set No. 4: Persons x (Raters Nested Within Tasks)
#'
#' Ratings of 10 persons on 3 tasks, each task scored by its own 4 raters, a
#' p x (r:t) design. Brennan's (2001) G study is in Table 3.4 and his
#' D studies in Tables 4.6 (random model) and 4.8 (random model and tasks
#' fixed). Appendix F gives the GENOVA run with tasks fixed.
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
#' @seealso [g_study()]
#'
#' @examples
#' head(brennan_p_x_r_in_t)
"brennan_p_x_r_in_t"

#' Brennan's Unbalanced Example: Items Nested Within Persons
#'
#' Scores from 2 to 9 of 9 persons on 3 to 5 items each, an unbalanced i:p
#' design. Brennan's (2001) G study is in Table 7.2 and his D study in
#' Section 7.2.1.
#'
#' @format A data frame with 37 rows and 3 variables:
#' \describe{
#'   \item{p}{Person, 1 to 9.}
#'   \item{i}{Item within person, 1 to 5; the same number denotes a different
#'     item for each person.}
#'   \item{score}{Score, 2 to 9.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 7.2.
#'
#' @seealso [g_study()], [d_study()]
#'
#' @examples
#' table(brennan_i_in_p_unbalanced$p)
"brennan_i_in_p_unbalanced"

#' Brennan's Unbalanced Example: Persons x (Items Nested Within h)
#'
#' Scores from 0 to 9 of 8 persons on 8 items, with 2, 4, and 2 items in the
#' three levels of h, an unbalanced p x (i:h) design. Brennan's (2001) G study
#' is in Table 7.5 and his D study in Section 7.2.2.
#'
#' @format A data frame with 64 rows and 4 variables:
#' \describe{
#'   \item{p}{Person, 1 to 8.}
#'   \item{h}{Level of h, 1 to 3.}
#'   \item{i}{Item, 1 to 8: items 1-2 are in level 1 of h, 3-6 in level 2, and
#'     7-8 in level 3.}
#'   \item{score}{Score, 0 to 9.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 7.5.
#'
#' @seealso [g_study()], [d_study()]
#'
#' @examples
#' head(brennan_p_x_i_in_h_unbalanced, 8)
"brennan_p_x_i_in_h_unbalanced"

#' Huynh's Example: Persons x Items With Missing Scores
#'
#' Dichotomous scores of 12 persons on 6 items, a p x i design in which each
#' person took 4 to 6 of the items. Brennan's (2001) G study is in Table 7.7
#' and his D study in Section 7.2.4.
#'
#' @format A data frame with 72 rows and 3 variables:
#' \describe{
#'   \item{p}{Person, 1 to 12.}
#'   \item{i}{Item, 1 to 6.}
#'   \item{score}{Item score, 0 or 1; `NA` for the 13 items a person did not
#'     take.}
#' }
#'
#' @source Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Table 7.7, from Huynh (1977).
#'
#' @seealso [g_study()], [d_study()]
#'
#' @examples
#' tapply(!is.na(brennan_p_x_i_missing$score), brennan_p_x_i_missing$p, sum)
"brennan_p_x_i_missing"

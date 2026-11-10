# gtheorytools

**G studies and D studies in R.**

gtheorytools estimates variance components for generalizability (G) studies
and uses them to plan decision (D) studies. It supports designs with any
number of crossed and nested facets, including random, finite, and fixed
facets, and reports standard errors for variance components and coefficients.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("zzg021/gtheorytools")
```

## Quick start

Start with Brennan's (2001) person-by-item data. Every person responds to the
same items, so persons are crossed with items. `p` and `f1` name the person
and item columns.

```r
library(gtheorytools)

g_result <- g_study(brennan_p_x_i, p = "p", f1 = "i", design = "crossed")
g_result$variance_components
```

Next, estimate the coefficients for a D study of persons with 12 items
sampled from an infinite universe of items. The sizes are named `F1`, the
facet given as `f1`.

```r
d_result <- d_study(g_result, universe_sizes = c(F1 = Inf), d_sizes = c(F1 = 12))
d_result$coefficient_results
```

The output includes universe-score variance (`tau`), relative error variance
(`delta`), absolute error variance (`Delta`), and their standard errors.
The generalizability coefficient (`G`) concerns relative decisions, such as
ranking persons. The dependability coefficient (`Phi`) concerns absolute
decisions, such as comparing a score with a cut score.

### Compare sample sizes

Use the same G study to compare item counts and select the candidate designs
whose estimated G coefficient is at least 0.80:

```r
comparison <- compare_d_studies(
  g_result,
  universe_sizes = c(F1 = Inf),
  d_sizes = list(F1 = c(4, 8, 12, 16))
)
comparison$scenarios

selection <- find_d_designs(comparison, coefficient = "G", target = 0.80)
selection$designs
```

## Designs

Supply data in long format with one row per score, a column for each facet,
and no replicated design cells. The score column is named `score` by
default; use `score_column` for a different name.

`g_study()` takes the one- and two-facet designs by keyword, with persons as
the object of measurement. `p`, `f1`, and `f2` name the columns, and P, F1,
and F2 stand for them in the component names. In the nested designs `f2` is
nested within `f1`.

| `design` | Arrangement | Variance components |
|:--|:--|:--|
| `"crossed"` | Persons crossed with F1 | P, F1, PxF1 |
| `"nested"` | F1 nested within persons | P, F1:P |
| `"fully_crossed"` | Persons crossed with F1 and F2 | P, F1, F2, PxF1, PxF2, F1xF2, PxF1xF2 |
| `"partially_nested"` | F2 nested within F1, crossed with persons | P, F1, PxF1, F2:F1, F2:F1xP |
| `"fully_nested"` | F2 nested within F1, nested within persons | P, F1:P, F2:F1:P |

`g_study()` fits linear mixed models by REML. `method = "anova"` gives the
ANOVA estimates from the expected mean squares, with urGENOVA's analogous
sums of squares for unbalanced data (Brennan, 2001, Sec. 7.1), and standard
errors equal to GENOVA's for balanced data (Eq. 6.1). A facet fixed in the
universe of admissible observations is named in `fixed`, such as
`fixed = "t"`; this requires complete, balanced data.

In `d_study()`, `universe_sizes` and `d_sizes` make each facet random,
finite, or fixed. `d_design` gives a more nested D-study design, with
`d_nesting` for the nesting order from a fully crossed G study: from
`"fully_crossed"` with items and occasions, `d_design = "partially_nested"`
and `d_nesting = "F1:F2"` give different items on each occasion. Without
`d_sizes`, `d_study()` repeats the G study, including the unequal numbers of
conditions and missing scores of unbalanced data (Sec. 7.2).

### Any design: design strings

`g_study_general()`, `d_study_general()`, and `compare_d_studies_general()`
take the design as a string of one-letter facet symbols, `x` for crossing,
`:` for nesting, and parentheses, with any number of facets. By default the
symbols are the data's column names; `facets` maps them to other columns, as
in `facets = c(p = "person", i = "item")`.

| `design` | Arrangement | Variance components |
|:--|:--|:--|
| `"p x i"` | Persons crossed with items | p, i, pi |
| `"i:p"` | Items nested within persons | p, i:p |
| `"p x i x o"` | Persons crossed with items and occasions | p, i, o, pi, po, io, pio |
| `"p x (r:t)"` | Raters nested within tasks, crossed with persons | p, t, r:t, pt, pr:t |
| `"r:t:p"` | Raters nested within tasks, nested within persons | p, t:p, r:t:p |
| `"p x (r:t) x o"` | Raters nested within tasks, crossed with persons and occasions | 11 components |

In `d_study_general()`, `object_of_measurement` names the facet whose levels
are measured, the sizes are named by facet symbol, and `d_design` is a
design string, such as `"p x (i:o)"` from a `"p x i x o"` G study:

```r
g_pio <- g_study_general(brennan_p_x_i_x_o, design = "p x i x o")
d_study_general(
  g_pio,
  object_of_measurement = "p",
  universe_sizes = c(i = Inf, o = Inf),
  d_sizes = c(i = 4, o = 2),
  d_design = "p x (i:o)"
)
```

## Example data

Seven synthetic data sets from Brennan (2001) are included, the last three
unbalanced.

| Data set | Design | Source |
|:--|:--|:--|
| `brennan_p_x_i` | `"p x i"` | Table 2.3 |
| `brennan_i_in_p` | `"i:p"` | Table 2.9 |
| `brennan_p_x_i_x_o` | `"p x i x o"` | Table 3.3 |
| `brennan_p_x_r_in_t` | `"p x (r:t)"` | Table 3.4 |
| `brennan_i_in_p_unbalanced` | `"i:p"`, 3 to 5 items per person | Table 7.2 |
| `brennan_p_x_i_in_h_unbalanced` | `"p x (i:h)"`, 2, 4, and 2 items | Table 7.5 |
| `brennan_p_x_i_missing` | `"p x i"`, 13 missing scores | Table 7.7 |

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*, *67*(1),
1–48. <https://doi.org/10.18637/jss.v067.i01>.

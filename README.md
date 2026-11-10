# gtheorytools

**G studies and D studies in R.**

gtheorytools estimates variance components for generalizability (G) studies
and uses them to plan decision (D) studies. It supports one- and two-facet
designs with crossed or nested facets, including random, finite, and fixed
facets, and reports standard errors for variance components and coefficients.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("zzg021/gtheorytools")
```

## Quick start

Start with Brennan's (2001) person-by-item data. Every person responds to the
same items, so the design is crossed.

```r
library(gtheorytools)

g_result <- g_study(
  brennan_p_x_i,
  p = "p",
  f1 = "i",
  design = "crossed"
)
g_result$variance_components
```

Next, estimate the coefficients for a D study with 12 items sampled from an
infinite universe of items. `F1` refers to the column supplied as `f1` in
the G study; in a two-facet study, `F2` refers to `f2`.

```r
d_result <- d_study(
  g_result,
  universe_sizes = c(F1 = Inf),
  d_sizes = c(F1 = 12)
)
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

## Supported designs

Supply data in long format with one row per score, columns identifying the
person and facets, and no replicated design cells. The score column is named
`score` by default; use `score_column` for a different name.

| Facets | `design` | Arrangement |
|:--|:--|:--|
| One | `"crossed"` | Persons crossed with a facet |
| One | `"nested"` | Facet conditions nested within persons |
| Two | `"fully_crossed"` | Persons crossed with both facets |
| Two | `"partially_nested"` | `f2` nested within `f1`, crossed with persons |
| Two | `"fully_nested"` | `f2` nested within `f1`, nested within persons |

`g_study()` fits linear mixed models by REML. `method = "anova"` gives the
ANOVA estimates from the expected mean squares, with urGENOVA's analogous
sums of squares for unbalanced data (Brennan, 2001, Sec. 7.1). Its standard
errors are those GENOVA reports for balanced data (Eq. 6.1) and, for
unbalanced data, unbiased normal-theory estimates that reduce to them. For
complete, balanced data, `fixed` fixes a facet of a fully crossed or
partially nested G study. In a D study, `universe_sizes` and `d_sizes`
specify random, finite, or fixed facets, and `d_design` can be more nested
than the G-study design. Without `d_sizes`, `d_study()` repeats the G study,
including the unequal numbers of conditions and missing scores of unbalanced
data (Sec. 7.2).

## Example data

Seven synthetic data sets from Brennan (2001) are included, the last three
unbalanced.

| Data set | Design | Source |
|:--|:--|:--|
| `brennan_p_x_i` | Persons crossed with items | Table 2.3 |
| `brennan_i_in_p` | Items nested within persons | Table 2.9 |
| `brennan_p_x_i_x_o` | Persons crossed with items and occasions | Table 3.3 |
| `brennan_p_x_r_in_t` | Raters nested within tasks, crossed with persons | Table 3.4 |
| `brennan_i_in_p_unbalanced` | Items nested within persons, 3 to 5 per person | Table 7.2 |
| `brennan_p_x_i_in_h_unbalanced` | Items nested within h, crossed with persons; 2, 4, and 2 items | Table 7.5 |
| `brennan_p_x_i_missing` | Persons crossed with items, 13 missing scores | Table 7.7 |

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*, *67*(1),
1–48. <https://doi.org/10.18637/jss.v067.i01>.

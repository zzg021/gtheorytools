# gtheorytools

**G studies and D studies in R.**

gtheorytools estimates variance components for generalizability (G) studies
and uses them to plan decision (D) studies. It supports one- and two-facet
designs with crossed or nested facets, including random, finite, and fixed
facets, and reports standard errors for variance components and coefficients.

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

Selection uses point estimates. Each comparison also reports standard errors
and any warnings or errors for the candidate D studies.

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

`g_study()` fits linear mixed models by REML. A facet can be fixed in the
G study for supported two-facet designs using `fixed`; this requires complete,
balanced data. In a D study, `universe_sizes` and `d_sizes` specify random,
finite, or fixed facets. D studies can also use a more nested design than
the G study. See `g_study()` and `d_study()` for the supported combinations
and definitions.

## Example data

Four synthetic data sets from Brennan (2001) are included. Their help pages
describe the variables and give examples.

| Data set | Design | Source |
|:--|:--|:--|
| `brennan_p_x_i` | Persons crossed with items | Table 2.3 |
| `brennan_i_in_p` | Items nested within persons | Table 2.9 |
| `brennan_p_x_i_x_o` | Persons crossed with items and occasions | Table 3.3 |
| `brennan_p_x_r_in_t` | Raters nested within tasks, crossed with persons | Table 3.4 |

## Learn more

The function reference includes arguments, returned values, and worked
examples for `g_study()`, `d_study()`, `compare_d_studies()`, and
`find_d_designs()`. In R, use `?g_study` to open the corresponding help page.

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*, *67*(1),
1–48. <https://doi.org/10.18637/jss.v067.i01>.

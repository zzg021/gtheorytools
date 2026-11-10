# urGENOVA's p x (i:h) run with empty cells, and the data columns of each
# effect, in the package's order with the residual last, for the calculations
# from n x n matrices.

# Long-format scores from a persons x conditions matrix, NA for a missing
# score; `conditions` gives the facet levels of each column.
long_scores <- function(scores, conditions) {
  result <- data.frame(
    p = rep(seq_len(nrow(scores)), times = ncol(scores)),
    conditions[rep(seq_len(ncol(scores)), each = nrow(scores)), , drop = FALSE],
    score = as.vector(scores),
    row.names = NULL
  )
  return(result)
}

# urGENOVA manual, Run 6: p x (i:h) with 3, 6, and 3 items in the three levels
# of h and six empty cells.
urgenova_p_x_i_in_h <- function() {
  scores <- rbind(
    c(5, 6, 5, 5, 5, 3, 4, NA, 6, 7, 3, 3), c(7, 1, 5, NA, 5, 3, 3, 3, 5, 5, 3, 0),
    c(3, 4, 3, 3, 5, 3, 3, 5, 6, NA, 1, 6), c(7, 5, 5, 3, 3, 1, 4, 3, 5, 3, 3, 5),
    c(9, 2, 9, 7, 7, NA, 3, 7, 2, 7, 5, 3), c(3, 4, 3, 5, 3, 3, 6, 3, 4, 5, 1, 2),
    c(7, 3, NA, 7, 7, 5, 5, 7, 5, 5, 5, 4), c(5, 8, 5, 7, 7, 5, 5, NA, 3, 2, 1, 1),
    c(9, 9, 8, 8, 6, 6, 6, 5, 5, 8, 1, 1), c(4, 4, 4, 3, 3, 5, 6, 5, 5, 7, 1, 1)
  )
  conditions <- data.frame(h = rep(1:3, c(3, 6, 3)), i = c(1:3, 1:6, 1:3))
  return(long_scores(scores, conditions))
}

effects_i_in_p <- list(p = "p", "i:p" = c("i", "p"))
effects_p_x_i <- list(p = "p", i = "i", pi = c("p", "i"))
effects_p_x_i_in_h <- list(
  p = "p", h = "h", "i:h" = c("i", "h"), ph = c("p", "h"), "pi:h" = c("p", "i", "h")
)

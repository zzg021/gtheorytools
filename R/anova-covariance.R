# Number of scores in each pair of cells of two facet sets, as a sparse matrix
# with the cells of the first set as rows: N_xy = C_x' C_y, where C_x marks the
# x cell of each score.
.count_table <- function(x_index, y_index, n_x, n_y) {
  pair <- (as.numeric(x_index) - 1) * n_y + y_index
  pairs <- unique(pair)
  counts <- tabulate(match(pair, pairs))
  result <- Matrix::sparseMatrix(
    i = (pairs - 1) %/% n_y + 1,
    j = (pairs - 1) %% n_y + 1,
    x = counts,
    dims = c(n_x, n_y)
  )
  return(result)
}

# Column-major positions, from 0, of the stored entries of a sparse matrix, in
# increasing order.
.entry_positions <- function(m) {
  column <- rep.int(as.numeric(seq_len(ncol(m)) - 1), diff(m@p))
  result <- m@i + column * nrow(m)
  return(result)
}

# Inner products <A_a, A_b> of sparse matrices of equal dimensions, as the
# cross-products of a sparse matrix with one column per A_a and one row per
# position that holds an entry.
.gram_matrix <- function(matrices) {
  n_matrices <- length(matrices)
  positions <- vector("list", n_matrices)
  values <- vector("list", n_matrices)
  for (a in seq_len(n_matrices)) {
    positions[[a]] <- .entry_positions(matrices[[a]])
    values[[a]] <- matrices[[a]]@x
  }
  all_positions <- unlist(positions)
  distinct_positions <- sort(unique(all_positions))
  stacked <- new(
    "dgCMatrix",
    i = findInterval(all_positions, distinct_positions) - 1L,
    p = c(0L, cumsum(lengths(positions))),
    x = unlist(values),
    Dim = c(length(distinct_positions), n_matrices)
  )
  result <- as.matrix(Matrix::crossprod(stacked))
  return(result)
}

# D_x^{-1/2} N_xy D_y^{-1/2} for the count table N_xy and the cell counts of x
# and y.
.scaled_count_table <- function(table, x_counts, y_counts) {
  y_count <- rep.int(y_counts, diff(table@p))
  table@x <- table@x / (sqrt(x_counts[table@i + 1L]) * sqrt(y_count))
  return(table)
}

# Sums of the entrywise products of a and b over consecutive blocks of columns
# of the given widths, for a and b of equal dimensions with every nonzero
# entry of b nonzero in a.
.block_inner_products <- function(a, b, block_widths) {
  position_a <- .entry_positions(a)
  position_b <- .entry_positions(b)
  products <- b@x * a@x[findInterval(position_b, position_a)]
  entries_per_block <- diff(b@p[c(0, cumsum(block_widths)) + 1])
  block <- rep.int(seq_along(block_widths), entries_per_block)
  result <- as.vector(rowsum(products, block, reorder = FALSE))
  return(result)
}

# tr(P_c Z_k Z_k' P_d Z_l Z_l') for all T terms c, d and effects k, l, with
# P_c the projection on the c cell means and Z_k the indicator of the k cell of
# each score. With N_xy the count table, D_x the diagonal of the x cell counts,
# and R_xk = D_x^{-1/2} N_xk, the trace is tr(R_ck R_dk' R_dl R_cl'), the
# inner product of R_ck R_dk' and R_cl R_dl', and also of R_cl' R_ck and
# R_dl' R_dk. For a pair of T terms, the traces over all pairs of components
# are therefore one Gram matrix of the R_ck R_dk'; for a pair of components,
# the traces over a set of T terms are one Gram matrix of the R_cl' R_ck. The
# second form suits T terms with many small cells. The residual has P = I and
# Z = I, so a trace that contains it is a shorter product of count tables.
.anova_traces <- function(dat, effect_names, facet_sets) {
  n_scores <- nrow(dat)
  n_effects <- length(effect_names)
  residual <- effect_names[[n_effects]]
  components <- effect_names[-n_effects]
  n_components <- n_effects - 1
  t_sets <- c("mean", effect_names)
  # The T terms other than the residual.
  partitions <- c("mean", components)
  n_partitions <- n_effects

  cell_index <- list(mean = rep(1L, n_scores))
  for (effect in components) {
    cell_index[[effect]] <- .cell_index(dat[facet_sets[[effect]]])
  }
  cell_counts <- lapply(cell_index, tabulate)
  n_cells <- lengths(cell_counts)
  component_cells <- n_cells[components]

  # N_xy for every pair of these T terms and R_xk for every component k, with
  # the numbers of nonzero entries in each column and each row of N_xk.
  tables <- list()
  scaled_tables <- list()
  x_per_k_cell <- list()
  k_per_x_cell <- list()
  for (k in components) {
    x_per_k_cell[[k]] <- matrix(0, nrow = n_cells[[k]], ncol = n_partitions)
  }
  for (j in seq_len(n_partitions)) {
    x <- partitions[[j]]
    tables[[x]] <- list()
    scaled_tables[[x]] <- list()
    k_per_x_cell[[x]] <- matrix(0, nrow = n_cells[[x]], ncol = n_components)
    for (y in partitions) {
      tables[[x]][[y]] <- .count_table(cell_index[[x]], cell_index[[y]], n_cells[[x]], n_cells[[y]])
    }
    for (a in seq_len(n_components)) {
      k <- components[[a]]
      table_xk <- tables[[x]][[k]]
      x_per_k_cell[[k]][, j] <- diff(table_xk@p)
      k_per_x_cell[[x]][, a] <- tabulate(table_xk@i + 1L, nbins = n_cells[[x]])
      table_xk@x <- table_xk@x / sqrt(cell_counts[[x]][table_xk@i + 1L])
      scaled_tables[[x]][[k]] <- table_xk
    }
  }
  # R_xk and N_kl side by side over the components.
  scaled_by_component <- list()
  for (x in partitions) {
    scaled_by_component[[x]] <- do.call(cbind, unname(scaled_tables[[x]][components]))
  }
  tables_by_component <- list()
  for (k in components) {
    tables_by_component[[k]] <- do.call(cbind, unname(tables[[k]][components]))
  }

  # Multiplications of nonzero entries in R_ck R_dk' over all components, for
  # each pair of T terms, and in R_xl' R_xk over all pairs of components, for
  # each T term.
  t_pair_cost <- matrix(0, nrow = n_partitions, ncol = n_partitions, dimnames = list(partitions, partitions))
  for (k in components) {
    t_pair_cost <- t_pair_cost + crossprod(x_per_k_cell[[k]])
  }
  component_pair_cost <- numeric(n_partitions)
  names(component_pair_cost) <- partitions
  for (x in partitions) {
    products <- crossprod(k_per_x_cell[[x]])
    component_pair_cost[[x]] <- sum(products[upper.tri(products, diag = TRUE)])
  }

  # T terms whose traces among themselves are computed by pairs of
  # components: in decreasing order of cells, a T term joins when it saves more
  # multiplications in the products for pairs of T terms than it adds in those
  # for pairs of components.
  by_components <- character(0)
  for (x in partitions[order(n_cells, decreasing = TRUE)]) {
    saved <- sum(t_pair_cost[x, c(by_components, x)])
    if (saved > component_pair_cost[[x]]) {
      by_components <- c(by_components, x)
    }
  }

  # The traces are unchanged when (c, k) and (d, l) are exchanged or when k
  # and l are exchanged.
  traces <- array(
    NA_real_,
    dim = c(n_effects + 1, n_effects, n_effects + 1, n_effects),
    dimnames = list(t_sets, effect_names, t_sets, effect_names)
  )
  # Traces for pairs of T terms, with the residual as a component entering as
  # S_cd = D_c^{-1/2} N_cd D_d^{-1/2}.
  for (i in seq_len(n_partitions)) {
    for (j in i:n_partitions) {
      c <- partitions[[i]]
      d <- partitions[[j]]
      if (c %in% by_components && d %in% by_components) {
        next
      }
      products <- vector("list", n_effects)
      for (a in seq_len(n_components)) {
        k <- components[[a]]
        products[[a]] <- Matrix::tcrossprod(scaled_tables[[c]][[k]], scaled_tables[[d]][[k]])
      }
      products[[n_effects]] <- .scaled_count_table(tables[[c]][[d]], cell_counts[[c]], cell_counts[[d]])
      gram <- .gram_matrix(products)
      traces[c, , d, ] <- gram
      traces[d, , c, ] <- gram
    }
  }
  # Traces for pairs of components, with the residual as a T term entering as
  # R_res,l' R_res,k = N_lk.
  x_sets <- c(by_components, residual)
  for (a in seq_len(n_components)) {
    for (b in a:n_components) {
      k <- components[[a]]
      l <- components[[b]]
      products <- vector("list", length(x_sets))
      for (j in seq_along(by_components)) {
        x <- by_components[[j]]
        products[[j]] <- Matrix::crossprod(scaled_tables[[x]][[l]], scaled_tables[[x]][[k]])
      }
      products[[length(x_sets)]] <- tables[[l]][[k]]
      gram <- .gram_matrix(products)
      traces[x_sets, k, x_sets, l] <- gram
      traces[x_sets, l, x_sets, k] <- gram
    }
  }

  # tr(P_c Z_k Z_k' P_d) = <R_ck, S_cd R_dk> and tr(P_c P_d) = ||S_cd||^2 for
  # c and d both computed by pairs of components.
  for (i in seq_along(by_components)) {
    for (j in i:length(by_components)) {
      c <- by_components[[i]]
      d <- by_components[[j]]
      s_cd <- .scaled_count_table(tables[[c]][[d]], cell_counts[[c]], cell_counts[[d]])
      values <- .block_inner_products(
        s_cd %*% scaled_by_component[[d]], scaled_by_component[[c]], component_cells
      )
      traces[c, components, d, residual] <- values
      traces[c, residual, d, components] <- values
      traces[d, components, c, residual] <- values
      traces[d, residual, c, components] <- values
      traces[c, residual, d, residual] <- sum(s_cd@x^2)
      traces[d, residual, c, residual] <- sum(s_cd@x^2)
    }
  }
  # tr(Z_k Z_k' P_d Z_l Z_l') = <R_dl, R_dk N_kl> for the other T terms d.
  for (d in setdiff(partitions, by_components)) {
    for (k in components) {
      values <- .block_inner_products(
        scaled_tables[[d]][[k]] %*% tables_by_component[[k]], scaled_by_component[[d]], component_cells
      )
      traces[residual, k, d, components] <- values
      traces[d, components, residual, k] <- values
    }
  }
  # tr(P_d Z_l Z_l') = ||R_dl||^2, tr(P_d) = the number of d cells, and
  # tr(Z_l Z_l') = the number of scores.
  for (d in partitions) {
    values <- numeric(n_components)
    for (b in seq_len(n_components)) {
      values[[b]] <- sum(scaled_tables[[d]][[components[[b]]]]@x^2)
    }
    traces[residual, residual, d, components] <- values
    traces[d, components, residual, residual] <- values
    traces[residual, components, d, residual] <- values
    traces[d, residual, residual, components] <- values
    traces[residual, residual, d, residual] <- n_cells[[d]]
    traces[d, residual, residual, residual] <- n_cells[[d]]
  }
  traces[residual, residual, residual, ] <- n_scores
  traces[residual, , residual, residual] <- n_scores
  return(traces)
}

# Unbiased estimate, under normality, of the covariance of the ANOVA estimates
# for unbalanced data. Each T term is a quadratic form y'P_c y, with P_c the
# projection on the c cell means, and SS(a) = sum_c L[a, c] T(c). With
# V = sum_k sigma^2_k Z_k Z_k', where Z_k marks the k cell of each score,
# Cov(SS(a), SS(b)) = 2 tr(Q_a V Q_b V) under normality for
# Q_a = sum_c L[a, c] P_c; the mean drops out because each SS is unchanged by
# a constant added to every score. The covariance of the estimates, W SS with
# W = K^{-1} diag(1 / df), is therefore sum_kl G[ij, kl] sigma^2_k sigma^2_l.
# Because E[theta_k theta_l] = sigma^2_k sigma^2_l + Cov(theta_k, theta_l), the
# products sigma^2_k sigma^2_l are estimated without bias by solving
# (I + G) p = theta theta', and G p estimates the covariance. For balanced
# data this is Brennan's (2001) Eq. 6.1, with df + 2.
.anova_covariance_unbalanced <- function(dat, effect_names, df, ems_coefficients, vc) {
  n_effects <- length(effect_names)
  facet_sets <- lapply(effect_names, .parse_effect)
  names(facet_sets) <- effect_names
  lower_effects <- .lower_effects(facet_sets)
  traces <- .anova_traces(dat, effect_names, facet_sets)

  # Coefficients of the T terms in each sum of squares.
  t_sets <- c("mean", effect_names)
  ss_weights <- matrix(0, nrow = n_effects, ncol = n_effects + 1, dimnames = list(effect_names, t_sets))
  for (effect in effect_names[order(lengths(facet_sets))]) {
    ss_weights[effect, effect] <- 1
    ss_weights[effect, "mean"] <- -1
    ss_weights[effect, ] <- ss_weights[effect, ] -
      colSums(ss_weights[lower_effects[[effect]], , drop = FALSE])
  }

  estimate_weights <- solve(ems_coefficients) %*% diag(1 / df, nrow = n_effects)
  g_matrix <- matrix(0, nrow = n_effects^2, ncol = n_effects^2)
  for (k in seq_len(n_effects)) {
    for (l in seq_len(n_effects)) {
      ss_covariance <- ss_weights %*% traces[, k, , l] %*% t(ss_weights)
      g_matrix[, (l - 1) * n_effects + k] <- as.vector(
        2 * estimate_weights %*% ss_covariance %*% t(estimate_weights)
      )
    }
  }

  products <- solve(diag(n_effects^2) + g_matrix, as.vector(vc %o% vc))
  vcov <- matrix(g_matrix %*% products, nrow = n_effects, ncol = n_effects)
  vcov <- (vcov + t(vcov)) / 2
  return(vcov)
}

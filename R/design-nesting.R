# A design is its nesting relation: for each facet, the facets it is nested
# within, transitively closed. Facets not related by nesting are crossed.

# Design strings for the design keywords of g_study(), written with the
# symbols of Brennan's (2001) data sets.
.design_keyword_strings <- c(
  crossed = "p x i",
  nested = "i:p",
  fully_crossed = "p x i x o",
  partially_nested = "p x (r:t)",
  fully_nested = "r:t:p"
)

# Read a design string such as "p x (r:t)": one-letter facet symbols, "x" or
# "X" for crossing, ":" for nesting, and parentheses. In A:B every facet of A
# is nested within every facet of B, and the chain A:B:C also nests B within
# C. A symbol may appear more than once, as b in "p x (a:b) x (c:d:b)"; every
# appearance is the same facet. Symbols are letters other than x and X, and
# their case is kept, so p and P are different facets. Returns the facets in
# order of first appearance, their nesting relation, and the string with
# " x " between crossed terms.
.parse_design <- function(design, argument_name) {
  .check_single_string(design, argument_name)
  if (design %in% names(.design_keyword_strings)) {
    stop(
      "\"", design, "\" is a design keyword of g_study() and d_study(). Here ",
      argument_name, " is a design string, such as \"",
      .design_keyword_strings[[design]], "\" with your facet symbols.",
      call. = FALSE
    )
  }

  compact <- gsub("[[:space:]]", "", design)
  compact <- gsub("X", "x", compact, fixed = TRUE)
  characters <- strsplit(compact, "")[[1]]
  position <- 1L
  nesting_pairs <- list()

  fail <- function() {
    stop(
      "Could not read ", argument_name, " = \"", design, "\". A design string ",
      "uses one-letter facet symbols, \"x\" for crossing, \":\" for nesting, ",
      "and parentheses, as in \"p x (r:t)\".",
      call. = FALSE
    )
  }
  next_character <- function() {
    if (position > length(characters)) {
      return("")
    }
    return(characters[[position]])
  }

  # A facet symbol, or a parenthesized design.
  read_factor <- function() {
    token <- next_character()
    if (token == "(") {
      position <<- position + 1L
      facets <- read_crossing()
      if (next_character() != ")") {
        fail()
      }
      position <<- position + 1L
      return(facets)
    }
    if (grepl("^[A-Za-z]$", token) && token != "x") {
      position <<- position + 1L
      return(token)
    }
    fail()
  }

  # A:B:C, each facet of a level nested within each facet of the next.
  read_chain <- function() {
    levels <- list(read_factor())
    while (next_character() == ":") {
      position <<- position + 1L
      levels[[length(levels) + 1L]] <- read_factor()
    }
    for (k in seq_len(length(levels) - 1L)) {
      for (child in levels[[k]]) {
        for (parent in levels[[k + 1L]]) {
          nesting_pairs[[length(nesting_pairs) + 1L]] <<- c(child, parent)
        }
      }
    }
    return(unique(unlist(levels)))
  }

  read_crossing <- function() {
    facets <- read_chain()
    while (next_character() == "x") {
      position <<- position + 1L
      facets <- union(facets, read_chain())
    }
    return(facets)
  }

  facets <- read_crossing()
  if (position <= length(characters)) {
    fail()
  }

  parents <- rep(list(character(0)), length(facets))
  names(parents) <- facets
  for (pair in nesting_pairs) {
    parents[[pair[[1]]]] <- union(parents[[pair[[1]]]], pair[[2]])
  }
  repeat {
    changed <- FALSE
    for (facet in facets) {
      closed <- union(parents[[facet]], unlist(parents[parents[[facet]]]))
      if (length(closed) > length(parents[[facet]])) {
        parents[[facet]] <- closed
        changed <- TRUE
      }
    }
    if (!changed) {
      break
    }
  }
  for (facet in facets) {
    if (facet %in% parents[[facet]]) {
      stop(
        argument_name, " = \"", design, "\" nests ", facet,
        " within itself through a chain of nested facets.",
        call. = FALSE
      )
    }
    parents[[facet]] <- facets[facets %in% parents[[facet]]]
  }

  return(list(
    facets = facets,
    parents = parents,
    design = gsub("x", " x ", compact, fixed = TRUE)
  ))
}

# The nesting indices of an effect are the facets its members are nested
# within; the other facets are its primary indices (Brennan 2001, p. 62).
.nesting_indices <- function(members, parents) {
  return(members[members %in% unlist(parents[members])])
}

.primary_indices <- function(members, parents) {
  return(setdiff(members, .nesting_indices(members, parents)))
}

# A primary index is at level 0, and a nesting index one level above the
# highest member nested within it: in ac:d:b, d is at level 1 and b at
# level 2, through c:d:b.
.nesting_levels <- function(members, parents) {
  level <- integer(length(members))
  names(level) <- members
  children <- lapply(members, function(facet) {
    return(members[vapply(parents[members], function(p) facet %in% p, logical(1))])
  })
  names(children) <- members

  for (step in seq_along(members)) {
    for (facet in members) {
      if (length(children[[facet]]) > 0) {
        level[[facet]] <- 1L + max(level[children[[facet]]])
      }
    }
  }
  return(level)
}

# Brennan's (2001, pp. 61-62) name of an effect: its primary indices, then
# its nesting indices of each level after a colon, each level in the order of
# the facets (pr:t, ac:d:b).
.effect_name <- function(members, parents) {
  level <- .nesting_levels(members, parents)
  ordered <- names(parents)[names(parents) %in% members]
  level_names <- character(0)
  for (k in seq(0L, max(level))) {
    level_names <- c(level_names, paste0(ordered[level[ordered] == k], collapse = ""))
  }
  result <- paste(level_names, collapse = ":")
  return(result)
}

# The effects of a design are the sets of facets that contain the parents of
# each member (Brennan 2001, pp. 61-62). They are ordered as in Brennan's
# tables: by number of primary indices, then by number of facets, with the
# residual, which contains every facet, last.
.design_effects <- function(parents) {
  facets <- names(parents)
  result <- list()

  for (n_members in seq_along(facets)) {
    for (members in combn(facets, n_members, simplify = FALSE)) {
      if (all(unlist(parents[members]) %in% members)) {
        result[[length(result) + 1]] <- members
      }
    }
  }

  n_primary <- vapply(
    result,
    function(members) length(.primary_indices(members, parents)),
    integer(1)
  )
  n_members <- vapply(result, length, integer(1))
  is_residual <- n_members == length(facets)
  result <- result[order(is_residual, n_primary, n_members)]

  names(result) <- vapply(
    result,
    .effect_name,
    character(1),
    parents = parents
  )
  return(result)
}

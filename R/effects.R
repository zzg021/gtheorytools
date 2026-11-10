# The facets of an effect, sorted, from its name: pr:t gives p, r, t.
.parse_effect <- function(effect_name) {
  facets <- strsplit(gsub(":", "", effect_name, fixed = TRUE), "")[[1]]
  return(sort(unique(facets)))
}

# The effects whose facets are all among `facets`.
.effects_within <- function(effect_names, facets) {
  is_within <- vapply(
    effect_names,
    function(effect_name) all(.parse_effect(effect_name) %in% facets),
    logical(1)
  )
  return(effect_names[is_within])
}

# An effect name written with the data's column names, with " x " between the
# indices of a level: pr:t as "person x rater:task".
.effect_label <- function(effect_name, facet_columns) {
  levels <- strsplit(effect_name, ":", fixed = TRUE)[[1]]
  level_labels <- vapply(
    levels,
    function(level) {
      symbols <- strsplit(level, "")[[1]]
      return(paste(unname(facet_columns[symbols]), collapse = " x "))
    },
    character(1)
  )
  return(paste(level_labels, collapse = ":"))
}

# How a facet is shown in printed output: its symbol, with the data column in
# parentheses when the two differ.
.facet_display <- function(facet, facet_columns) {
  column <- facet_columns[[facet]]
  if (column == facet) {
    return(facet)
  }
  return(paste0(facet, " (", column, ")"))
}

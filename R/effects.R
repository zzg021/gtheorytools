.effect <- function(...) {
  return(sort(unique(c(...))))
}

# The facets of an effect, sorted, ignoring crossing and nesting.
.parse_effect <- function(effect_name) {
  facets <- unlist(strsplit(effect_name, "[:x]", perl = TRUE))
  return(.effect(facets[facets != ""]))
}

# Convert a G-theory effect name to the interaction syntax used by lme4.
.lmer_effect_name <- function(effect_name) {
  return(gsub("x", ":", effect_name, fixed = TRUE))
}

# An effect name written with the data's column names: PxF1 as "p x i".
.effect_label <- function(effect_name, facet_columns) {
  token_match <- gregexpr("F1|F2|P|[:x]", effect_name, perl = TRUE)
  tokens <- regmatches(effect_name, token_match)[[1]]
  label_parts <- character(length(tokens))

  for (i in seq_along(tokens)) {
    token <- tokens[[i]]

    if (token %in% names(facet_columns)) {
      label_parts[[i]] <- unname(facet_columns[[token]])
    } else if (token == "x") {
      label_parts[[i]] <- " x "
    } else {
      label_parts[[i]] <- ":"
    }
  }

  return(paste0(label_parts, collapse = ""))
}

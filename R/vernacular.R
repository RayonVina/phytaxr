# R/vernacular.R
# CRUD interface for the internal vernacular-to-taxon dictionary.
# The underlying object (.VERNACULAR_TO_TAXON) lives in R/sysdata.rda and is
# loaded as a locked binding by R -- it cannot be modified with <<-.
# Session-level mutations are stored in .vern_env$dict (a mutable environment).
# Changes are lost on restart. To make them permanent, edit data-raw/vernacular.R
# and re-run it.

# Private mutable environment, initialised once when the package loads.
.vern_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .vern_env$dict <- .VERNACULAR_TO_TAXON
}

# Internal accessor used by all functions below.
.vdict <- function() .vern_env$dict

#' List all vernacular-to-taxon mappings
#'
#' @return A named character vector (name = common name, value = taxon).
#' @export
vernacular_list <- function() {
  .vdict()
}

#' Look up a vernacular name
#'
#' @param name Character. Common name to look up (case-insensitive).
#' @return The mapped taxon name, or `NA_character_` if not found.
#' @export
vernacular_lookup <- function(name) {
  unname(.vdict()[tolower(trimws(name))])
}

#' Add a vernacular mapping (session only)
#'
#' Adds a new entry to the in-memory dictionary for the current session.
#' To make the change permanent, also add it to `data-raw/vernacular.R`
#' and re-run that script.
#'
#' @param common Character. Common name (stored as lowercase).
#' @param taxon  Character. Valid taxonomic name to map to.
#' @return The updated dictionary, invisibly.
#' @export
vernacular_add <- function(common, taxon) {
  key <- tolower(trimws(common))
  d   <- .vdict()
  if (key %in% names(d)) {
    message(sprintf(
      "'%s' already exists -> '%s'. Use vernacular_update() to change it.",
      key, d[[key]]
    ))
    return(invisible(d))
  }
  d[[key]]        <- taxon
  .vern_env$dict  <- d
  message(sprintf("Added: '%s' -> '%s'", key, taxon))
  invisible(.vdict())
}

#' Update an existing vernacular mapping (session only)
#'
#' @param common Character. Common name to update (case-insensitive).
#' @param taxon  Character. New taxonomic name.
#' @return The updated dictionary, invisibly.
#' @export
vernacular_update <- function(common, taxon) {
  key <- tolower(trimws(common))
  d   <- .vdict()
  if (!key %in% names(d)) {
    message(sprintf("'%s' not found. Use vernacular_add() to create it.", key))
    return(invisible(d))
  }
  old            <- d[[key]]
  d[[key]]       <- taxon
  .vern_env$dict <- d
  message(sprintf("Updated: '%s': '%s' -> '%s'", key, old, taxon))
  invisible(.vdict())
}

#' Remove a vernacular mapping (session only)
#'
#' @param common Character. Common name to remove (case-insensitive).
#' @return The updated dictionary, invisibly.
#' @export
vernacular_remove <- function(common) {
  key <- tolower(trimws(common))
  d   <- .vdict()
  if (!key %in% names(d)) {
    message(sprintf("'%s' not found in dictionary.", key))
    return(invisible(d))
  }
  .vern_env$dict <- d[names(d) != key]
  message(sprintf("Removed: '%s'", key))
  invisible(.vdict())
}

#' Search vernacular mappings by pattern
#'
#' @param pattern Character. Regex pattern to search in common names or
#'   taxon values.
#' @param where One of `"keys"` (search common names), `"values"` (search
#'   taxon names), or `"both"`. Default `"both"`.
#' @return A named character vector of matching entries.
#' @export
vernacular_search <- function(pattern, where = c("both", "keys", "values")) {
  where <- match.arg(where)
  d     <- .vdict()
  hits  <- switch(
    where,
    keys   = grepl(pattern, names(d), ignore.case = TRUE),
    values = grepl(pattern, d,        ignore.case = TRUE),
    both   = grepl(pattern, names(d), ignore.case = TRUE) |
             grepl(pattern, d,        ignore.case = TRUE)
  )
  d[hits]
}

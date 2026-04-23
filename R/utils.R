#' String concatenation operator
#'
#' Shorthand for `paste0(a, b)`.
#'
#' @param a A character string.
#' @param b A character string.
#'
#' @return A character string: `paste0(a, b)`.
#'
#' @export
`%&%` <- function(a, b) paste0(a, b)

# ----------------------------------------

#' Null-coalescing operator
#'
#' Returns `lhs` if not `NULL`, otherwise returns `rhs`.
#'
#' @param lhs An object.
#' @param rhs Fallback value returned when `lhs` is `NULL`.
#'
#' @return `lhs` if not `NULL`, else `rhs`.
#'
#' @name op-null-default
#' @aliases %||%
#' @keywords internal
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

# ----------------------------------------

#' Calculate string similarity between two taxon names
#'
#' Computes a normalised similarity score (0--1) between two strings using
#' Levenshtein edit distance. Used internally to validate fuzzy WoRMS matches.
#'
#' @param str1 A character string.
#' @param str2 A character string.
#'
#' @return A numeric value between 0 (no similarity) and 1 (identical).
#'   Returns 0 if either input is `NA` or an empty string.
#'
#' @importFrom stringdist stringdist
#'
#' @keywords internal
calculate_similarity <- function(str1, str2) {
  if (is.na(str1) || is.na(str2) || str1 == "" || str2 == "") {
    return(0)
  }
  lev_dist <- stringdist::stringdist(
    tolower(str1),
    tolower(str2),
    method = "lv"
  )
  max_len <- max(nchar(str1), nchar(str2))
  if (max_len == 0) {
    return(1)
  }
  1 - lev_dist / max_len
}

# ----------------------------------------

#' Query WoRMS for records by name (internal, isolated for testability)
#'
#' Thin wrapper around [worrms::wm_records_name()] with a timeout guard.
#' Exists as a standalone function so that tests can mock it independently
#' of `worms_with_timeout()` without triggering eager evaluation issues.
#'
#' @param genus_name A character string with the genus name to query.
#' @param timeout_sec Timeout in seconds. Default `10`.
#'
#' @return A data frame of WoRMS records, or `NULL` on failure/timeout.
#'
#' @importFrom worrms wm_records_name
#'
#' @keywords internal
.wm_records_name_safe <- function(genus_name, timeout_sec = 10) {
  worms_with_timeout(
    worrms::wm_records_name(genus_name, fuzzy = FALSE, marine_only = FALSE),
    timeout_sec = timeout_sec
  )
}

# ----------------------------------------

#' Check whether a genus name exists in WoRMS
#'
#' Queries WoRMS to verify that a given genus name is a valid accepted
#' genus. Used internally to guard against accepting fuzzy matches whose
#' genus component does not exist in WoRMS.
#'
#' Delegates the actual API call to `.wm_records_name_safe()`, which is
#' isolated to allow mocking in unit tests without eager-evaluation issues.
#'
#' @param genus_name A character string with the genus name to validate.
#' @param timeout_sec Timeout in seconds passed to the API wrapper.
#'   Default `10`.
#'
#' @return `TRUE` if WoRMS returns at least one record ranked as genus,
#'   `FALSE` otherwise (including when `genus_name` is `NA` or empty).
#'
#' @keywords internal
verify_genus_exists <- function(genus_name, timeout_sec = 10) {
  if (is.na(genus_name) || genus_name == "") {
    return(FALSE)
  }
  result <- .wm_records_name_safe(genus_name, timeout_sec = timeout_sec)
  if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
    genus_matches <- result[tolower(result$rank) == "genus", ]
    return(nrow(genus_matches) > 0)
  }
  FALSE
}

# ----------------------------------------

#' Convert a character vector to a minimal pipeline-compatible data frame
#'
#' Takes a character vector of taxon names and wraps it in a one-column
#' data frame with the internal column name `taxon_clean`, which is the
#' working column used throughout the resolution pipeline. This is used
#' internally to allow pipeline functions to accept bare character vectors
#' in addition to data frames.
#'
#' @param x A character vector of taxon names. `NA` values are preserved.
#'
#' @return A data frame with a single column `taxon_clean`.
#'
#' @keywords internal
names_to_df <- function(x) {
  if (!is.character(x)) {
    stop("`x` must be a character vector.")
  }
  data.frame(taxon_clean = x, stringsAsFactors = FALSE)
}

# ----------------------------------------

#' Resolve and normalise the taxon column inside a data frame
#'
#' Validates that the user-supplied column name exists in the data frame
#' and, if it differs from `taxon_clean`, renames it to `taxon_clean` so
#' the rest of the pipeline can work with a single known column name.
#' The original column name is stored in the attribute `"original_col"`
#' of the returned data frame so it can be restored later if needed.
#'
#' @param df A data frame.
#' @param col A character string with the name of the column that contains
#'   taxon names. Default: `"taxon_clean"`.
#'
#' @return `df` with the taxon column renamed to `taxon_clean` (if
#'   necessary), and attribute `"original_col"` set to `col`.
#'
#' @keywords internal
resolve_col <- function(df, col = "taxon_clean") {
  if (!col %in% names(df)) {
    stop(sprintf(
      "Column '%s' not found in data frame. Available columns: %s",
      col,
      paste(names(df), collapse = ", ")
    ))
  }
  if (col != "taxon_clean") {
    names(df)[names(df) == col] <- "taxon_clean"
  }
  attr(df, "original_col") <- col
  df
}

# ----------------------------------------

#' Format pipeline output for direct query mode
#'
#' Selects and returns only the most relevant columns from a resolved
#' data frame when pipeline functions are called with a bare character
#' vector instead of a full data frame. Keeps the output readable in
#' an interactive session without exposing the 30+ taxonomy rank columns
#' that are mostly `NA` at query time.
#'
#' @param df A data frame returned by any pipeline resolution function.
#' @param original_col Character. The original taxon column name supplied
#'   by the user (stored in `attr(df, "original_col")`). Used to rename
#'   the output column back to what the user expects. Default: `"taxon_clean"`.
#'
#' @return A [tibble::tibble()] with columns: the taxon column (renamed to
#'   `original_col`), `matched_name`, `accepted_name`, `taxonomic_status`,
#'   `resolution_method`, and `resolution_notes`.
#'
#' @importFrom tibble as_tibble
#'
#' @keywords internal
format_query_result <- function(df, original_col = "taxon_clean") {
  keep <- c(
    "taxon_clean",
    "matched_name",
    "accepted_name",
    "taxonomic_status",
    "resolution_method",
    "resolution_notes"
  )
  # Keep only columns that actually exist (pipeline may not have all yet)
  keep <- keep[keep %in% names(df)]
  out <- df[, keep, drop = FALSE]

  # Restore original column name
  if (original_col != "taxon_clean" && "taxon_clean" %in% names(out)) {
    names(out)[names(out) == "taxon_clean"] <- original_col
  }
  tibble::as_tibble(out)
}

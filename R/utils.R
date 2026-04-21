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

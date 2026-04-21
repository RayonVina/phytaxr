#' Execute a WoRMS API call with a timeout
#'
#' Wraps a `worrms` API call in a timeout guard. Returns `NULL` silently
#' if the call times out or errors, instead of crashing the pipeline.
#'
#' @param expr An expression calling a `worrms` function.
#' @param timeout_sec Timeout in seconds. Default: `10`.
#'
#' @return The result of `expr`, or `NULL` on timeout or error.
#'
#' @keywords internal
worms_with_timeout <- function(expr, timeout_sec = 10) {
  tryCatch(
    R.utils::withTimeout(expr, timeout = timeout_sec, onTimeout = "silent"),
    error = function(e) NULL,
    TimeoutException = function(e) NULL
  )
}

# ----------------------------------------

#' Robust WoRMS name search with retries
#'
#' Searches WoRMS for a taxon name using either exact match or fuzzy
#' matching. Retries on failure and respects a rate-limit delay between
#' attempts.
#'
#' @param name A character string with the taxon name to search.
#' @param method Either `"exact"` (default) or `"fuzzy"`.
#' @param max_retries Maximum number of attempts. Default: `2`.
#' @param timeout_sec Timeout per attempt in seconds. Default: `10`.
#' @param rate_delay Numeric vector of length 2: min and max seconds to
#'   wait between attempts. Default: `c(0.5, 1)`.
#'
#' @return A named list with elements `AphiaID`, `scientificName`,
#'   `status`, `valid_name`, and `valid_AphiaID`, or `NULL` if not found.
#'
#' @importFrom worrms wm_name2id wm_records_name wm_record
#'
#' @keywords internal
worms_robust_search <- function(
  name,
  method = "exact",
  max_retries = 2,
  timeout_sec = 10,
  rate_delay = c(0.5, 1)
) {
  if (is.na(name) || trimws(name) == "") {
    return(NULL)
  }

  for (attempt in seq_len(max_retries)) {
    Sys.sleep(stats::runif(1, rate_delay[1], rate_delay[2]))

    aphia_id <- if (method == "exact") {
      worms_with_timeout(worrms::wm_name2id(name), timeout_sec)
    } else {
      records <- worms_with_timeout(
        worrms::wm_records_name(name, fuzzy = TRUE, marine_only = FALSE),
        timeout_sec
      )
      if (!is.null(records) && nrow(records) > 0) records$AphiaID[1] else NA
    }

    if (!is.null(aphia_id) && !is.na(aphia_id)) {
      full_record <- worms_with_timeout(
        worrms::wm_record(aphia_id),
        timeout_sec
      )
      if (!is.null(full_record)) {
        return(list(
          AphiaID = aphia_id,
          scientificName = full_record$scientificname,
          status = full_record$status,
          valid_name = full_record$valid_name,
          valid_AphiaID = full_record$valid_AphiaID
        ))
      }
    }
  }
  NULL
}

# ----------------------------------------

#' Look up a taxon name in a pre-loaded dictionary
#'
#' Searches a named list (dictionary) for a taxon name and returns
#' resolution metadata. Returns a tibble with `found = FALSE` if the
#' name is not present.
#'
#' @param search_term Character. The taxon name to look up.
#' @param lookup A named list where each name is a taxon string and each
#'   value is a one-row data frame with resolution columns.
#' @param include_notes Logical. Whether to include the `marginalia`
#'   column. Default: `TRUE`.
#'
#' @return A one-row [tibble::tibble()] with columns `matched_aphiaid`,
#'   `matched_name`, `accepted_name`, `accepted_aphiaid`,
#'   `taxonomic_status`, `resolution_method`, `aphiaid`, `marginalia`,
#'   and `found`.
#'
#' @importFrom tibble tibble
#'
#' @export
lookup_taxon_info <- function(search_term, lookup, include_notes = TRUE) {
  result <- lookup[[search_term]]

  if (is.null(result) || nrow(result) == 0) {
    return(tibble::tibble(
      matched_aphiaid = NA_integer_,
      matched_name = NA_character_,
      accepted_name = NA_character_,
      accepted_aphiaid = NA_integer_,
      taxonomic_status = NA_character_,
      resolution_method = NA_character_,
      aphiaid = NA_integer_,
      marginalia = NA_character_,
      found = FALSE
    ))
  }

  tibble::tibble(
    matched_aphiaid = as.integer(result$matched_aphiaid[1]),
    matched_name = as.character(result$matched_name[1]),
    accepted_name = as.character(result$accepted_name[1]),
    accepted_aphiaid = as.integer(result$accepted_aphiaid[1]),
    taxonomic_status = as.character(result$taxonomic_status[1]),
    resolution_method = as.character(result$resolution_method[1]),
    aphiaid = as.integer(result$aphiaid[1]),
    marginalia = if (include_notes && !is.null(result$marginalia[1])) {
      as.character(result$marginalia[1])
    } else {
      NA_character_
    },
    found = TRUE
  )
}

# ----------------------------------------

#' Look up full taxonomy from an AphiaID dictionary
#'
#' Given an AphiaID, retrieves the full taxonomic classification from a
#' pre-loaded named list keyed by AphiaID.
#'
#' @param aphiaid_val Integer or numeric. The AphiaID to look up.
#' @param lookup A named list where each name is an AphiaID (as character)
#'   and each value is a one-row data frame with taxonomic rank columns.
#'
#' @return A one-row [tibble::tibble()] with all taxonomic rank columns
#'   from kingdom to forma, plus `rank` and `taxonomic_status`. All
#'   columns are `NA_character_` if the AphiaID is not found.
#'
#' @importFrom tibble tibble
#'
#' @export
lookup_taxonomy_info <- function(aphiaid_val, lookup) {
  .empty <- function() {
    tibble::tibble(
      kingdom = NA_character_,
      subkingdom = NA_character_,
      infrakingdom = NA_character_,
      phylum = NA_character_,
      subphylum = NA_character_,
      infraphylum = NA_character_,
      parvphylum = NA_character_,
      gigaclass = NA_character_,
      superclass = NA_character_,
      class = NA_character_,
      subclass = NA_character_,
      infraclass = NA_character_,
      subterclass = NA_character_,
      superorder = NA_character_,
      order = NA_character_,
      suborder = NA_character_,
      infraorder = NA_character_,
      parvorder = NA_character_,
      superfamily = NA_character_,
      family = NA_character_,
      subfamily = NA_character_,
      tribe = NA_character_,
      genus = NA_character_,
      subgenus = NA_character_,
      section = NA_character_,
      subsection = NA_character_,
      species = NA_character_,
      subspecies = NA_character_,
      variety = NA_character_,
      forma = NA_character_,
      rank = NA_character_,
      taxonomic_status = NA_character_
    )
  }

  if (is.na(aphiaid_val)) {
    return(.empty())
  }

  result <- lookup[[as.character(aphiaid_val)]]
  if (is.null(result) || nrow(result) == 0) {
    return(.empty())
  }

  tibble::tibble(
    kingdom = as.character(result$kingdom[1]),
    subkingdom = as.character(result$subkingdom[1]),
    infrakingdom = as.character(result$infrakingdom[1]),
    phylum = as.character(result$phylum[1]),
    subphylum = as.character(result$subphylum[1]),
    infraphylum = as.character(result$infraphylum[1]),
    parvphylum = as.character(result$parvphylum[1]),
    gigaclass = as.character(result$gigaclass[1]),
    superclass = as.character(result$superclass[1]),
    class = as.character(result$class[1]),
    subclass = as.character(result$subclass[1]),
    infraclass = as.character(result$infraclass[1]),
    subterclass = as.character(result$subterclass[1]),
    superorder = as.character(result$superorder[1]),
    order = as.character(result$order[1]),
    suborder = as.character(result$suborder[1]),
    infraorder = as.character(result$infraorder[1]),
    parvorder = as.character(result$parvorder[1]),
    superfamily = as.character(result$superfamily[1]),
    family = as.character(result$family[1]),
    subfamily = as.character(result$subfamily[1]),
    tribe = as.character(result$tribe[1]),
    genus = as.character(result$genus[1]),
    subgenus = as.character(result$subgenus[1]),
    section = as.character(result$section[1]),
    subsection = as.character(result$subsection[1]),
    species = as.character(result$species[1]),
    subspecies = as.character(result$subspecies[1]),
    variety = as.character(result$variety[1]),
    forma = as.character(result$forma[1]),
    rank = as.character(result$rank[1]),
    taxonomic_status = as.character(result$taxonomic_status[1])
  )
}

# ----------------------------------------

#' Ensure the resolution schema is present in a data frame
#'
#' Checks for all columns required by the automatic resolution pipeline
#' and adds any that are missing, initialised to \code{NA} with the correct
#' type. Existing columns and their values are never modified.
#'
#' This function is called internally at the start of each resolution
#' function (\code{search_worms_priority()}, \code{search_worms_taxamatch()},
#' \code{search_gbif_strict()}, \code{resolve_taxonomic_status()},
#' \code{search_worms_fuzzy_minor()}, \code{get_taxonomy()}, and
#' \code{process_fuzzy_batch()}), so the user does not need to call it
#' manually. When all columns are already present the function returns
#' \code{df} unchanged with virtually no overhead.
#'
#' @param df A data frame, typically the output of the cleaning pipeline.
#'
#' @return The input \code{df} with any missing resolution or taxonomic rank
#'   columns added and initialised to \code{NA}.
#'
#' @keywords internal
ensure_resolution_schema <- function(df) {
  resolution_cols <- list(
    matched_aphiaid   = NA_integer_,
    matched_name      = NA_character_,
    accepted_name     = NA_character_,
    accepted_aphiaid  = NA_integer_,
    taxonomic_status  = NA_character_,
    resolution_method = NA_character_,
    resolution_notes  = NA_character_,
    rank              = NA_character_,
    flag_for_removal  = FALSE
  )

  tax_ranks <- c(
    "kingdom", "subkingdom", "infrakingdom", "phylum", "subphylum",
    "infraphylum", "parvphylum", "gigaclass", "superclass", "class",
    "subclass", "infraclass", "subterclass", "superorder", "order",
    "suborder", "infraorder", "parvorder", "superfamily", "family",
    "subfamily", "tribe", "genus", "subgenus", "section", "subsection",
    "species", "subspecies", "variety", "forma"
  )
  tax_cols <- setNames(
    replicate(length(tax_ranks), NA_character_, simplify = FALSE),
    tax_ranks
  )

  all_cols <- c(resolution_cols, tax_cols)
  missing  <- setdiff(names(all_cols), names(df))
  if (length(missing) > 0) {
    df[missing] <- all_cols[missing]
  }
  df
}

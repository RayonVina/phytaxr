# R/strict.R
# Automatic strict resolution pipeline -- Step 5 of the taxonomy pipeline.
# Mirrors §5 of master_taxonomy.R exactly.
# Functions are called sequentially on the unique_taxon_clean data frame.

# ============================================================================
# UTILITY FUNCTIONS (§5.2)
# ============================================================================

#' Call the GBIF species match API
#'
#' Queries the GBIF `/species/match` endpoint for a taxon name. Returns the
#' parsed JSON response if a `usageKey` is found, or `NULL` otherwise.
#'
#' @param name Character. Scientific name to query.
#' @param max_retries Integer. Number of attempts before giving up. Default `2`.
#' @param timeout_sec Timeout per attempt in seconds. Default `10`.
#'
#' @return A parsed list from the GBIF API, or `NULL`.
#'
#' @importFrom httr GET timeout user_agent status_code content
#'
#' @keywords internal
call_gbif_api <- function(name, max_retries = 2, timeout_sec = 10) {
  base_url <- "https://api.gbif.org/v1/species/match"
  for (attempt in seq_len(max_retries)) {
    Sys.sleep(stats::runif(1, 0.3, 0.7))
    tryCatch(
      {
        response <- httr::GET(
          url = base_url,
          query = list(name = name),
          httr::timeout(timeout_sec),
          httr::user_agent("R-PhytoMAP-TaxonomyResolver/1.0")
        )
        if (httr::status_code(response) == 200) {
          content <- httr::content(response, as = "parsed")
          if (!is.null(content$usageKey)) return(content)
        }
      },
      error = function(e) NULL
    )
    if (attempt < max_retries) Sys.sleep(0.5)
  }
  NULL
}

# ----------------------------------------

#' Cross-validate a name against WoRMS
#'
#' Checks whether a scientific name (e.g. returned by GBIF) exists in WoRMS.
#' Used as a validation gate in the GBIF strict step.
#'
#' @param gbif_matched_name Character. Name to check in WoRMS.
#'
#' @return Logical. `TRUE` if the name is found in WoRMS, `FALSE` otherwise.
#'
#' @importFrom worrms wm_name2id
#'
#' @keywords internal
cross_validate_with_worms <- function(gbif_matched_name) {
  if (is.na(gbif_matched_name) || gbif_matched_name == "") {
    return(FALSE)
  }
  worms_check <- worms_with_timeout(worrms::wm_name2id(gbif_matched_name))
  !is.null(worms_check) && !is.na(worms_check)
}

# ============================================================================
# PIPELINE FUNCTIONS (§5.3)
# ============================================================================

#' Step 5.1 -- WoRMS exact search
#'
#' For each unresolved row in `df`, attempts an exact name lookup in WoRMS
#' using [worms_robust_search()]. Fills `matched_aphiaid`, `matched_name`,
#' `accepted_name`, `taxonomic_status`, and `resolution_method`.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame with at least column `taxon_clean`.
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#'
#' @return The input `df` with resolved rows updated in-place.
#'
#' @importFrom progress progress_bar
#'
#' @export
search_worms_priority <- function(df, col = "taxon_clean") {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.1 WoRMS Priority (Exact) ---\n")

  unresolved_indices <- which(is.na(df$matched_aphiaid))
  if (length(unresolved_indices) == 0) {
    cat("  All entries already resolved\n\n")
    return(df)
  }

  cat(sprintf("  Processing: %d entries\n", length(unresolved_indices)))

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = length(unresolved_indices),
    clear = FALSE,
    width = 80
  )

  resolved_count <- 0
  for (i in unresolved_indices) {
    result <- worms_robust_search(df$taxon_clean[i], method = "exact")
    if (!is.null(result)) {
      df$matched_aphiaid[i] <- as.integer(result$AphiaID)
      df$matched_name[i] <- result$scientificName
      df$accepted_name[i] <- ifelse(
        !is.null(result$valid_name) && !is.na(result$valid_name),
        result$valid_name,
        result$scientificName
      )
      df$taxonomic_status[i] <- result$status
      df$resolution_method[i] <- "worms_exact"
      if (
        !is.null(result$valid_AphiaID) &&
          !is.na(result$valid_AphiaID) &&
          result$AphiaID != result$valid_AphiaID
      ) {
        df$accepted_aphiaid[i] <- as.integer(result$valid_AphiaID)
        df$resolution_notes[i] <- paste0(
          "Synonym: ",
          result$scientificName,
          " \u2192 ",
          result$valid_name
        )
      }
      resolved_count <- resolved_count + 1
    }
    pb$tick()
  }

  cat(sprintf("  Resolved: %d\n\n", resolved_count))
  df
}

# ----------------------------------------

#' Step 5.2 -- WoRMS Taxamatch
#'
#' For each unresolved row with a name of >= 6 characters (excluding
#' "unidentified / unknown / not classified"), queries the WoRMS Taxamatch
#' service. Applies the best match when a valid record is returned.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame as returned by [search_worms_priority()].
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#'
#' @return Updated `df`.
#'
#' @importFrom progress progress_bar
#' @importFrom stringr str_detect
#' @importFrom worrms wm_records_taxamatch
#'
#' @export
search_worms_taxamatch <- function(df, col = "taxon_clean") {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.2 WoRMS Taxamatch ---\n")

  unresolved_indices <- which(is.na(df$matched_aphiaid))
  if (length(unresolved_indices) == 0) {
    cat("  No unresolved entries\n\n")
    return(df)
  }

  unresolved_names <- df$taxon_clean[unresolved_indices]
  candidate_names <- unresolved_names[
    nchar(unresolved_names) >= 6 &
      !stringr::str_detect(
        tolower(unresolved_names),
        "\\bunidentified|\\bunknown|\\bnot\\s+classified"
      )
  ]

  if (length(candidate_names) == 0) {
    cat("  No suitable candidates\n\n")
    return(df)
  }

  cat(sprintf("  Processing: %d candidates\n", length(candidate_names)))

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = length(candidate_names),
    clear = FALSE,
    width = 80
  )

  resolved_count <- 0
  for (name in candidate_names) {
    Sys.sleep(stats::runif(1, 0.5, 1))
    result_list <- worms_with_timeout(
      worrms::wm_records_taxamatch(name = name, marine_only = FALSE)
    )

    if (
      !is.null(result_list) &&
        length(result_list) > 0 &&
        nrow(result_list[[1]]) > 0
    ) {
      best_match <- result_list[[1]][1, ]
      if (
        !is.na(best_match$scientificname) &&
          nchar(best_match$scientificname) >= 3 &&
          !is.na(best_match$AphiaID)
      ) {
        matching_indices <- which(
          df$taxon_clean == name & is.na(df$matched_aphiaid)
        )
        for (idx in matching_indices) {
          df$matched_aphiaid[idx] <- as.integer(best_match$AphiaID)
          df$matched_name[idx] <- best_match$scientificname
          df$accepted_name[idx] <- best_match$valid_name %||%
            best_match$scientificname
          df$taxonomic_status[idx] <- best_match$status
          df$resolution_method[idx] <- "worms_taxamatch"
          if (
            !is.null(best_match$valid_AphiaID) &&
              !is.na(best_match$valid_AphiaID) &&
              best_match$AphiaID != best_match$valid_AphiaID
          ) {
            df$accepted_aphiaid[idx] <- as.integer(best_match$valid_AphiaID)
            df$resolution_notes[idx] <- paste0(
              "Taxamatch synonym: ",
              best_match$scientificname,
              " \u2192 ",
              best_match$valid_name
            )
          } else if (name != best_match$scientificname) {
            df$resolution_notes[idx] <- paste0(
              "Taxamatch correction: ",
              name,
              " \u2192 ",
              best_match$scientificname
            )
          }
        }
        resolved_count <- resolved_count + length(matching_indices)
      }
    }
    pb$tick()
  }

  cat(sprintf("  Resolved: %d\n\n", resolved_count))
  df
}

# ----------------------------------------

#' Step 5.3 -- GBIF strict (confidence >= 99, EXACT, cross-validated with WoRMS)
#'
#' For each remaining unresolved row, queries GBIF and accepts the match only
#' when confidence is >= 99, matchType is EXACT, and the returned name is
#' found in WoRMS.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame as returned by [search_worms_taxamatch()].
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#'
#' @return Updated `df`.
#'
#' @importFrom progress progress_bar
#'
#' @export
search_gbif_strict <- function(df, col = "taxon_clean") {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.3 GBIF Strict ---\n")

  unresolved_indices <- which(is.na(df$matched_aphiaid))
  if (length(unresolved_indices) == 0) {
    cat("  No unresolved entries\n\n")
    return(df)
  }

  cat(sprintf("  Processing: %d entries\n", length(unresolved_indices)))

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = length(unresolved_indices),
    clear = FALSE,
    width = 80
  )

  resolved_count <- 0
  for (i in unresolved_indices) {
    gbif_result <- call_gbif_api(df$taxon_clean[i])
    if (
      !is.null(gbif_result) &&
        !is.null(gbif_result$confidence) &&
        as.integer(gbif_result$confidence) >= 99 &&
        !is.null(gbif_result$matchType) &&
        gbif_result$matchType == "EXACT"
    ) {
      matched_name <- gbif_result$scientificName
      if (cross_validate_with_worms(matched_name)) {
        df$matched_aphiaid[i] <- NA_integer_
        df$matched_name[i] <- matched_name
        df$accepted_name[i] <- gbif_result$canonicalName %||% matched_name
        df$taxonomic_status[i] <- gbif_result$taxonomicStatus %||% "unknown"
        df$resolution_method[i] <- "gbif_strict"
        df$resolution_notes[i] <- paste0(
          "GBIF usageKey: ",
          gbif_result$usageKey
        )
        resolved_count <- resolved_count + 1
      }
    }
    pb$tick()
  }

  cat(sprintf("  Resolved: %d\n\n", resolved_count))
  df
}

# ----------------------------------------

#' Step 5.4 -- Taxonomic status resolution
#'
#' For all rows resolved via WoRMS, fetches the full WoRMS record and, if the
#' status is not "accepted", follows the `valid_AphiaID` pointer to fill
#' `accepted_aphiaid` and `accepted_name`.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame as returned by [search_gbif_strict()].
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#'
#' @return Updated `df`.
#'
#' @importFrom progress progress_bar
#' @importFrom worrms wm_record
#'
#' @export
resolve_taxonomic_status <- function(df, col = "taxon_clean") {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.4 Taxonomic Status Resolution ---\n")

  resolved_indices <- which(!is.na(df$matched_aphiaid))
  if (length(resolved_indices) == 0) {
    cat("  No resolved records to process\n\n")
    return(df)
  }

  needs_check <- resolved_indices[
    is.na(df$accepted_aphiaid[resolved_indices]) |
      is.na(df$accepted_name[resolved_indices])
  ]

  if (length(needs_check) == 0) {
    cat("  All resolved records already have accepted name\n\n")
    return(df)
  }

  cat(sprintf("  Checking: %d records\n", length(needs_check)))

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = length(needs_check),
    clear = FALSE,
    width = 80
  )

  for (i in needs_check) {
    record <- worms_with_timeout(worrms::wm_record(df$matched_aphiaid[i]))
    if (!is.null(record)) {
      df$taxonomic_status[i] <- record$status
      if (
        !is.null(record$valid_AphiaID) &&
          !is.na(record$valid_AphiaID) &&
          record$AphiaID != record$valid_AphiaID
      ) {
        df$accepted_aphiaid[i] <- as.integer(record$valid_AphiaID)
        valid_record <- worms_with_timeout(
          worrms::wm_record(record$valid_AphiaID)
        )
        if (!is.null(valid_record)) {
          df$accepted_name[i] <- valid_record$scientificname
          if (is.na(df$resolution_notes[i]) || df$resolution_notes[i] == "") {
            df$resolution_notes[i] <- paste0(
              "Synonym resolved: ",
              record$scientificname,
              " \u2192 ",
              valid_record$scientificname
            )
          }
        }
      } else {
        if (is.na(df$accepted_name[i])) {
          df$accepted_name[i] <- record$scientificname
        }
      }
    }
    pb$tick()
  }

  cat("  Done\n\n")
  df
}

# ----------------------------------------

#' Step 5.5 -- WoRMS fuzzy search (minor corrections)
#'
#' For each remaining unresolved row whose name has >= 6 characters and
#' contains at least two words (excluding names with "classified"), queries
#' WoRMS Taxamatch and accepts the result only when the Levenshtein
#' similarity exceeds `min_similarity_score` and the genus is confirmed to
#' exist in WoRMS.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame as returned by [resolve_taxonomic_status()].
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#' @param fuzzy_config Named list with elements `max_levenshtein_distance`
#'   (integer), `min_similarity_score` (numeric 0-1), and
#'   `require_genus_match` (logical, default `TRUE`).
#' @param ncores Integer. Number of parallel workers. Default `2`.
#'
#' @return Updated `df`.
#'
#' @importFrom dplyr mutate filter
#' @importFrom furrr future_map2 furrr_options
#' @importFrom future plan multisession sequential
#' @importFrom progress progress_bar
#' @importFrom purrr map
#' @importFrom stringr str_detect str_count word
#' @importFrom worrms wm_records_taxamatch
#'
#' @export
search_worms_fuzzy_minor <- function(
  df,
  col = "taxon_clean",
  fuzzy_config = list(
    max_levenshtein_distance = 3,
    min_similarity_score = 0.85,
    require_genus_match = TRUE
  ),
  ncores = 2
) {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.5 WoRMS Fuzzy (Minor Corrections) ---\n")

  unresolved_indices <- which(is.na(df$matched_aphiaid))
  if (length(unresolved_indices) == 0) {
    cat("  No unresolved entries\n\n")
    return(df)
  }

  candidates_df <- df[unresolved_indices, ] |>
    dplyr::mutate(search_idx = unresolved_indices) |>
    dplyr::filter(
      nchar(taxon_clean) >= 6,
      !stringr::str_detect(tolower(taxon_clean), "classified"),
      stringr::str_count(taxon_clean, "\\S+") >= 2
    )

  if (nrow(candidates_df) == 0) {
    cat("  No suitable candidates\n\n")
    return(df)
  }

  cat(sprintf("  Candidates: %d\n", nrow(candidates_df)))
  cat("  Fetching WoRMS data...\n")

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = nrow(candidates_df),
    clear = FALSE,
    width = 80
  )

  candidates_df <- candidates_df |>
    dplyr::mutate(
      worms_result = purrr::map(taxon_clean, function(name) {
        Sys.sleep(stats::runif(1, 0.5, 1))
        pb$tick()
        worms_with_timeout(
          worrms::wm_records_taxamatch(name = name, marine_only = FALSE)
        )
      })
    )

  cat("  Calculating similarities...\n")

  future::plan(future::multisession, workers = ncores)

  candidates_df <- candidates_df |>
    dplyr::mutate(
      match_info = furrr::future_map2(
        taxon_clean,
        worms_result,
        function(original_name, result_list) {
          if (!is.null(result_list) && length(result_list) > 0) {
            for (result_df in result_list) {
              if (nrow(result_df) > 0) {
                best_match <- result_df[1, ]
                if (
                  !is.na(best_match$scientificname) &&
                    !is.na(best_match$AphiaID)
                ) {
                  similarity <- calculate_similarity(
                    original_name,
                    best_match$scientificname
                  )
                  if (similarity >= fuzzy_config$min_similarity_score) {
                    genus_ok <- if (fuzzy_config$require_genus_match) {
                      verify_genus_exists(stringr::word(original_name, 1))
                    } else {
                      TRUE
                    }
                    if (genus_ok) {
                      return(list(
                        found = TRUE,
                        aphiaid = as.integer(best_match$AphiaID),
                        name = best_match$scientificname,
                        valid_name = best_match$valid_name,
                        valid_id = as.integer(best_match$valid_AphiaID),
                        status = best_match$status,
                        similarity = similarity
                      ))
                    }
                  }
                }
              }
            }
          }
          list(found = FALSE)
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    )

  future::plan(future::sequential)

  resolved_count <- 0
  for (row_i in seq_len(nrow(candidates_df))) {
    info <- candidates_df$match_info[[row_i]]
    if (isTRUE(info$found)) {
      idx <- candidates_df$search_idx[row_i]
      df$matched_aphiaid[idx] <- info$aphiaid
      df$matched_name[idx] <- info$name
      df$accepted_name[idx] <- info$valid_name %||% info$name
      df$taxonomic_status[idx] <- info$status
      df$resolution_method[idx] <- "worms_fuzzy_minor"
      df$resolution_notes[idx] <- sprintf(
        "Fuzzy correction (similarity=%.2f): %s \u2192 %s",
        info$similarity,
        candidates_df$taxon_clean[row_i],
        info$name
      )
      if (
        !is.null(info$valid_id) &&
          !is.na(info$valid_id) &&
          info$aphiaid != info$valid_id
      ) {
        df$accepted_aphiaid[idx] <- info$valid_id
      }
      resolved_count <- resolved_count + 1
    }
  }

  cat(sprintf("  Resolved: %d\n\n", resolved_count))
  df
}

# ============================================================================
# TAXONOMY RETRIEVAL (§5.6)
# ============================================================================

#' Step 5.6 -- Retrieve full WoRMS taxonomy
#'
#' For every row with a non-`NA` `matched_aphiaid` (or `accepted_aphiaid`),
#' fetches the full WoRMS classification and populates the taxonomic rank
#' columns (`kingdom` through `forma`) plus `rank`.
#'
#' Calls [ensure_resolution_schema()] on entry so that the function can be
#' used directly on the output of the cleaning pipeline without requiring the
#' user to pre-initialise the resolution columns.
#'
#' @param df A data frame with at least `matched_aphiaid` and
#'   `accepted_aphiaid` columns.
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#'
#' @return Updated `df` with taxonomic rank columns populated.
#'
#' @importFrom progress progress_bar
#' @importFrom worrms wm_classification wm_record
#'
#' @export
get_taxonomy <- function(df, col = "taxon_clean") {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  cat("--- 5.6 Taxonomy Retrieval ---\n")

  tax_ranks <- c(
    "kingdom",
    "subkingdom",
    "infrakingdom",
    "phylum",
    "subphylum",
    "infraphylum",
    "parvphylum",
    "gigaclass",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "subterclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "subgenus",
    "section",
    "subsection",
    "species",
    "subspecies",
    "variety",
    "forma"
  )

  resolved_indices <- which(!is.na(df$matched_aphiaid))
  if (length(resolved_indices) == 0) {
    cat("  No resolved records\n\n")
    return(df)
  }

  cat(sprintf("  Retrieving: %d records\n", length(resolved_indices)))

  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent (:current/:total) :eta",
    total = length(resolved_indices),
    clear = FALSE,
    width = 80
  )

  for (i in resolved_indices) {
    target_id <- if (!is.na(df$accepted_aphiaid[i])) {
      df$accepted_aphiaid[i]
    } else {
      df$matched_aphiaid[i]
    }

    if (!is.na(target_id)) {
      taxonomy <- worms_with_timeout(worrms::wm_classification(target_id))
      if (!is.null(taxonomy) && nrow(taxonomy) > 0) {
        for (j in seq_len(nrow(taxonomy))) {
          rank_name <- tolower(taxonomy$rank[j])
          if (rank_name %in% tax_ranks) {
            df[[rank_name]][i] <- taxonomy$scientificname[j]
          }
        }
      }
      record <- worms_with_timeout(worrms::wm_record(target_id))
      if (!is.null(record)) {
        df$rank[i] <- record$rank
        if (is.na(df$taxonomic_status[i]) || df$taxonomic_status[i] == "") {
          df$taxonomic_status[i] <- record$status
        }
      }
    }
    pb$tick()
  }

  cat("  Done\n\n")
  df
}

# ============================================================================
# PIPELINE ORCHESTRATOR (§5.4)
# ============================================================================

#' Run the full Step 5 automatic resolution pipeline
#'
#' Executes the six pipeline steps sequentially on `df` and prints a
#' resolution summary at the end. Equivalent to the §5.4 block in
#' `master_taxonomy.R`.
#'
#' Calls [ensure_resolution_schema()] on entry so that the pipeline can be
#' applied directly to the output of the cleaning pipeline.
#'
#' @param df A data frame produced by the cleaning pipeline (Steps 1-4),
#'   with at least a `taxon_clean` column.
#' @param col Character. Name of the column containing the cleaned taxon
#'   names. Default `"taxon_clean"`.
#' @param fuzzy_config Named list passed to [search_worms_fuzzy_minor()].
#' @param ncores Integer. Parallel workers for the fuzzy step. Default `2`.
#'
#' @return The fully resolved `df`.
#'
#' @export
run_resolution_pipeline <- function(
  df,
  col = "taxon_clean",
  fuzzy_config = list(
    max_levenshtein_distance = 3,
    min_similarity_score = 0.85,
    require_genus_match = TRUE
  ),
  ncores = 2
) {
  query_mode <- is.character(df)
  if (query_mode) {
    df <- names_to_df(df)
  }
  df <- resolve_col(df, col)
  df <- ensure_resolution_schema(df)
  df <- search_worms_priority(df)
  df <- search_worms_taxamatch(df)
  df <- search_gbif_strict(df)
  df <- resolve_taxonomic_status(df)
  df <- search_worms_fuzzy_minor(
    df,
    fuzzy_config = fuzzy_config,
    ncores = ncores
  )
  df <- get_taxonomy(df)
  if (query_mode) {
    return(format_query_result(df, attr(df, "original_col") %||% col))
  }

  cat("=== Resolution Summary ===\n")
  cat(sprintf("  Total:      %d\n", nrow(df)))
  cat(sprintf("  Resolved:   %d\n", sum(!is.na(df$matched_aphiaid))))
  cat(sprintf("  Unresolved: %d\n", sum(is.na(df$matched_aphiaid))))
  if ("resolution_method" %in% names(df)) {
    tbl <- table(df$resolution_method, useNA = "ifany")
    for (nm in names(tbl)) {
      label <- if (is.na(nm)) "unresolved" else nm
      cat(sprintf("    %-30s %d\n", label, tbl[[nm]]))
    }
  }
  df
}

# R/fuzzy.R
# Semi-automatic fuzzy resolution -- Step 6 of the taxonomy pipeline.
# Provides vocabulary builders, WoRMS REST wrappers, edit-distance helpers,
# the cascading suggestion orchestrator, and the interactive batch loop.

# ============================================================================
# VOCABULARY BUILDERS
# ============================================================================

#' Build a genus vocabulary from resolved taxa
#'
#' Extracts all unique genus names (>= 4 characters) from the `genus` column
#' of a resolved data frame. Used to power edit-distance genus correction in
#' the semi-automatic fuzzy pipeline.
#'
#' @param resolved_df A data frame with a `genus` column (character).
#'
#' @return A sorted character vector of unique genus names.
#'
#' @importFrom dplyr filter pull
#'
#' @export
build_genus_vocabulary <- function(resolved_df) {
  resolved_df |>
    dplyr::filter(!is.na(genus), nchar(genus) >= 4) |>
    dplyr::pull(genus) |>
    unique() |>
    sort()
}

# ----------------------------------------

#' Build an epithet vocabulary from resolved taxa
#'
#' Extracts all unique specific epithets (>= 4 characters) from resolved
#' species names. Used to power token-level edit-distance correction in the
#' semi-automatic fuzzy pipeline.
#'
#' @param resolved_df A data frame with `genus` and `species` columns.
#'
#' @return A sorted character vector of unique specific epithets.
#'
#' @importFrom dplyr filter mutate pull
#' @importFrom stringr str_remove
#'
#' @export
build_epithet_vocabulary <- function(resolved_df) {
  resolved_df |>
    dplyr::filter(!is.na(species), !is.na(genus)) |>
    dplyr::mutate(
      epithet = stringr::str_remove(
        tolower(species),
        paste0("^", tolower(genus), "\\s+")
      )
    ) |>
    dplyr::filter(nchar(epithet) >= 4) |>
    dplyr::pull(epithet) |>
    unique() |>
    sort()
}

# ============================================================================
# WoRMS REST WRAPPERS (marine_only = false, extant_only = false)
# ============================================================================
# worrms does not expose extant_only (added to WoRMS REST API Feb 2025).
# These wrappers replace wm_records_name() / wm_records_taxamatch() wherever
# fossil-only taxa (e.g. Cyclococcolithina) must be included in searches.

#' Query WoRMS REST API by name
#'
#' Direct HTTP query to the WoRMS `/AphiaRecordsByName` endpoint with
#' `marine_only = false` and `extant_only = false`, bypassing the `worrms`
#' package which does not expose the `extant_only` parameter.
#'
#' @param name Scientific name to search.
#' @param fuzzy Logical. Whether to enable fuzzy (LIKE) matching. Default `TRUE`.
#' @param timeout_sec Request timeout in seconds. Default `10`.
#'
#' @return A [tibble::tibble()] with WoRMS record fields, or `NULL` if no
#'   results or the request fails.
#'
#' @importFrom httr modify_url GET timeout http_error content
#' @importFrom jsonlite fromJSON
#' @importFrom tibble as_tibble
#'
#' @keywords internal
worms_rest_byname <- function(name, fuzzy = TRUE, timeout_sec = 10) {
  url <- httr::modify_url(
    "https://www.marinespecies.org",
    path = paste0(
      "/rest/AphiaRecordsByName/",
      utils::URLencode(name, reserved = FALSE)
    ),
    query = list(
      like = "true",
      fuzzy = tolower(as.character(fuzzy)),
      marine_only = "false",
      extant_only = "false",
      offset = 1
    )
  )

  resp <- tryCatch(
    httr::GET(url, httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::http_error(resp)) {
    return(NULL)
  }

  content <- tryCatch(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    error = function(e) NULL
  )
  if (is.null(content) || content == "null" || content == "[]") {
    return(NULL)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(content, flatten = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || !is.data.frame(parsed) || nrow(parsed) == 0) {
    return(NULL)
  }
  tibble::as_tibble(parsed)
}

# ----------------------------------------

#' Query WoRMS REST API with Taxamatch
#'
#' Direct HTTP query to the WoRMS `/AphiaRecordsByMatchNames` endpoint with
#' `marine_only = false` and `extant_only = false`.
#'
#' @param name Scientific name to match.
#' @param timeout_sec Request timeout in seconds. Default `10`.
#'
#' @return A [tibble::tibble()] with WoRMS record fields, or `NULL`.
#'
#' @importFrom httr modify_url GET timeout http_error content
#' @importFrom jsonlite fromJSON
#' @importFrom tibble as_tibble
#'
#' @keywords internal
worms_rest_taxamatch <- function(name, timeout_sec = 10) {
  url <- httr::modify_url(
    "https://www.marinespecies.org",
    path = "/rest/AphiaRecordsByMatchNames",
    query = list(
      "scientificnames[]" = name,
      marine_only = "false",
      extant_only = "false"
    )
  )

  resp <- tryCatch(
    httr::GET(url, httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::http_error(resp)) {
    return(NULL)
  }

  content <- tryCatch(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    error = function(e) NULL
  )
  if (is.null(content) || content == "null" || content == "[[null]]") {
    return(NULL)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(content, flatten = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || length(parsed) == 0) {
    return(NULL)
  }

  result <- tryCatch(tibble::as_tibble(parsed[[1]]), error = function(e) NULL)
  if (is.null(result) || nrow(result) == 0) {
    return(NULL)
  }
  result
}

# ============================================================================
# NAME NORMALISATION
# ============================================================================

#' Normalise a taxonomic name before fuzzy search
#'
#' Handles two universal pre-search cases:
#' 1. **Vernacular names**: mapped to valid taxonomic names via the internal
#'    `.VERNACULAR_TO_TAXON` dictionary (case-insensitive).
#' 2. **Concatenated genus + epithet**: e.g. `"Chaetocerosdecipiens"` ->
#'    `"Chaetoceros decipiens"`, detected against the provided genus
#'    vocabulary.
#'
#' @param text Character. The taxon name to normalise.
#' @param genus_vocab Character vector of known genus names (from
#'   [build_genus_vocabulary()]).
#'
#' @return The normalised name, or the original if no rule applies.
#'
#' @keywords internal
normalize_taxonomic_name <- function(text, genus_vocab) {
  if (is.na(text)) {
    return(text)
  }

  # 1. Vernacular lookup
  text_lower <- tolower(text)
  if (text_lower %in% names(.vdict())) {
    return(.vdict()[[text_lower]])
  }

  # 2. Concatenated genus+epithet (no space, >= 10 chars)
  if (nchar(text) >= 10 && !grepl(" ", text)) {
    for (g in genus_vocab) {
      if (nchar(g) < 4) {
        next
      }
      pattern <- paste0("^(", g, ")([a-z].+)$")
      if (grepl(pattern, text, ignore.case = FALSE)) {
        return(sub(pattern, "\\1 \\2", text))
      }
    }
  }
  text
}

# ============================================================================
# EDIT-DISTANCE GENUS CORRECTION
# ============================================================================

#' Correct genus via edit-distance against a vocabulary
#'
#' Finds genus names in `vocab` within Levenshtein distance <= `max_dist`
#' of the input genus, builds full candidate names (genus + epithet), and
#' validates each against WoRMS via the REST API.
#'
#' @param taxon_clean Character. The cleaned taxon name to correct.
#' @param vocab Character vector of known genus names.
#' @param max_dist Integer. Maximum Levenshtein distance. Default `3`.
#' @param max_candidates Integer. Maximum genus candidates to test. Default `15`.
#' @param silent Logical. Suppress progress messages. Default `FALSE`.
#' @param timeout_sec Timeout for WoRMS REST calls. Default `10`.
#'
#' @return A [tibble::tibble()] with columns `variant`, `aphiaid`, `dist`,
#'   `method`, or `NULL` if no valid match found.
#'
#' @importFrom stringdist stringdist
#' @importFrom dplyr tibble arrange slice_head
#' @importFrom purrr map compact
#' @importFrom tibble tibble
#'
#' @keywords internal
search_by_edit_distance <- function(
  taxon_clean,
  vocab,
  max_dist = 3,
  max_candidates = 15,
  silent = FALSE,
  timeout_sec = 10
) {
  words <- strsplit(taxon_clean, " ")[[1]]
  genus <- words[1]
  epithet <- if (length(words) >= 2) {
    paste(words[-1], collapse = " ")
  } else {
    NA_character_
  }

  if (nchar(genus) < 4) {
    return(NULL)
  }

  dists <- stringdist::stringdist(genus, vocab, method = "lv")
  hits_idx <- which(dists > 0 & dists <= max_dist)
  if (length(hits_idx) == 0) {
    return(NULL)
  }

  hits_df <- dplyr::tibble(
    corrected_genus = vocab[hits_idx],
    dist = dists[hits_idx]
  ) |>
    dplyr::arrange(dist) |>
    dplyr::slice_head(n = max_candidates)

  if (!silent) {
    cat(sprintf(
      " -> Edit-distance: %d genus candidates (dist <= %d)\n",
      nrow(hits_df),
      max_dist
    ))
  }

  results <- purrr::map(hits_df$corrected_genus, function(cg) {
    candidate <- if (!is.na(epithet)) paste(cg, epithet) else cg

    fz_exact <- worms_with_timeout(
      worms_rest_byname(candidate, fuzzy = FALSE),
      timeout_sec = 5
    )
    if (!is.null(fz_exact) && is.data.frame(fz_exact) && nrow(fz_exact) > 0) {
      return(dplyr::tibble(
        variant = candidate,
        aphiaid = fz_exact$AphiaID[1],
        dist = hits_df$dist[hits_df$corrected_genus == cg][1],
        method = "edit_dist_exact"
      ))
    }

    tm <- worms_with_timeout(
      worms_rest_taxamatch(candidate),
      timeout_sec = 8
    )
    if (!is.null(tm) && is.data.frame(tm) && nrow(tm) > 0) {
      return(dplyr::tibble(
        variant = candidate,
        aphiaid = tm$AphiaID[1],
        dist = hits_df$dist[hits_df$corrected_genus == cg][1],
        method = "edit_dist_taxamatch"
      ))
    }
    NULL
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  if (nrow(results) == 0) {
    return(NULL)
  }
  results
}

# ----------------------------------------

#' Correct genus and epithet independently via token edit-distance
#'
#' Splits `taxon_clean` into genus + epithet tokens and corrects each
#' independently against the provided vocabularies. Builds all combinations
#' of corrected genus x corrected epithet and validates against WoRMS.
#'
#' @param taxon_clean Character. The cleaned binomial (or trinomial) name.
#' @param genus_vocab Character vector of known genus names.
#' @param epithet_vocab Character vector of known specific epithets.
#' @param max_dist Integer. Maximum per-token Levenshtein distance. Default `3`.
#' @param max_candidates Integer. Maximum epithet candidates to combine. Default `15`.
#' @param timeout_sec Timeout for WoRMS REST calls. Default `10`.
#'
#' @return A list of WoRMS record rows (tibbles), or `NULL`.
#'
#' @importFrom stringdist stringdist
#' @importFrom dplyr tibble arrange slice_head mutate bind_rows
#' @importFrom tidyr crossing
#'
#' @keywords internal
search_by_token_editdist <- function(
  taxon_clean,
  genus_vocab,
  epithet_vocab,
  max_dist = 3,
  max_candidates = 15,
  timeout_sec = 10
) {
  words <- strsplit(taxon_clean, " ")[[1]]
  n_words <- length(words)
  if (n_words < 2) {
    return(NULL)
  }

  genus <- words[1]
  epithet <- words[2]
  suffix <- if (n_words >= 3) {
    paste(words[-(1:2)], collapse = " ")
  } else {
    NA_character_
  }

  cat(sprintf(
    " -> L1c: token edit-distance (genus: '%s', epithet: '%s')...\n",
    genus,
    epithet
  ))

  # Genus candidates
  genus_dists <- stringdist::stringdist(
    tolower(genus),
    tolower(genus_vocab),
    method = "lv"
  )
  genus_exact_idx <- match(tolower(genus), tolower(genus_vocab))
  genus_exact <- !is.na(genus_exact_idx) && genus_dists[genus_exact_idx] == 0

  genus_candidates <- if (genus_exact) {
    dplyr::tibble(corrected = genus, dist = 0)
  } else {
    hits <- which(genus_dists <= max_dist)
    if (length(hits) == 0) {
      return(NULL)
    }
    dplyr::tibble(corrected = genus_vocab[hits], dist = genus_dists[hits]) |>
      dplyr::arrange(dist) |>
      dplyr::slice_head(n = 5)
  }

  # Epithet candidates
  epithet_dists <- stringdist::stringdist(
    tolower(epithet),
    epithet_vocab,
    method = "lv"
  )
  epithet_hits <- which(epithet_dists > 0 & epithet_dists <= max_dist)

  epithet_candidates <- if (length(epithet_hits) > 0) {
    dplyr::tibble(
      corrected = epithet_vocab[epithet_hits],
      dist = epithet_dists[epithet_hits]
    ) |>
      dplyr::arrange(dist) |>
      dplyr::slice_head(n = max_candidates)
  } else {
    dplyr::tibble(corrected = epithet, dist = 0)
  }

  cat(sprintf(
    " -> L1c: %d genus x %d epithet candidate(s)\n",
    nrow(genus_candidates),
    nrow(epithet_candidates)
  ))

  combos <- tidyr::crossing(
    g = genus_candidates$corrected,
    e = epithet_candidates$corrected
  ) |>
    dplyr::mutate(
      candidate = paste(g, e),
      candidate = if (!is.na(suffix)) paste(candidate, suffix) else candidate
    )

  results <- list()
  for (i in seq_len(nrow(combos))) {
    cand <- combos$candidate[i]
    if (tolower(cand) == tolower(taxon_clean)) {
      next
    }

    fz_exact <- worms_with_timeout(
      worms_rest_byname(cand, fuzzy = FALSE),
      timeout_sec = 5
    )
    if (!is.null(fz_exact) && is.data.frame(fz_exact) && nrow(fz_exact) > 0) {
      rec <- fz_exact[1, ]
      rec$ocr_corrected_from <- taxon_clean
      rec$ocr_variant <- cand
      if (!is.na(suffix)) {
        rec$is_binomial_match <- TRUE
        rec$binomial_suffix <- suffix
      }
      results <- c(results, list(rec))
      next
    }

    tm <- worms_with_timeout(worms_rest_taxamatch(cand), timeout_sec = 8)
    if (!is.null(tm) && is.data.frame(tm) && nrow(tm) > 0) {
      tm$ocr_corrected_from <- taxon_clean
      tm$ocr_variant <- cand
      if (!is.na(suffix)) {
        tm$is_binomial_match <- TRUE
        tm$binomial_suffix <- suffix
      }
      results <- c(results, split(tm, seq_len(nrow(tm))))
    }
  }

  if (length(results) == 0) {
    return(NULL)
  }
  results
}

# ============================================================================
# CASCADING SUGGESTION ORCHESTRATOR
# ============================================================================

#' Search WoRMS with cascading fuzzy strategies
#'
#' Orchestrates a multi-level cascade of WoRMS search strategies for a single
#' taxon name that could not be resolved automatically:
#'
#' - **L1**: Full name, taxamatch + REST fuzzy
#' - **L1b**: Edit-distance genus correction (local vocabulary)
#' - **L1c**: Token-level edit-distance genus + epithet (local vocabulary)
#' - **L1d**: WoRMS fuzzy genus reconstruction (remote, genera not in vocab)
#' - **L1e**: Dehyphenated epithet
#' - **L2**: Truncate last token stem
#' - **L3**: Truncate last + second token (trinomials)
#' - **L4**: Drop suffix, search binomial (trinomials)
#' - **L5**: Drop suffix + truncate epithet (trinomials)
#' - **L6**: Genus only (always, shown as LTR fallback)
#'
#' Returns consolidated, ranked suggestions grouped by confidence tier.
#'
#' @param taxon_clean Character. Cleaned taxon name to search.
#' @param original_taxon Character. Original verbatim name (for display).
#' @param genus_vocab Character vector from [build_genus_vocabulary()].
#' @param epithet_vocab Character vector from [build_epithet_vocabulary()].
#' @param min_similarity Minimum similarity score for "good" suggestions.
#'   Default `0.85`.
#' @param max_suggestions Maximum number of good suggestions to return.
#'   Default `15`.
#' @param max_longshot_suggestions Maximum number of longshot suggestions.
#'   Default `10`.
#' @param longshot_threshold Minimum similarity for longshot bucket. Default `0.30`.
#' @param edit_max_dist Maximum Levenshtein distance for edit-distance steps.
#'   Default `3`.
#' @param edit_max_candidates Maximum genus candidates in edit-distance steps.
#'   Default `15`.
#' @param timeout_sec Timeout for WoRMS REST calls in seconds. Default `15`.
#'
#' @return A named list with elements `good`, `longshot`, `genus_fallback`,
#'   and `binomial_fallback` (each a tibble or `NULL`).
#'
#' @importFrom dplyr filter mutate arrange bind_rows distinct if_else
#' @importFrom stringr str_split str_count str_detect str_replace_all
#' @importFrom purrr map_dbl
#' @importFrom worrms wm_record
#'
#' @export
search_worms_fuzzy_suggestions <- function(
  taxon_clean,
  original_taxon,
  genus_vocab,
  epithet_vocab,
  min_similarity = 0.85,
  max_suggestions = 15,
  max_longshot_suggestions = 10,
  longshot_threshold = 0.30,
  edit_max_dist = 3,
  edit_max_candidates = 15,
  timeout_sec = 15
) {
  if (is.na(taxon_clean) || trimws(taxon_clean) == "") {
    return(list(
      good = NULL,
      longshot = NULL,
      genus_fallback = NULL,
      binomial_fallback = NULL
    ))
  }

  # Pre-normalise
  normalized <- normalize_taxonomic_name(taxon_clean, genus_vocab)
  if (normalized != taxon_clean) {
    cat(sprintf(" -> Normalized: '%s' -> '%s'\n", taxon_clean, normalized))
    taxon_clean <- normalized
  }

  cat(sprintf("\nSearching WoRMS for '%s'...\n", taxon_clean))

  words <- strsplit(taxon_clean, " ")[[1]]
  n_words <- length(words)
  genus <- words[1]
  suffix <- if (n_words >= 3) {
    paste(words[-(1:2)], collapse = " ")
  } else {
    NA_character_
  }

  all_records <- list()

  stem_of <- function(token) {
    substr(token, 1, max(5L, floor(nchar(token) * 0.7)))
  }

  worms_query <- function(
    query,
    tag_binomial = FALSE,
    tag_suffix = NA_character_
  ) {
    recs <- list()
    tm <- worms_with_timeout(worms_rest_taxamatch(query), timeout_sec = 8)
    if (!is.null(tm) && is.data.frame(tm) && nrow(tm) > 0) {
      recs <- c(recs, split(tm, seq_len(nrow(tm))))
    }
    fz <- worms_with_timeout(
      worms_rest_byname(query, fuzzy = TRUE),
      timeout_sec = 8
    )
    if (!is.null(fz) && is.data.frame(fz) && nrow(fz) > 0) {
      recs <- c(recs, split(fz, seq_len(nrow(fz))))
    }
    if (length(recs) == 0) {
      return(NULL)
    }
    out <- dplyr::bind_rows(recs) |>
      dplyr::filter(
        !is.na(AphiaID),
        !is.na(scientificname),
        scientificname != ""
      )
    if (tag_binomial && !is.na(tag_suffix)) {
      out <- out |>
        dplyr::mutate(
          query_tokens = stringr::str_count(trimws(query), "\\S+"),
          returned_tokens = stringr::str_count(trimws(scientificname), "\\S+"),
          is_binomial_match = returned_tokens <= query_tokens,
          binomial_suffix = dplyr::if_else(
            is_binomial_match,
            tag_suffix,
            NA_character_
          )
        ) |>
        dplyr::select(-query_tokens, -returned_tokens)
    }
    if (nrow(out) == 0) {
      return(NULL)
    }
    split(out, seq_len(nrow(out)))
  }

  has_good <- function(recs) {
    if (length(recs) == 0) {
      return(FALSE)
    }
    combined <- dplyr::bind_rows(recs) |>
      dplyr::filter(
        !is.na(AphiaID),
        !is.na(scientificname),
        scientificname != ""
      ) |>
      dplyr::distinct(AphiaID, .keep_all = TRUE)
    if (nrow(combined) == 0) {
      return(FALSE)
    }
    sims <- purrr::map_dbl(
      combined$scientificname,
      ~ calculate_similarity(.x, taxon_clean)
    )
    any(sims >= min_similarity, na.rm = TRUE)
  }

  cat(" -> L1: full name...\n")
  l1 <- worms_query(taxon_clean)
  if (!is.null(l1)) {
    all_records <- c(all_records, l1)
  }

  if (n_words >= 2 && length(genus_vocab) > 0) {
    cat(sprintf(" -> L1b: edit-distance genus (input: '%s')...\n", genus))
    ed <- search_by_edit_distance(
      taxon_clean,
      genus_vocab,
      max_dist = edit_max_dist,
      max_candidates = edit_max_candidates,
      silent = TRUE,
      timeout_sec = timeout_sec
    )
    if (!is.null(ed) && nrow(ed) > 0) {
      cat(sprintf(" -> L1b: %d match(es)\n", nrow(ed)))
      for (i in seq_len(nrow(ed))) {
        rec <- worms_with_timeout(worrms::wm_record(ed$aphiaid[i]))
        if (!is.null(rec)) {
          rec$ocr_corrected_from <- taxon_clean
          rec$ocr_variant <- ed$variant[i]
          all_records <- c(all_records, list(rec))
        }
      }
    }
  }

  if (
    !has_good(all_records) &&
      n_words >= 2 &&
      length(genus_vocab) > 0 &&
      length(epithet_vocab) > 0
  ) {
    lc <- search_by_token_editdist(
      taxon_clean,
      genus_vocab,
      epithet_vocab,
      max_dist = edit_max_dist,
      max_candidates = edit_max_candidates,
      timeout_sec = timeout_sec
    )
    if (!is.null(lc)) all_records <- c(all_records, lc)
  }

  if (!has_good(all_records) && n_words >= 2) {
    cat(sprintf(" -> L1d: WoRMS fuzzy genus lookup '%s'...\n", genus))
    genus_fz <- worms_with_timeout(
      worms_rest_byname(genus, fuzzy = TRUE),
      timeout_sec = 8
    )
    if (!is.null(genus_fz) && is.data.frame(genus_fz) && nrow(genus_fz) > 0) {
      genus_hits <- genus_fz |>
        dplyr::filter(
          tolower(rank) == "genus",
          tolower(scientificname) != tolower(genus)
        ) |>
        dplyr::distinct(scientificname) |>
        dplyr::slice_head(n = 5)

      if (nrow(genus_hits) > 0) {
        cat(sprintf(" -> L1d: %d genus candidate(s)\n", nrow(genus_hits)))
        for (g in genus_hits$scientificname) {
          candidate <- paste(g, paste(words[-1], collapse = " "))
          cat(sprintf(" -> L1d: trying '%s'...\n", candidate))
          l1d_recs <- worms_query(candidate)
          if (!is.null(l1d_recs)) {
            l1d_tagged <- dplyr::bind_rows(l1d_recs) |>
              dplyr::mutate(
                ocr_corrected_from = taxon_clean,
                ocr_variant = candidate
              )
            all_records <- c(
              all_records,
              split(l1d_tagged, seq_len(nrow(l1d_tagged)))
            )
          }
        }
      }
    }
  }

  if (stringr::str_detect(taxon_clean, "\\w+-\\w+")) {
    dehyphenated <- stringr::str_replace_all(taxon_clean, "-(\\w)", " \\1")
    if (dehyphenated != taxon_clean) {
      cat(sprintf(" -> L1e: dehyphenated: '%s'...\n", dehyphenated))
      l1e <- worms_query(dehyphenated)
      if (!is.null(l1e)) all_records <- c(all_records, l1e)
    }
  }

  if (!has_good(all_records)) {
    last_stem <- stem_of(words[n_words])
    words_l2 <- replace(words, n_words, last_stem)
    query_l2 <- paste(words_l2, collapse = " ")
    if (query_l2 != taxon_clean) {
      cat(sprintf(" -> L2: truncate last token: '%s'...\n", query_l2))
      l2 <- worms_query(
        query_l2,
        tag_binomial = !is.na(suffix),
        tag_suffix = suffix
      )
      if (!is.null(l2)) all_records <- c(all_records, l2)
    }

    if (!has_good(all_records) && n_words >= 3) {
      second_stem <- stem_of(words[2])
      words_l3 <- replace(words_l2, 2, second_stem)
      query_l3 <- paste(words_l3, collapse = " ")
      if (query_l3 != query_l2) {
        cat(sprintf(" -> L3: truncate 2nd + last token: '%s'...\n", query_l3))
        l3 <- worms_query(query_l3, tag_binomial = TRUE, tag_suffix = suffix)
        if (!is.null(l3)) all_records <- c(all_records, l3)
      }

      if (!has_good(all_records) && !is.na(suffix)) {
        binomial <- paste(words[1:2], collapse = " ")
        cat(sprintf(
          " -> L4: binomial '%s' [suffix -> note: '%s']...\n",
          binomial,
          suffix
        ))
        l4 <- worms_query(binomial, tag_binomial = TRUE, tag_suffix = suffix)
        if (!is.null(l4)) {
          all_records <- c(all_records, l4)
        }

        if (!has_good(all_records)) {
          query_l5 <- paste(genus, stem_of(words[2]))
          cat(sprintf(
            " -> L5: stem binomial '%s' [suffix -> note: '%s']...\n",
            query_l5,
            suffix
          ))
          l5 <- worms_query(query_l5, tag_binomial = TRUE, tag_suffix = suffix)
          if (!is.null(l5)) all_records <- c(all_records, l5)
        }
      }
    }
  }

  cat(sprintf(" -> L6: genus '%s'...\n", genus))
  gr <- worms_with_timeout(
    worms_rest_byname(genus, fuzzy = TRUE),
    timeout_sec = 8
  )
  if (!is.null(gr) && is.data.frame(gr) && nrow(gr) > 0) {
    all_records <- c(all_records, split(gr, seq_len(nrow(gr))))
  }

  if (length(all_records) == 0) {
    cat(" No matches\n")
    return(list(
      good = NULL,
      longshot = NULL,
      genus_fallback = NULL,
      binomial_fallback = NULL
    ))
  }

  cat(sprintf(" -> Consolidating %d records...\n", length(all_records)))

  combined <- dplyr::bind_rows(all_records)
  if (!"ocr_variant" %in% names(combined)) {
    combined$ocr_variant <- NA_character_
  }
  if (!"ocr_corrected_from" %in% names(combined)) {
    combined$ocr_corrected_from <- NA_character_
  }
  if (!"is_binomial_match" %in% names(combined)) {
    combined$is_binomial_match <- FALSE
  }
  if (!"binomial_suffix" %in% names(combined)) {
    combined$binomial_suffix <- NA_character_
  }

  unique_records <- combined |>
    dplyr::filter(
      !is.na(AphiaID),
      !is.na(scientificname),
      scientificname != ""
    ) |>
    dplyr::distinct(AphiaID, .keep_all = TRUE) |>
    dplyr::mutate(
      is_ocr_match = !is.na(ocr_variant),
      is_binomial_match = dplyr::if_else(
        is.na(is_binomial_match),
        FALSE,
        is_binomial_match
      ),
      binomial_suffix = dplyr::if_else(
        is.na(binomial_suffix),
        NA_character_,
        binomial_suffix
      )
    )

  cat(sprintf(" -> %d unique records\n", nrow(unique_records)))

  unique_records <- unique_records |>
    dplyr::mutate(
      similarity = purrr::map_dbl(
        scientificname,
        ~ calculate_similarity(.x, taxon_clean)
      ),
      adjusted_similarity = dplyr::case_when(
        is_ocr_match ~ pmin(similarity + 0.15, 1.0),
        is_binomial_match ~ pmin(similarity + 0.10, 1.0),
        TRUE ~ similarity
      ),
      epithet_jw = purrr::map_dbl(scientificname, function(nm) {
        w1 <- strsplit(tolower(nm), "[\\s-]+")[[1]]
        w2 <- strsplit(tolower(taxon_clean), "[\\s-]+")[[1]]
        if (length(w1) >= 2 && length(w2) >= 2) {
          1 - stringdist::stringdist(w1[2], w2[2], method = "jw")
        } else {
          0
        }
      })
    ) |>
    dplyr::arrange(
      dplyr::desc(adjusted_similarity),
      dplyr::desc(epithet_jw),
      status
    )

  # Genus/higher-rank fallback
  SUPRA_RANKS <- c(
    "genus",
    "subgenus",
    "family",
    "subfamily",
    "tribe",
    "superfamily",
    "order",
    "suborder",
    "class",
    "subclass",
    "phylum",
    "kingdom"
  )

  genus_fallback <- NULL
  gfb <- unique_records |>
    dplyr::filter(
      tolower(rank) %in% SUPRA_RANKS,
      tolower(trimws(scientificname)) == tolower(trimws(genus))
    ) |>
    dplyr::arrange(dplyr::desc(adjusted_similarity)) |>
    dplyr::slice_head(n = 1)
  if (nrow(gfb) > 0) {
    genus_fallback <- gfb
  }

  fallback_ids <- if (!is.null(genus_fallback)) {
    genus_fallback$AphiaID
  } else {
    integer(0)
  }

  good <- unique_records |>
    dplyr::filter(adjusted_similarity >= min_similarity) |>
    dplyr::slice_head(n = max_suggestions)

  longshot <- unique_records |>
    dplyr::filter(
      adjusted_similarity >= longshot_threshold,
      adjusted_similarity < min_similarity,
      !(AphiaID %in% fallback_ids)
    ) |>
    dplyr::slice_head(n = max_longshot_suggestions)

  if (nrow(good) == 0 && nrow(longshot) == 0) {
    cat(" Below threshold -- showing all candidates\n")
    longshot <- unique_records |>
      dplyr::filter(tolower(rank) != "genus", !(AphiaID %in% fallback_ids)) |>
      dplyr::slice_head(n = max_longshot_suggestions)
  }

  # Binomial parent fallback (trinomials only)
  binomial_fallback <- NULL
  if (!is.na(suffix)) {
    binomial_exact <- paste(words[1:2], collapse = " ")
    bfb <- unique_records |>
      dplyr::filter(
        tolower(rank) == "species",
        grepl(
          paste0("^", binomial_exact, "$"),
          scientificname,
          ignore.case = TRUE
        )
      ) |>
      dplyr::slice_head(n = 1)
    if (nrow(bfb) > 0) binomial_fallback <- bfb
  }

  list(
    good = if (nrow(good) > 0) good else NULL,
    longshot = if (nrow(longshot) > 0) longshot else NULL,
    genus_fallback = genus_fallback,
    binomial_fallback = binomial_fallback
  )
}

# ============================================================================
# INTERACTIVE PROMPTS
# ============================================================================

#' Present fuzzy suggestions to the user and collect a choice
#'
#' Displays ranked WoRMS candidates grouped by tier (good, binomial parent,
#' longshot, genus fallback) and prompts the user to accept one, enter an
#' AphiaID manually, flag for removal, request expert review, or skip.
#'
#' @param suggestions Named list returned by [search_worms_fuzzy_suggestions()].
#' @param original_taxon Character. Original verbatim taxon name.
#' @param taxon_clean Character. Cleaned taxon name.
#'
#' @return A named list with `action` (one of `"accept"`, `"skip"`,
#'   `"cancel"`, `"expert_review"`, `"flag_for_removal"`) plus, when
#'   `action == "accept"`, fields `aphiaid`, `type`, `name`, and
#'   `binomial_suffix`.
#'
#' @export
prompt_fuzzy_suggestions <- function(suggestions, original_taxon, taxon_clean) {
  repeat {
    cat("\n")
    cat(strrep("=", 70), "\n", sep = "")
    cat(sprintf("Original: '%s'\n", original_taxon))
    if (original_taxon != taxon_clean) {
      cat(sprintf("Cleaned:  '%s'\n", taxon_clean))
    }
    cat(strrep("-", 70), "\n", sep = "")

    options <- list()
    counter <- 1
    has_suggestions <- FALSE

    # Good suggestions
    if (!is.null(suggestions$good)) {
      has_suggestions <- TRUE
      cat("Good suggestions:\n")
      for (i in seq_len(nrow(suggestions$good))) {
        row <- suggestions$good[i, ]
        ocr_tag <- if (!is.na(row$is_ocr_match) && row$is_ocr_match) {
          sprintf(" [Edit-dist -> '%s']", row$ocr_variant)
        } else {
          ""
        }
        bin_tag <- if (!is.na(row$is_binomial_match) && row$is_binomial_match) {
          sprintf(" [Binomial -- add '%s' as note]", row$binomial_suffix)
        } else {
          ""
        }
        cat(sprintf(
          " [%d] %s\n     Status: %-12s Sim: %.2f AphiaID: %d%s%s\n",
          counter,
          row$scientificname,
          row$status,
          row$adjusted_similarity,
          row$AphiaID,
          ocr_tag,
          bin_tag
        ))
        options[[as.character(counter)]] <- list(
          aphiaid = row$AphiaID,
          type = "good",
          name = row$scientificname,
          binomial_suffix = if (
            !is.na(row$is_binomial_match) && row$is_binomial_match
          ) {
            row$binomial_suffix
          } else {
            NA_character_
          }
        )
        counter <- counter + 1
      }
    }

    # Binomial parent
    if (!is.null(suggestions$binomial_fallback)) {
      has_suggestions <- TRUE
      row <- suggestions$binomial_fallback[1, ]
      w <- strsplit(taxon_clean, " ")[[1]]
      stored_suffix <- if (!is.na(row$binomial_suffix)) {
        row$binomial_suffix
      } else if (length(w) >= 3) {
        paste(w[-(1:2)], collapse = " ")
      } else {
        NA_character_
      }
      sfx_tag <- if (!is.na(stored_suffix)) {
        sprintf(" [Parent -- add '%s' as note]", stored_suffix)
      } else {
        ""
      }
      cat("Binomial parent:\n")
      cat(sprintf(
        " [%d] %s\n     Status: %-12s AphiaID: %d%s\n",
        counter,
        row$scientificname,
        row$status,
        row$AphiaID,
        sfx_tag
      ))
      options[[as.character(counter)]] <- list(
        aphiaid = row$AphiaID,
        type = "good",
        name = row$scientificname,
        binomial_suffix = stored_suffix
      )
      counter <- counter + 1
    }

    # Longshot
    if (!is.null(suggestions$longshot)) {
      has_suggestions <- TRUE
      cat("Long-shot suggestions:\n")
      for (i in seq_len(nrow(suggestions$longshot))) {
        row <- suggestions$longshot[i, ]
        bin_tag <- if (!is.na(row$is_binomial_match) && row$is_binomial_match) {
          sprintf(" [Binomial -- add '%s' as note]", row$binomial_suffix)
        } else {
          ""
        }
        cat(sprintf(
          " [%d] %s\n     Status: %-12s Sim: %.2f AphiaID: %d%s\n",
          counter,
          row$scientificname,
          row$status,
          row$adjusted_similarity,
          row$AphiaID,
          bin_tag
        ))
        options[[as.character(counter)]] <- list(
          aphiaid = row$AphiaID,
          type = "longshot",
          name = row$scientificname,
          binomial_suffix = if (
            !is.na(row$is_binomial_match) && row$is_binomial_match
          ) {
            row$binomial_suffix
          } else {
            NA_character_
          }
        )
        counter <- counter + 1
      }
    }

    # Genus fallback
    if (!is.null(suggestions$genus_fallback)) {
      has_suggestions <- TRUE
      row <- suggestions$genus_fallback[1, ]
      cat("Genus-level fallback:\n")
      cat(sprintf(
        " [%d] %s\n     Status: %-12s AphiaID: %d [LTR]\n",
        counter,
        row$scientificname,
        row$status,
        row$AphiaID
      ))
      options[[as.character(counter)]] <- list(
        aphiaid = row$AphiaID,
        type = "genus_fallback",
        name = row$scientificname,
        binomial_suffix = NA_character_
      )
      counter <- counter + 1
    }

    if (!has_suggestions) {
      cat("No suggestions found.\n")
    }

    cat(strrep("-", 70), "\n", sep = "")
    cat(
      "Options: [number] accept | [m]anual | [e]xpert review | [f]lag for removal | [s]kip | [c]ancel\n"
    )
    cat("Choice: ")
    choice <- tolower(trimws(readline()))

    if (choice == "c") {
      return(list(action = "cancel"))
    } else if (choice == "s") {
      return(list(action = "skip"))
    } else if (choice == "e") {
      return(list(action = "expert_review"))
    } else if (choice == "f") {
      return(list(action = "flag_for_removal"))
    } else if (choice == "m") {
      result <- prompt_manual_aphiaid()
      if (result$action != "retry") return(result)
    } else if (choice %in% names(options)) {
      sel <- options[[choice]]
      return(list(
        action = "accept",
        aphiaid = sel$aphiaid,
        type = sel$type,
        name = sel$name,
        binomial_suffix = sel$binomial_suffix
      ))
    } else {
      cat(" Invalid choice. Please try again.\n")
      Sys.sleep(0.5)
    }
  }
}

# ----------------------------------------

#' Prompt for a manual AphiaID or name entry
#'
#' Allows the user to enter an AphiaID (integer) or a scientific name directly.
#' Supports fuzzy fallback when exact name lookup fails.
#'
#' @return A named list with `action` (`"accept"`, `"skip"`, or `"retry"`)
#'   plus `aphiaid` and `type` when accepted.
#'
#' @importFrom worrms wm_record wm_name2id wm_records_name
#'
#' @export
prompt_manual_aphiaid <- function() {
  repeat {
    cat("Enter AphiaID (integer) or scientific name [s = skip]: ")
    input <- trimws(readline())

    if (tolower(input) == "s") {
      return(list(action = "skip"))
    }

    aphia_id <- suppressWarnings(as.integer(input))
    if (!is.na(aphia_id)) {
      cat(sprintf(" Fetching WoRMS record for AphiaID %d...\n", aphia_id))
      record <- worms_with_timeout(worrms::wm_record(aphia_id))
      if (!is.null(record)) {
        cat(sprintf(
          " Found: %s | Status: %s | AphiaID: %d\n",
          record$scientificname,
          record$status,
          record$AphiaID
        ))
        cat(" Accept? [y/n]: ")
        if (tolower(trimws(readline())) == "y") {
          return(list(
            action = "accept",
            aphiaid = aphia_id,
            type = "manual",
            binomial_suffix = NA_character_
          ))
        } else {
          cat(" Cancelled. Try again.\n")
          next
        }
      } else {
        cat(sprintf(
          " No WoRMS record found for AphiaID %d. Try again.\n",
          aphia_id
        ))
        next
      }
    }

    # Name lookup
    cat(sprintf(" Searching WoRMS for '%s'...\n", input))
    result_id <- worms_with_timeout(worrms::wm_name2id(input))
    if (!is.null(result_id) && !is.na(result_id)) {
      record <- worms_with_timeout(worrms::wm_record(result_id))
      if (!is.null(record)) {
        cat(sprintf(
          " Found: %s | Status: %s | AphiaID: %d\n",
          record$scientificname,
          record$status,
          record$AphiaID
        ))
        cat(" Accept? [y/n]: ")
        if (tolower(trimws(readline())) == "y") {
          return(list(
            action = "accept",
            aphiaid = result_id,
            type = "manual",
            binomial_suffix = NA_character_
          ))
        } else {
          cat(" Cancelled. Try again.\n")
        }
        next
      }
    }

    # Fuzzy fallback
    cat(" Exact match not found. Trying fuzzy...\n")
    fuzzy_r <- worms_with_timeout(
      worrms::wm_records_name(input, fuzzy = TRUE, marine_only = FALSE)
    )
    if (!is.null(fuzzy_r) && nrow(fuzzy_r) > 0) {
      n_show <- min(nrow(fuzzy_r), 5)
      for (i in seq_len(n_show)) {
        cat(sprintf(
          " [%d] %s | %s | AphiaID: %d\n",
          i,
          fuzzy_r$scientificname[i],
          fuzzy_r$status[i],
          fuzzy_r$AphiaID[i]
        ))
      }
      cat(sprintf(" Select [1-%d] or [s]kip: ", n_show))
      sel <- trimws(readline())
      if (tolower(sel) == "s") {
        return(list(action = "retry"))
      }
      sel_n <- suppressWarnings(as.integer(sel))
      if (!is.na(sel_n) && sel_n >= 1 && sel_n <= n_show) {
        return(list(
          action = "accept",
          aphiaid = fuzzy_r$AphiaID[sel_n],
          type = "manual",
          binomial_suffix = NA_character_
        ))
      } else {
        cat(" Invalid choice.\n")
      }
    } else {
      cat(" No matches found. Try again.\n")
    }
  }
}

# ----------------------------------------

#' Prompt for resolution notes after accepting a match
#'
#' Shows the matched and accepted names, then asks the user for optional
#' free-text notes. Returns `"retry"` if the user wants to undo the selection.
#'
#' @param matched_name Character. The name as matched in WoRMS.
#' @param accepted_name Character. The currently accepted name.
#'
#' @return A named list with `action` (`"ok"` or `"retry"`) and `notes`
#'   (character or `NA_character_`).
#'
#' @export
prompt_resolution_notes <- function(matched_name, accepted_name) {
  cat(sprintf(" Matched: '%s' -> '%s'\n", matched_name, accepted_name))
  cat(" Add notes [Enter to skip] or [r]etry: ")
  notes <- trimws(readline())
  if (tolower(notes) == "r") {
    return(list(action = "retry"))
  }
  list(action = "ok", notes = if (notes == "") NA_character_ else notes)
}

# ============================================================================
# WoRMS RECORD RETRIEVAL AND APPLICATION
# ============================================================================

#' Fetch a full WoRMS record including classification
#'
#' Given an AphiaID, retrieves the matched record, resolves it to the accepted
#' name if it is a synonym, and returns the full taxonomic classification.
#'
#' @param aphia_id Integer. The AphiaID to fetch.
#' @param timeout_sec Timeout in seconds. Default `15`.
#'
#' @return A named list with `aphiaid`, `valid_aphiaid`, `matched_name`,
#'   `accepted_name`, `status`, `rank`, and `taxonomy` (a named character
#'   vector of all rank columns), or `NULL` if the record cannot be fetched.
#'
#' @importFrom worrms wm_record wm_classification
#'
#' @export
fetch_worms_full_record <- function(aphia_id, timeout_sec = 15) {
  if (is.na(aphia_id) || is.null(aphia_id)) {
    return(NULL)
  }

  record <- worms_with_timeout(worrms::wm_record(aphia_id), timeout_sec)
  if (is.null(record)) {
    return(NULL)
  }

  target_id <- if (
    !is.null(record$valid_AphiaID) &&
      !is.na(record$valid_AphiaID) &&
      record$AphiaID != record$valid_AphiaID
  ) {
    record$valid_AphiaID
  } else {
    aphia_id
  }

  classification <- worms_with_timeout(
    worrms::wm_classification(target_id),
    timeout_sec
  )
  accepted_record <- if (target_id != aphia_id) {
    worms_with_timeout(worrms::wm_record(target_id), timeout_sec)
  } else {
    record
  }

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
  taxonomy <- stats::setNames(rep(NA_character_, length(tax_ranks)), tax_ranks)
  if (!is.null(classification) && nrow(classification) > 0) {
    for (j in seq_len(nrow(classification))) {
      rank_lc <- tolower(classification$rank[j])
      if (rank_lc %in% tax_ranks) {
        taxonomy[[rank_lc]] <- classification$scientificname[j]
      }
    }
  }

  list(
    aphiaid = aphia_id,
    valid_aphiaid = if (target_id != aphia_id) {
      as.integer(target_id)
    } else {
      NA_integer_
    },
    matched_name = record$scientificname,
    accepted_name = if (!is.null(accepted_record)) {
      accepted_record$scientificname
    } else {
      record$scientificname
    },
    status = record$status,
    taxonomy = taxonomy,
    rank = if (!is.null(accepted_record)) accepted_record$rank else record$rank
  )
}

# ----------------------------------------

#' Apply a WoRMS resolution result to a data frame row
#'
#' Writes all resolution fields (IDs, names, status, method, full taxonomy)
#' into the specified row of the working data frame.
#'
#' @param df A data frame with the resolution and taxonomy columns.
#' @param idx Integer. Row index to update.
#' @param worms_data Named list returned by [fetch_worms_full_record()].
#' @param resolution_type Character. One of `"good"`, `"longshot"`,
#'   `"genus_fallback"`, or `"manual"`.
#'
#' @return The updated data frame.
#'
#' @export
apply_resolution <- function(df, idx, worms_data, resolution_type) {
  method_name <- switch(
    resolution_type,
    good = "manual_worms_fuzzy",
    longshot = "manual_worms_longshot",
    genus_fallback = "manual_worms_genus_ltr",
    manual = "manual_worms_direct",
    "manual_worms_fuzzy"
  )

  df$matched_aphiaid[idx] <- as.integer(worms_data$aphiaid)
  df$matched_name[idx] <- worms_data$matched_name
  df$accepted_name[idx] <- worms_data$accepted_name
  df$taxonomic_status[idx] <- worms_data$status
  df$resolution_method[idx] <- method_name
  df$rank[idx] <- worms_data$rank

  if (!is.na(worms_data$valid_aphiaid)) {
    df$accepted_aphiaid[idx] <- worms_data$valid_aphiaid
  }

  for (tax_rank in names(worms_data$taxonomy)) {
    if (tax_rank %in% names(df) && !is.na(worms_data$taxonomy[[tax_rank]])) {
      df[[tax_rank]][idx] <- worms_data$taxonomy[[tax_rank]]
    }
  }
  df
}

# ----------------------------------------

#' Save fuzzy resolution progress to disk
#'
#' Saves the current state of `df` (and optionally the genus/epithet
#' vocabularies) to `checkpoint_file` as an RDS file. Progress can be
#' resumed in a new session by passing the same `checkpoint_file` path
#' to [process_fuzzy_batch()].
#'
#' @param df The working data frame to save.
#' @param checkpoint_file Path to the RDS file. Set to `NULL` to disable.
#'
#' @return `TRUE` invisibly on success, `FALSE` on error.
#'
#' @keywords internal
save_progress <- function(
  df,
  checkpoint_file = NULL,
  genus_vocab = NULL,
  epithet_vocab = NULL
) {
  tryCatch(
    {
      if (!is.null(checkpoint_file)) {
        saveRDS(
          list(
            df = df,
            checkpoint_file = checkpoint_file,
            genus_vocab = genus_vocab,
            epithet_vocab = epithet_vocab
          ),
          file = checkpoint_file
        )
      }
      invisible(TRUE)
    },
    error = function(e) {
      cat(sprintf(" Failed to save progress: %s\n", e$message))
      invisible(FALSE)
    }
  )
}

# ============================================================================
# BATCH PROCESSING LOOP
# ============================================================================

#' Run the interactive semi-automatic fuzzy resolution loop
#'
#' Iterates over all unresolved entries in `df` (where `matched_aphiaid` is
#' `NA` and `resolution_method` is not `"expert_review"`), presents cascading
#' WoRMS suggestions via [search_worms_fuzzy_suggestions()], collects user
#' choices via [prompt_fuzzy_suggestions()], and writes accepted resolutions
#' back to `df`. Progress is saved to `.GlobalEnv` after each accepted entry.
#'
#' @param df A data frame with the resolution schema (as produced by the
#'   cleaning + automatic resolution pipeline).
#' @param genus_vocab Character vector from [build_genus_vocabulary()].
#' @param epithet_vocab Character vector from [build_epithet_vocabulary()].
#' @param batch_size Integer. Number of entries to process before pausing for
#'   continuation confirmation. Default `10`.
#' @param min_similarity Minimum similarity for "good" suggestions. Default `0.85`.
#' @param max_suggestions Maximum good suggestions to display. Default `15`.
#' @param max_longshot_suggestions Maximum longshot suggestions. Default `10`.
#' @param longshot_threshold Minimum similarity for longshot bucket. Default `0.30`.
#' @param edit_max_dist Maximum Levenshtein distance for edit-distance steps.
#'   Default `3`.
#' @param edit_max_candidates Maximum genus candidates in edit-distance steps.
#'   Default `15`.
#' @param timeout_sec Timeout for WoRMS calls in seconds. Default `15`.
#' @param checkpoint_file Character. Path to an RDS file used to save and
#'   resume progress between sessions. After each accepted resolution the
#'   current state is written to this file. If the file already exists when
#'   the function is called, it is loaded automatically and processing
#'   resumes from where it left off. Set to `NULL` to disable checkpointing.
#'   Default `"phytaxr_step6_checkpoint.rds"`.
#'
#' @return The updated data frame.
#'
#' @importFrom stringr str_split
#'
#' @export
process_fuzzy_batch <- function(
  df,
  genus_vocab,
  epithet_vocab,
  batch_size = 10,
  checkpoint_file = "phytaxr_step6_checkpoint.rds",
  min_similarity = 0.85,
  max_suggestions = 15,
  max_longshot_suggestions = 10,
  longshot_threshold = 0.30,
  edit_max_dist = 3,
  edit_max_candidates = 15,
  timeout_sec = 15
) {
  df <- ensure_resolution_schema(df)
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    cat(sprintf("Checkpoint found: %s -- resuming.\n", checkpoint_file))
    ckpt <- readRDS(checkpoint_file)
    df <- ckpt$df
    if (!is.null(ckpt$genus_vocab)) {
      genus_vocab <- ckpt$genus_vocab
    }
    if (!is.null(ckpt$epithet_vocab)) epithet_vocab <- ckpt$epithet_vocab
  }

  unresolved_idx <- which(
    is.na(df$matched_aphiaid) &
      (is.na(df$resolution_method) | df$resolution_method != "expert_review")
  )
  total_unresolved <- length(unresolved_idx)

  if (total_unresolved == 0) {
    cat("No unresolved entries.\n")
    return(df)
  }

  cat(sprintf(
    "Unresolved entries: %d | Batch size: %d\n\n",
    total_unresolved,
    batch_size
  ))

  batch_start <- 1

  while (batch_start <= total_unresolved) {
    batch_end <- min(batch_start + batch_size - 1, total_unresolved)
    batch_indices <- unresolved_idx[batch_start:batch_end]

    cat(sprintf(
      "\n+-- BATCH %d-%d of %d --+\n",
      batch_start,
      batch_end,
      total_unresolved
    ))

    for (i in seq_along(batch_indices)) {
      idx <- batch_indices[i]
      original <- df$taxon[idx]
      cleaned <- df$taxon_clean[idx]

      cat(sprintf(
        "\nEntry %d of %d\n%s\n",
        batch_start + i - 1,
        total_unresolved,
        strrep("=", 70)
      ))

      suggestions <- search_worms_fuzzy_suggestions(
        cleaned,
        original,
        genus_vocab = genus_vocab,
        epithet_vocab = epithet_vocab,
        min_similarity = min_similarity,
        max_suggestions = max_suggestions,
        max_longshot_suggestions = max_longshot_suggestions,
        longshot_threshold = longshot_threshold,
        edit_max_dist = edit_max_dist,
        edit_max_candidates = edit_max_candidates,
        timeout_sec = timeout_sec
      )
      user_choice <- prompt_fuzzy_suggestions(suggestions, original, cleaned)

      if (user_choice$action == "cancel") {
        cat("Cancelled. Stopping batch processing.\n")
        return(df)
      } else if (user_choice$action == "skip") {
        next
      } else if (user_choice$action == "flag_for_removal") {
        df$flag_for_removal[idx] <- TRUE
        cat(" Flagged for removal.\n")
        save_progress(df, checkpoint_file)
        next
      } else if (user_choice$action == "expert_review") {
        df$resolution_method[idx] <- "expert_review"
        df$resolution_notes[idx] <- "To be reviewed by an expert"
        cat(" Marked for expert review.\n")
        save_progress(df, checkpoint_file, genus_vocab, epithet_vocab)
      } else if (user_choice$action == "accept") {
        post_retry_action <- NULL

        repeat {
          worms_data <- fetch_worms_full_record(
            user_choice$aphiaid,
            timeout_sec
          )
          if (is.null(worms_data)) {
            cat(" Could not fetch WoRMS record.\n")
            break
          }

          df <- apply_resolution(df, idx, worms_data, user_choice$type)
          cat(sprintf(
            " Assigned: %s (AphiaID: %d)\n",
            worms_data$accepted_name,
            user_choice$aphiaid
          ))

          # Auto-note synonymy
          if (!is.na(worms_data$valid_aphiaid)) {
            syn_note <- sprintf(
              "%s -> %s",
              worms_data$matched_name,
              worms_data$accepted_name
            )
            df$resolution_notes[idx] <- if (!is.na(df$resolution_notes[idx])) {
              paste(df$resolution_notes[idx], syn_note, sep = "; ")
            } else {
              syn_note
            }
            cat(sprintf(" Synonymy noted: %s\n", syn_note))
          }

          # LTR genus assignment prompt
          if (user_choice$type == "genus_fallback") {
            epithet_words <- strsplit(cleaned, " ")[[1]]
            unresolved_tokens <- if (length(epithet_words) >= 2) {
              paste(epithet_words[-1], collapse = " ")
            } else {
              NA_character_
            }
            if (!is.na(unresolved_tokens)) {
              cat(sprintf(" Unresolved tokens: '%s'\n", unresolved_tokens))
            }
            cat(
              " Mark as [l]ower taxonomic resolution or [g]enus-level? [l]/g: "
            )
            ltr_choice <- tolower(trimws(readline()))
            if (ltr_choice != "g") {
              lts_note <- if (!is.na(unresolved_tokens)) {
                sprintf("LTR | unresolved: '%s'", unresolved_tokens)
              } else {
                "LTR"
              }
              df$resolution_notes[idx] <- if (
                !is.na(df$resolution_notes[idx])
              ) {
                paste(df$resolution_notes[idx], lts_note, sep = "; ")
              } else {
                lts_note
              }
              cat(sprintf(" LTR note: %s\n", lts_note))
            } else {
              df$resolution_method[idx] <- "manual_worms_genus"
              inf_note <- if (!is.na(unresolved_tokens)) {
                sprintf("genus-level | informal: '%s'", unresolved_tokens)
              } else {
                "genus-level"
              }
              df$resolution_notes[idx] <- if (
                !is.na(df$resolution_notes[idx])
              ) {
                paste(df$resolution_notes[idx], inf_note, sep = "; ")
              } else {
                inf_note
              }
              cat(" Genus-level assignment.\n")
            }
          }

          # Auto-note binomial suffix
          if (
            !is.null(user_choice$binomial_suffix) &&
              !is.na(user_choice$binomial_suffix)
          ) {
            sfx_note <- sprintf(
              "infraspecific epithet: '%s'",
              user_choice$binomial_suffix
            )
            df$resolution_notes[idx] <- if (!is.na(df$resolution_notes[idx])) {
              paste(df$resolution_notes[idx], sfx_note, sep = "; ")
            } else {
              sfx_note
            }
          }

          notes_result <- prompt_resolution_notes(
            worms_data$matched_name,
            worms_data$accepted_name
          )

          if (notes_result$action == "retry") {
            cat(" Retry assignation\n")
            df$matched_aphiaid[idx] <- NA_integer_
            df$matched_name[idx] <- NA_character_
            df$accepted_name[idx] <- NA_character_
            df$taxonomic_status[idx] <- NA_character_
            df$resolution_method[idx] <- NA_character_
            df$resolution_notes[idx] <- NA_character_
            df$rank[idx] <- NA_character_
            user_choice <- prompt_fuzzy_suggestions(
              suggestions,
              original,
              cleaned
            )
            if (user_choice$action != "accept") {
              post_retry_action <- user_choice$action
              break
            }
            next
          }

          if (!is.na(notes_result$notes)) {
            df$resolution_notes[idx] <- if (!is.na(df$resolution_notes[idx])) {
              paste(df$resolution_notes[idx], notes_result$notes, sep = "; ")
            } else {
              notes_result$notes
            }
            cat(" Notes added.\n")
          }

          save_progress(df, checkpoint_file)
          break
        }

        # Propagate cancel/skip/etc. from retry
        if (!is.null(post_retry_action)) {
          if (post_retry_action == "cancel") {
            cat("Cancelled. Stopping batch processing.\n")
            return(df)
          } else if (post_retry_action == "skip") {
            next
          } else if (post_retry_action == "flag_for_removal") {
            df$flag_for_removal[idx] <- TRUE
            cat(" Flagged for removal.\n")
            save_progress(df, checkpoint_file)
          } else if (post_retry_action == "expert_review") {
            df$resolution_method[idx] <- "expert_review"
            df$resolution_notes[idx] <- "To be reviewed by an expert"
            cat(" Marked for expert review.\n")
            save_progress(df, checkpoint_file)
          }
        }
      }
    }

    if (batch_end < total_unresolved) {
      cat(sprintf(
        "\n%s\nBatch complete. %d entries remaining.\nContinue? [y/n]: ",
        strrep("-", 70),
        total_unresolved - batch_end
      ))
      if (tolower(trimws(readline())) != "y") {
        cat("Stopping.\n")
        break
      }
      batch_start <- batch_end + 1
    } else {
      break
    }
  }

  df
}

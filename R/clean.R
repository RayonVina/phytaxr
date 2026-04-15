#' Normalize special characters in taxon names
#'
#' Converts taxon names to UTF-8, replaces typographic variants (curly quotes,
#' em-dashes, invisible spaces, etc.), collapses whitespace, and initializes
#' the `taxon_clean`, `tax_epithet`, and `uncertain` output columns.
#'
#' @param df A data frame with a `taxon` column (character).
#' @return The input data frame with columns `taxon_clean`, `tax_epithet`
#'   (NA), and `uncertain` (FALSE).
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace_all str_squish
#' @importFrom stringi stri_trans_nfd stri_trans_nfc
#' @export
normalize_characters <- function(df) {
  df |>
    dplyr::mutate(
      taxon_clean = taxon,
      taxon_clean = iconv(taxon_clean, from = "", to = "UTF-8", sub = "byte"),
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        c(
          "\"" = "",
          "'" = "",
          "\u201C" = "",
          "\u201D" = "",
          "\u2018" = "",
          "\u2019" = "",
          "\u00AB" = "",
          "\u00BB" = "",
          "`" = "",
          "\u2013" = "-",
          "\u2014" = "-",
          "\u2212" = "-",
          "\uFF0F" = "/",
          "\u2044" = "/",
          "\\\\" = "/",
          "\u00A0" = " ",
          "\u2009" = " ",
          "\u2008" = " ",
          "\u202F" = " ",
          "\u201A" = ",",
          "\u201E" = ",",
          "\u2264" = "<=",
          "\u2265" = ">=",
          "\u00C2\u00B5" = "\u00B5",
          "\u039C" = "\u00B5",
          "\u03BC" = "\u00B5",
          "_" = " "
        )
      ),
      taxon_clean = stringr::str_replace_all(taxon_clean, "[\\t\\r\\n]+", " "),
      taxon_clean = stringr::str_replace_all(taxon_clean, "\\s*/\\s*", " / "),
      taxon_clean = stringr::str_squish(taxon_clean),
      taxon_clean = stringi::stri_trans_nfc(stringi::stri_trans_nfd(
        taxon_clean
      )),
      tax_epithet = NA_character_,
      uncertain = FALSE
    )
}

#' Extract rank prefixes (O., C., F., P.)
#'
#' Detects and removes rank abbreviation prefixes such as `O.`, `C.`, `F.`,
#' `P.` from `taxon_clean`, storing rank information in `tax_epithet`.
#'
#' @param df A data frame processed by [normalize_characters()].
#' @return The input data frame with updated `taxon_clean`, `tax_epithet`,
#'   and `uncertain`.
#' @importFrom dplyr mutate select case_when if_else
#' @importFrom stringr str_detect str_extract str_remove str_trim str_squish
#' @export
process_taxonomic_prefixes <- function(df) {
  rank_pattern <- "^(O|C|F|P)\\.\\s+"
  df |>
    dplyr::mutate(
      has_rank_prefix = stringr::str_detect(
        taxon_clean,
        stringr::regex(rank_pattern, ignore_case = FALSE)
      ),
      rank_abbrev = dplyr::if_else(
        has_rank_prefix,
        stringr::str_extract(
          taxon_clean,
          stringr::regex("^(O|C|F|P)", ignore_case = FALSE)
        ),
        NA_character_
      ),
      rank_name = dplyr::case_when(
        rank_abbrev == "O" ~ "Order",
        rank_abbrev == "C" ~ "Class",
        rank_abbrev == "F" ~ "Family",
        rank_abbrev == "P" ~ "Phylum",
        TRUE ~ NA_character_
      ),
      prefix_info = dplyr::if_else(
        has_rank_prefix & !is.na(rank_name),
        paste(rank_name, "uncertain"),
        NA_character_
      ),
      tax_epithet = dplyr::case_when(
        has_rank_prefix & !is.na(prefix_info) & !is.na(tax_epithet) ~
          stringr::str_squish(paste(tax_epithet, prefix_info)),
        has_rank_prefix & !is.na(prefix_info) ~ prefix_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = dplyr::if_else(
        has_rank_prefix,
        stringr::str_trim(stringr::str_remove(
          taxon_clean,
          stringr::regex(rank_pattern, ignore_case = FALSE)
        )),
        taxon_clean
      ),
      uncertain = dplyr::if_else(has_rank_prefix, TRUE, uncertain)
    ) |>
    dplyr::select(-has_rank_prefix, -rank_abbrev, -rank_name, -prefix_info)
}

#' Extract incertae sedis and qualification markers
#'
#' Detects and removes qualifiers such as `cf.`, `aff.`, `sensu lato`,
#' `type`, `complex`, `group`, `?`, and parenthetical content from
#' `taxon_clean`, storing them in `tax_epithet`.
#'
#' @param df A data frame processed by [process_taxonomic_prefixes()].
#' @return The input data frame with updated `taxon_clean`, `tax_epithet`,
#'   and `uncertain`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_detect str_extract str_extract_all str_remove
#'   str_remove_all str_replace_all str_squish regex
#' @importFrom purrr map_chr
#' @export
process_incertae_entries <- function(df) {
  qualif_pattern <- "\\b(?:cf|aff|nr|agg|indet|incerta|incertae sedis|sensu\\s+lato|sensu\\s+stricto|type|complex|group|juv|et)\\.?(?=\\b)"
  df |>
    dplyr::mutate(
      uncertain = stringr::str_detect(taxon_clean, qualif_pattern) |
        stringr::str_detect(taxon_clean, "\\?"),
      paren_info = stringr::str_extract_all(taxon_clean, "\\([^)]*\\)") |>
        purrr::map_chr(~ stringr::str_squish(paste(.x, collapse = " "))) |>
        stringr::str_remove_all("[()]"),
      type_info = stringr::str_extract(
        taxon_clean,
        stringr::regex("\\btype\\b.*$", ignore_case = TRUE)
      ),
      complex_info = stringr::str_extract(
        taxon_clean,
        stringr::regex("\\bcomplex\\b.*$", ignore_case = TRUE)
      ),
      group_info = stringr::str_extract(
        taxon_clean,
        stringr::regex("\\bgroup\\b.*$", ignore_case = TRUE)
      ),
      tax_epithet = dplyr::case_when(
        paren_info != "" & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          paren_info
        )),
        paren_info != "" ~ paren_info,
        !is.na(type_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          type_info
        )),
        !is.na(type_info) ~ type_info,
        !is.na(complex_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          complex_info
        )),
        !is.na(complex_info) ~ complex_info,
        !is.na(group_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          group_info
        )),
        !is.na(group_info) ~ group_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = taxon_clean |>
        stringr::str_remove_all(stringr::regex(
          qualif_pattern,
          ignore_case = TRUE
        )) |>
        stringr::str_remove_all(stringr::regex(
          "incertae\\s+sedis",
          ignore_case = TRUE
        )) |>
        stringr::str_remove_all(stringr::regex(
          "sensu\\s+lato|sensu\\s+stricto",
          ignore_case = TRUE
        )) |>
        stringr::str_remove_all(stringr::regex(
          "\\btype\\b.*$",
          ignore_case = TRUE
        )) |>
        stringr::str_remove_all(stringr::regex(
          "\\bcomplex\\b.*$",
          ignore_case = TRUE
        )) |>
        stringr::str_remove_all(stringr::regex(
          "\\bgroup\\b.*$",
          ignore_case = TRUE
        )) |>
        stringr::str_replace_all("\\?", "") |>
        stringr::str_remove_all("\\s*\\([^)]*\\)") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-paren_info, -type_info, -complex_info, -group_info)
}

#' Extract sp./spp. designations
#'
#' Moves `sp.`, `spp.`, `sp1`, etc. from `taxon_clean` to `tax_epithet`.
#'
#' @param df A data frame processed by [process_incertae_entries()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract str_remove str_squish regex
#' @export
process_sp_entries <- function(df) {
  sp_pattern <- stringr::regex(
    "\\b(?:sp\\d+\\.?|sp\\.?|spp\\.?)\\b",
    ignore_case = TRUE
  )
  df |>
    dplyr::mutate(
      sp_num = stringr::str_extract(taxon_clean, sp_pattern),
      tax_epithet = dplyr::case_when(
        !is.na(sp_num) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          sp_num
        )),
        !is.na(sp_num) ~ sp_num,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove(taxon_clean, sp_pattern) |>
        stringr::str_squish()
    ) |>
    dplyr::select(-sp_num)
}

#' Extract bracket annotations (`-[...]`)
#'
#' Moves bracket-style annotations of the form `` -[...] `` from `taxon_clean`
#' to `tax_epithet`.
#'
#' @param df A data frame processed by [process_sp_entries()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract str_remove str_squish
#' @export
process_bracket_entries <- function(df) {
  df |>
    dplyr::mutate(
      bracket_info = stringr::str_extract(taxon_clean, "\\s*-\\[.*"),
      tax_epithet = dplyr::case_when(
        !is.na(bracket_info) & is.na(tax_epithet) ~ stringr::str_squish(
          bracket_info
        ),
        !is.na(bracket_info) ~ stringr::str_squish(paste(
          tax_epithet,
          bracket_info
        )),
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove(taxon_clean, "\\s*-\\[.*") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-bracket_info)
}

#' Extract size information
#'
#' Detects and moves numeric size expressions (e.g. `10-20 \u00B5m`, `>5 mm`)
#' from `taxon_clean` to `tax_epithet`.
#'
#' @param df A data frame processed by [process_bracket_entries()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract_all str_remove_all str_replace_all str_squish regex
#' @importFrom purrr map_chr
#' @export
move_size_to_epithet <- function(df) {
  size_pattern <- paste0(
    "(?:<=|>=|<|>|\u2248|~|\u2264|\u2265)?\\s*",
    "\\d+(?:\\.\\d+)?",
    "(?:\\s*(?:-|\u2013|\u2014)\\s*\\d+(?:\\.\\d+)?)?",
    "\\s*(?:\u00B5m|um|mm)?",
    "(?:\\s*(?:cell\\s*)?(?:width|length|diameter|height))?"
  )
  df |>
    dplyr::mutate(
      size_info = stringr::str_extract_all(
        taxon_clean,
        stringr::regex(size_pattern, ignore_case = TRUE, multiline = TRUE)
      ) |>
        purrr::map_chr(~ stringr::str_squish(paste(.x, collapse = " "))),
      tax_epithet = dplyr::case_when(
        size_info != "" & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          size_info
        )),
        size_info != "" ~ size_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove_all(
        taxon_clean,
        stringr::regex(size_pattern, ignore_case = TRUE, multiline = TRUE)
      ) |>
        stringr::str_replace_all("\u00C2", "") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-size_info)
}

#' Extract remaining parenthetical epithets
#'
#' Moves any remaining content in parentheses from `taxon_clean` to
#' `tax_epithet`.
#'
#' @param df A data frame processed by [move_size_to_epithet()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract str_remove_all str_squish
#' @export
process_epithet_entries <- function(df) {
  df |>
    dplyr::mutate(
      tax_epithet = dplyr::coalesce(
        tax_epithet,
        stringr::str_extract(taxon_clean, "\\(([^)]+)\\)") |>
          stringr::str_remove_all("[()]")
      ),
      taxon_clean = stringr::str_remove_all(taxon_clean, "\\s*\\([^)]*\\)") |>
        stringr::str_squish()
    )
}

#' Normalize infraspecific rank notation
#'
#' Ensures consistent formatting of `var.`, `f.`, `subsp.`, `ssp.`, `cv.`
#' in `taxon_clean` (space before, dot after, space after dot).
#'
#' @param df A data frame processed by [process_epithet_entries()].
#' @return The input data frame with normalized `taxon_clean`.
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace_all str_squish regex
#' @export
normalize_infraspecific_ranks <- function(df) {
  df |>
    dplyr::mutate(
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        stringr::regex("\\b(var|f|subsp|ssp|cv)(?=\\.)", ignore_case = FALSE),
        " \\1"
      ) |>
        stringr::str_squish(),
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        stringr::regex(
          "\\b(var|f|subsp|ssp|cv)(?!\\.)\\s+",
          ignore_case = FALSE
        ),
        "\\1. "
      ),
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        stringr::regex(
          "\\b(var|f|subsp|ssp|cv)\\.(?!\\s)",
          ignore_case = FALSE
        ),
        "\\1. "
      ) |>
        stringr::str_squish()
    )
}

#' Remove stray dots from taxon names
#'
#' Removes dots from `taxon_clean` while protecting dots that belong to
#' infraspecific rank abbreviations (`var.`, `f.`, `subsp.`, etc.).
#'
#' @param df A data frame processed by [normalize_infraspecific_ranks()].
#' @return The input data frame with dots removed from `taxon_clean`.
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace_all str_squish regex
#' @export
remove_dots <- function(df) {
  df |>
    dplyr::mutate(
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        stringr::regex("\\b(var|f|subsp|ssp|cv)\\.", ignore_case = FALSE),
        "\\1\u00A7\u00A7\u00A7"
      ),
      taxon_clean = stringr::str_replace_all(taxon_clean, "\\.", ""),
      taxon_clean = stringr::str_replace_all(
        taxon_clean,
        "\u00A7\u00A7\u00A7",
        "."
      ) |>
        stringr::str_squish()
    )
}

#' Move reproductive and morphological structure terms to tax_epithet
#'
#' Detects terms like `cysts`, `spores`, `colonial`, `benthic`, etc. and
#' moves them from `taxon_clean` to `tax_epithet`.
#'
#' @param df A data frame processed by [remove_dots()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract_all str_remove_all str_replace_all str_squish regex
#' @importFrom purrr map_chr
#' @export
move_reproductive_structures <- function(df) {
  repro_pattern <- "\\b(cysts?|gametes?|resting|rest|spores?|filaments?|cells?|sphere|naked|fusiform|round|small|large|tiny|micro|nano|symbiont|with\\s+symbiont|epiphytic|solitary|colonial|chains?|clustered|green|blue|olive|brown|golden|planktonic|benthic|cluster|coccoid)\\b"
  df |>
    dplyr::mutate(
      reprod_info = stringr::str_extract_all(
        taxon_clean,
        stringr::regex(repro_pattern, ignore_case = TRUE)
      ) |>
        purrr::map_chr(~ stringr::str_squish(paste(.x, collapse = " "))),
      tax_epithet = dplyr::case_when(
        reprod_info != "" & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          reprod_info
        )),
        reprod_info != "" ~ reprod_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove_all(
        taxon_clean,
        stringr::regex(repro_pattern, ignore_case = TRUE)
      ) |>
        stringr::str_replace_all("\\s+", " ") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-reprod_info)
}

#' Move uncertainty descriptor terms to tax_epithet
#'
#' Detects terms like `unidentified`, `undetermined`, `unclassified`, etc.
#' and moves them from `taxon_clean` to `tax_epithet`, also setting
#' `uncertain = TRUE`.
#'
#' @param df A data frame processed by [move_reproductive_structures()].
#' @return The input data frame with updated `taxon_clean`, `tax_epithet`,
#'   and `uncertain`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract_all str_remove_all str_replace_all str_squish regex
#' @importFrom purrr map_chr
#' @export
move_uncertainty_descriptors <- function(df) {
  uncertainty_pattern <- "\\b(unid|unidentified|unknown|not\\s+classified|not\\s+classfied|unclassified|unclassifies|uncertain|undetermined|undeterminable|undertermined|unresolved|unclear|dubious|questionable|other\\s+forms|catalogued)\\b"
  df |>
    dplyr::mutate(
      uncertain_info = stringr::str_extract_all(
        taxon_clean,
        stringr::regex(uncertainty_pattern, ignore_case = TRUE)
      ) |>
        purrr::map_chr(~ stringr::str_squish(paste(.x, collapse = " "))),
      tax_epithet = dplyr::case_when(
        uncertain_info != "" & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          uncertain_info
        )),
        uncertain_info != "" ~ uncertain_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove_all(
        taxon_clean,
        stringr::regex(uncertainty_pattern, ignore_case = TRUE)
      ) |>
        stringr::str_replace_all("\\s+", " ") |>
        stringr::str_squish(),
      uncertain = ifelse(uncertain_info != "", TRUE, uncertain)
    ) |>
    dplyr::select(-uncertain_info)
}

#' Move morphological descriptor terms to tax_epithet
#'
#' Detects terms like `centric`, `pennate`, `thecate`, `armored`, etc. and
#' moves them from `taxon_clean` to `tax_epithet`.
#'
#' @param df A data frame processed by [move_uncertainty_descriptors()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract_all str_remove_all str_replace_all str_squish regex
#' @importFrom purrr map_chr
#' @export
move_morphological_descriptors <- function(df) {
  morpho_pattern <- "\\b(centric|pennate|fusiform|round|spherical|elongated|curved|straight|branched|unbranched|motile|non-motile|thecate|athecate|armored|armoured|unarmored|unarmoured|gymnodinioid)\\b"
  df |>
    dplyr::mutate(
      morpho_info = stringr::str_extract_all(
        taxon_clean,
        stringr::regex(morpho_pattern, ignore_case = TRUE)
      ) |>
        purrr::map_chr(~ stringr::str_squish(paste(.x, collapse = " "))),
      tax_epithet = dplyr::case_when(
        morpho_info != "" & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          morpho_info
        )),
        morpho_info != "" ~ morpho_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove_all(
        taxon_clean,
        stringr::regex(morpho_pattern, ignore_case = TRUE)
      ) |>
        stringr::str_replace_all("\\s+", " ") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-morpho_info)
}

#' Move 'with ...' clauses to tax_epithet
#'
#' Detects `with ...` trailing clauses in `taxon_clean` and moves them to
#' `tax_epithet`.
#'
#' @param df A data frame processed by [move_morphological_descriptors()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract str_remove str_squish regex
#' @export
move_with_descriptors <- function(df) {
  df |>
    dplyr::mutate(
      with_info = stringr::str_extract(
        taxon_clean,
        stringr::regex("\\bwith.*$", ignore_case = TRUE)
      ),
      tax_epithet = dplyr::case_when(
        !is.na(with_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          with_info
        )),
        !is.na(with_info) ~ with_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove(
        taxon_clean,
        stringr::regex("\\bwith.*$", ignore_case = TRUE)
      ) |>
        stringr::str_squish()
    ) |>
    dplyr::select(-with_info)
}

#' Move forma (f.) designations to tax_epithet
#'
#' Detects ` f ` (standalone forma marker without dot) in `taxon_clean`
#' and moves the trailing content to `tax_epithet`.
#'
#' @param df A data frame processed by [move_with_descriptors()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when
#' @importFrom stringr str_extract str_remove str_squish
#' @export
move_formia_to_epithet <- function(df) {
  df |>
    dplyr::mutate(
      forma_info = stringr::str_extract(taxon_clean, "(?<=\\b f \\b)[^\\s].*$"),
      tax_epithet = dplyr::case_when(
        !is.na(forma_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          forma_info
        )),
        !is.na(forma_info) ~ stringr::str_squish(forma_info),
        TRUE ~ tax_epithet
      ),
      taxon_clean = stringr::str_remove(taxon_clean, "\\s+f\\s+(?!\\.).*$") |>
        stringr::str_squish()
    ) |>
    dplyr::select(-forma_info)
}

#' Move comma-separated trailing content to tax_epithet
#'
#' Splits `taxon_clean` at the first comma: keeps the left part in
#' `taxon_clean` and moves the right part to `tax_epithet`.
#'
#' @param df A data frame processed by [move_formia_to_epithet()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_trim str_replace str_extract
#' @export
move_commas_to_epithet <- function(df) {
  df |>
    dplyr::mutate(
      before = stringr::str_trim(stringr::str_replace(taxon_clean, ",.*$", "")),
      after = stringr::str_trim(stringr::str_extract(taxon_clean, "(?<=,).*$")),
      tax_epithet = dplyr::case_when(
        !is.na(after) &
          after != "" &
          !is.na(tax_epithet) ~ stringr::str_squish(paste(tax_epithet, after)),
        !is.na(after) & after != "" ~ after,
        TRUE ~ tax_epithet
      ),
      taxon_clean = before,
      before = NULL,
      after = NULL
    )
}

#' Move author strings to tax_epithet
#'
#' Detects trailing author citations (e.g. `Hustedt 1930`, `(Cleve) Gran`)
#' in `taxon_clean` and moves them to `tax_epithet`.
#'
#' @param df A data frame processed by [move_commas_to_epithet()].
#' @return The input data frame with updated `taxon_clean` and `tax_epithet`.
#' @importFrom dplyr mutate select case_when if_else
#' @importFrom stringr str_detect str_extract str_remove str_squish regex
#' @export
move_authors_to_epithet <- function(df) {
  author_pattern <- stringr::regex(
    paste0(
      "\\s+",
      "(?:",
      "(?:\\([^)]+\\)\\s*)?",
      "(?:[A-Z\\p{Lu}][a-z\\p{Ll}][\\p{L}\\p{M}'-]*",
      "|[A-Z\\p{Lu}][\\p{L}\\p{M}'-]*\\.)",
      "(?:\\s*(?:ex|et|&|and|von|v\\.|de|del|van|der)\\s*",
      "(?:[A-Z\\p{Lu}][a-z\\p{Ll}][\\p{L}\\p{M}'-]*",
      "|[A-Z\\p{Lu}][\\p{L}\\p{M}'-]*\\.?))*",
      "\\s*(?:III|II|Jr\\.)?",
      "\\s*(?:,?\\s*\\d{4})?",
      ")+$"
    ),
    ignore_case = FALSE
  )
  df |>
    dplyr::mutate(
      author_info = dplyr::if_else(
        stringr::str_detect(taxon_clean, "/"),
        NA_character_,
        stringr::str_extract(taxon_clean, author_pattern)
      ),
      author_info = dplyr::if_else(
        !is.na(author_info),
        stringr::str_squish(author_info),
        author_info
      ),
      tax_epithet = dplyr::case_when(
        !is.na(author_info) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          author_info
        )),
        !is.na(author_info) ~ author_info,
        TRUE ~ tax_epithet
      ),
      taxon_clean = dplyr::if_else(
        !is.na(author_info),
        stringr::str_remove(taxon_clean, author_pattern) |>
          stringr::str_squish(),
        taxon_clean
      )
    ) |>
    dplyr::select(-author_info)
}

#' Remove residual sp/spp/ssp tokens from tax_epithet
#'
#' Cleans up redundant `sp`, `spp`, `ssp` tokens that may remain in
#' `tax_epithet` after pipeline processing.
#'
#' @param df A data frame processed by [move_authors_to_epithet()].
#' @return The input data frame with cleaned `tax_epithet`.
#' @importFrom dplyr mutate
#' @importFrom stringr str_remove_all str_squish regex
#' @export
remove_sp_tokens <- function(df) {
  df |>
    dplyr::mutate(
      tax_epithet = stringr::str_remove_all(
        tax_epithet,
        stringr::regex("\\b(?:sp|spp|ssp)\\.?\\b", ignore_case = TRUE)
      ) |>
        stringr::str_squish(),
      tax_epithet = dplyr::na_if(tax_epithet, "")
    )
}

#' Split entries with separators (/, &, +, or) into individual rows
#'
#' Detects taxon names containing `/`, `&`, `+`, or ` or ` and unfolds
#' them into separate rows, recording the aggregate partner in `tax_epithet`.
#'
#' @param df A data frame processed by [remove_sp_tokens()].
#' @return A data frame with potentially more rows than the input, with
#'   updated `taxon_clean`, `tax_epithet`, and `uncertain`.
#' @importFrom dplyr mutate select group_by ungroup distinct case_when if_else n
#' @importFrom stringr str_detect str_replace_all str_trim str_extract
#'   str_squish str_count str_c word
#' @importFrom tidyr uncount
#' @export
split_separator_entries <- function(df) {
  df |>
    dplyr::mutate(original_id = dplyr::row_number()) |>
    dplyr::mutate(
      has_slash = dplyr::if_else(
        !is.na(taxon_clean) &
          stringr::str_detect(taxon_clean, "/") &
          !stringr::str_detect(taxon_clean, "-\\["),
        TRUE,
        FALSE
      ),
      has_ampersand = dplyr::if_else(
        !is.na(taxon_clean) & stringr::str_detect(taxon_clean, "&"),
        TRUE,
        FALSE
      ),
      has_plus = dplyr::if_else(
        !is.na(taxon_clean) & stringr::str_detect(taxon_clean, "\\+"),
        TRUE,
        FALSE
      ),
      has_or = dplyr::if_else(
        !is.na(taxon_clean) & stringr::str_detect(taxon_clean, "\\bor\\b"),
        TRUE,
        FALSE
      ),
      has_separator = has_slash | has_ampersand | has_plus | has_or,
      weight = dplyr::if_else(has_separator, 2L, 1L)
    ) |>
    dplyr::mutate(
      taxon_clean = dplyr::case_when(
        has_slash ~ stringr::str_replace_all(taxon_clean, "\\s*/\\s*", " / "),
        has_ampersand ~ stringr::str_replace_all(
          taxon_clean,
          "\\s*&\\s*",
          " & "
        ),
        has_plus ~ stringr::str_replace_all(taxon_clean, "\\s*\\+\\s*", " + "),
        has_or ~ stringr::str_replace_all(taxon_clean, "\\s+or\\s+", " or "),
        TRUE ~ taxon_clean
      )
    ) |>
    tidyr::uncount(weights = weight, .remove = FALSE) |>
    dplyr::group_by(original_id) |>
    dplyr::mutate(copy_id = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      left_part = dplyr::case_when(
        has_slash ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "^[^/]+"
        )),
        has_ampersand ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "^[^&]+"
        )),
        has_plus ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "^[^+]+"
        )),
        has_or ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "^.*?(?=\\s+or\\s+)"
        )),
        TRUE ~ NA_character_
      ),
      right_part = dplyr::case_when(
        has_slash ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "(?<=/).*"
        )),
        has_ampersand ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "(?<=&).*"
        )),
        has_plus ~ stringr::str_trim(stringr::str_extract(
          taxon_clean,
          "(?<=\\+).*"
        )),
        has_or ~ stringr::str_trim(stringr::str_replace(
          taxon_clean,
          "^.*?\\s+or\\s+",
          ""
        )),
        TRUE ~ NA_character_
      ),
      n_left_items = dplyr::if_else(
        has_separator,
        stringr::str_count(left_part, "\\S+"),
        NA_integer_
      )
    ) |>
    dplyr::mutate(
      taxon_clean_final = dplyr::case_when(
        !has_separator ~ taxon_clean,
        n_left_items == 1 & copy_id == 1L ~ left_part,
        n_left_items == 1 & copy_id == 2L ~ right_part,
        n_left_items > 1 & copy_id == 1L ~ left_part,
        n_left_items > 1 & copy_id == 2L ~ dplyr::if_else(
          stringr::str_detect(right_part, "\\s"),
          right_part,
          stringr::str_c(stringr::word(left_part, 1), " ", right_part)
        ),
        TRUE ~ taxon_clean
      ),
      agg_partner = dplyr::case_when(
        !has_separator ~ NA_character_,
        copy_id == 1L & n_left_items == 1 ~ paste("Agg:", right_part),
        copy_id == 2L & n_left_items == 1 ~ paste("Agg:", left_part),
        copy_id == 1L & n_left_items > 1 ~ paste(
          "Agg:",
          dplyr::if_else(
            stringr::str_detect(right_part, "\\s"),
            right_part,
            stringr::str_c(stringr::word(left_part, 1), " ", right_part)
          )
        ),
        copy_id == 2L & n_left_items > 1 ~ paste("Agg:", left_part),
        TRUE ~ NA_character_
      ),
      tax_epithet = dplyr::case_when(
        !is.na(agg_partner) & !is.na(tax_epithet) ~ stringr::str_squish(paste(
          tax_epithet,
          agg_partner
        )),
        !is.na(agg_partner) ~ agg_partner,
        TRUE ~ tax_epithet
      ),
      taxon_clean = taxon_clean_final,
      uncertain = dplyr::if_else(has_separator, TRUE, uncertain)
    ) |>
    dplyr::select(
      -original_id,
      -left_part,
      -right_part,
      -n_left_items,
      -copy_id,
      -has_separator,
      -has_slash,
      -has_ampersand,
      -has_plus,
      -has_or,
      -weight,
      -taxon_clean_final,
      -agg_partner
    ) |>
    dplyr::distinct()
}

#' Map common English names to accepted taxonomic names
#'
#' Replaces informal names and morphotypes (e.g. `diatom`, `ciliate`,
#' `gymnodinioid`) with the corresponding accepted higher taxon name.
#'
#' @param df A data frame processed by [split_separator_entries()].
#' @return The input data frame with updated `taxon_clean`, `tax_epithet`,
#'   and `uncertain`.
#' @importFrom dplyr mutate select case_when if_else
#' @importFrom stringr str_detect str_extract str_remove_all str_squish
#'   str_to_lower regex
#' @export
process_generic_taxa <- function(df) {
  common_names_map <- c(
    "dinoflagellate" = "Dinoflagellata",
    "dinoflagellates" = "Dinoflagellata",
    "gymnodinioid" = "Gymnodiniales",
    "gymnodinioids" = "Gymnodiniales",
    "diatom" = "Bacillariophyceae",
    "diatoms" = "Bacillariophyceae",
    "centric" = "Bacillariophyceae",
    "centrics" = "Bacillariophyceae",
    "pennate" = "Bacillariophyceae",
    "pennates" = "Bacillariophyceae",
    "penate" = "Bacillariophyceae",
    "penates" = "Bacillariophyceae",
    "Ciliater" = "Ciliophora",
    "ciliate" = "Ciliophora",
    "ciliates" = "Ciliophora",
    "ciliater" = "Ciliophora"
  )
  key_terms <- names(common_names_map)
  valid_genera_exceptions <- c("Diatoma", "Diatomaceae", "Gymnodinium")
  uncertainty_prefixes <- c(
    "unid",
    "unidentified",
    "unclassified",
    "undetermined",
    "unarmored",
    "unarmoured",
    "other",
    "naked",
    "thecate",
    "thecae",
    "armoured",
    "armored",
    "empty",
    "benthic",
    "planktonic",
    "not classified",
    "classified",
    "inc",
    "tintinnid",
    "tintinnids"
  )
  size_descriptors <- c(
    "micro",
    "nano",
    "pico",
    "small",
    "medium",
    "large",
    "tiny",
    "big"
  )
  morpho_descriptors <- c(
    "fusiform",
    "round",
    "spherical",
    "elongated",
    "curved",
    "straight"
  )

  df |>
    dplyr::mutate(
      is_valid_genus = stringr::str_detect(
        taxon_clean,
        stringr::regex(
          paste0("^(", paste(valid_genera_exceptions, collapse = "|"), ")\\b"),
          ignore_case = TRUE
        )
      ),
      has_generic_term = stringr::str_detect(
        taxon_clean,
        stringr::regex(
          paste0("\\b(", paste(key_terms, collapse = "|"), ")\\b"),
          ignore_case = TRUE
        )
      ),
      should_process = has_generic_term & !is_valid_genus & !is.na(taxon_clean),
      main_term = dplyr::if_else(
        should_process,
        stringr::str_extract(
          taxon_clean,
          stringr::regex(
            paste0("\\b(", paste(key_terms, collapse = "|"), ")\\b"),
            ignore_case = TRUE
          )
        ),
        NA_character_
      ),
      correct_taxon = dplyr::if_else(
        should_process & !is.na(main_term),
        common_names_map[stringr::str_to_lower(main_term)],
        NA_character_
      ),
      should_be_uncertain = should_process &
        (stringr::str_detect(
          taxon_clean,
          stringr::regex(
            paste0("\\b(", paste(uncertainty_prefixes, collapse = "|"), ")\\b"),
            ignore_case = TRUE
          )
        ) |
          stringr::str_detect(taxon_clean, "\\?") |
          stringr::str_detect(taxon_clean, "\\bsp\\b") |
          stringr::str_detect(
            taxon_clean,
            stringr::regex(
              paste0("\\b(", paste(size_descriptors, collapse = "|"), ")\\b"),
              ignore_case = TRUE
            )
          ) |
          stringr::str_detect(
            taxon_clean,
            stringr::regex(
              paste0("\\b(", paste(morpho_descriptors, collapse = "|"), ")\\b"),
              ignore_case = TRUE
            )
          ) |
          stringr::str_detect(taxon_clean, "\\d+\\s*[-\u2013]\\s*\\d+") |
          stringr::str_detect(taxon_clean, "\\d+\\s*[\u00B5u]m") |
          stringr::str_detect(taxon_clean, "[<>\\u2264\\u2265]\\s*\\d+") |
          stringr::str_detect(taxon_clean, "\\([^)]+\\)") |
          stringr::str_detect(taxon_clean, "\\s+-\\s+")),
      taxon_clean = dplyr::if_else(
        should_process & !is.na(correct_taxon),
        correct_taxon,
        taxon_clean
      ),
      uncertain = dplyr::if_else(should_be_uncertain, TRUE, uncertain)
    ) |>
    dplyr::select(
      -is_valid_genus,
      -has_generic_term,
      -should_process,
      -main_term,
      -correct_taxon,
      -should_be_uncertain
    )
}

#' Remove trailing hyphens and isolated dashes
#'
#' Final cleanup: removes trailing ` - ` constructs and isolated hyphens
#' left over after other cleaning steps.
#'
#' @param df A data frame processed by [process_generic_taxa()].
#' @return The input data frame with cleaned `taxon_clean`.
#' @importFrom dplyr mutate
#' @importFrom stringr str_remove str_replace_all str_squish
#' @export
clean_trailing_hyphens <- function(df) {
  df |>
    dplyr::mutate(
      taxon_clean = stringr::str_remove(taxon_clean, "\\s*-\\s*$") |>
        stringr::str_squish(),
      taxon_clean = stringr::str_replace_all(taxon_clean, "\\s+-\\s+", " ") |>
        stringr::str_squish()
    )
}

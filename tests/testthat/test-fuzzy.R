# build_genus_vocabulary() -----------------------------------------------

test_that("build_genus_vocabulary() returns sorted unique genera", {
  df <- data.frame(
    genus = c("Thalassiosira", "Chaetoceros", "Thalassiosira", NA),
    stringsAsFactors = FALSE
  )
  result <- build_genus_vocabulary(df)
  expect_type(result, "character")
  expect_equal(result, c("Chaetoceros", "Thalassiosira"))
})

test_that("build_genus_vocabulary() excludes short names (< 4 chars)", {
  df <- data.frame(
    genus = c("Tha", "Chaetoceros", "Ab"),
    stringsAsFactors = FALSE
  )
  result <- build_genus_vocabulary(df)
  expect_equal(result, "Chaetoceros")
})

test_that("build_genus_vocabulary() returns empty vector when all NA", {
  df <- data.frame(genus = NA_character_, stringsAsFactors = FALSE)
  result <- build_genus_vocabulary(df)
  expect_length(result, 0)
})

# build_epithet_vocabulary() ---------------------------------------------

test_that("build_epithet_vocabulary() extracts epithets correctly", {
  df <- data.frame(
    genus   = c("Thalassiosira", "Chaetoceros"),
    species = c("Thalassiosira weissflogii", "Chaetoceros debilis"),
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_type(result, "character")
  expect_true("weissflogii" %in% result)
  expect_true("debilis" %in% result)
})

test_that("build_epithet_vocabulary() excludes short epithets (< 4 chars)", {
  df <- data.frame(
    genus   = "Genus",
    species = "Genus sp",
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_length(result, 0)
})

test_that("build_epithet_vocabulary() returns empty vector when species is NA", {
  df <- data.frame(
    genus   = "Thalassiosira",
    species = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_length(result, 0)
})

# normalize_taxonomic_name() ---------------------------------------------

test_that("normalize_taxonomic_name() returns NA for NA input", {
  result <- phytaxr:::normalize_taxonomic_name(NA_character_, character(0))
  expect_true(is.na(result))
})

test_that("normalize_taxonomic_name() splits concatenated genus+epithet", {
  vocab <- c("Chaetoceros")
  result <- phytaxr:::normalize_taxonomic_name("Chaetocerosdebilis", vocab)
  expect_equal(result, "Chaetoceros debilis")
})

test_that("normalize_taxonomic_name() returns original if no rule applies", {
  vocab <- c("Thalassiosira")
  result <- phytaxr:::normalize_taxonomic_name("Skeletonema costatum", vocab)
  expect_equal(result, "Skeletonema costatum")
})

# apply_resolution() -----------------------------------------------------
# Tests for apply_resolution() and save_progress() removed: they depended on
# initialize_taxonomy_df() which was removed from the package.

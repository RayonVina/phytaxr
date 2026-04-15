test_that("calculate_similarity() returns 1 for identical strings", {
  expect_equal(phytaxr:::calculate_similarity("Skeletonema", "Skeletonema"), 1)
})

test_that("calculate_similarity() returns 0 for NA inputs", {
  expect_equal(phytaxr:::calculate_similarity(NA, "Skeletonema"), 0)
  expect_equal(phytaxr:::calculate_similarity("Skeletonema", NA), 0)
})

test_that("calculate_similarity() is case-insensitive", {
  expect_equal(
    phytaxr:::calculate_similarity("skeletonema", "Skeletonema"),
    1
  )
})

test_that("calculate_similarity() returns partial score for similar strings", {
  score <- phytaxr:::calculate_similarity("Skeletonema", "Skeletoneme")
  expect_gt(score, 0.8)
  expect_lt(score, 1)
})

# ---------------------------------------------------------------------------
# verify_genus_exists() -- mocked via .wm_records_name_safe
# ---------------------------------------------------------------------------

test_that("verify_genus_exists() returns FALSE for NA or empty input", {
  expect_false(phytaxr:::verify_genus_exists(NA_character_))
  expect_false(phytaxr:::verify_genus_exists(""))
})

test_that("verify_genus_exists() returns TRUE when WoRMS finds a genus", {
  fake_result <- data.frame(
    rank            = "Genus",
    scientificname  = "Skeletonema",
    AphiaID         = 149151L,
    stringsAsFactors = FALSE
  )
  with_mocked_bindings(
    .wm_records_name_safe = function(genus_name, timeout_sec = 10) fake_result,
    .package = "phytaxr",
    expect_true(phytaxr:::verify_genus_exists("Skeletonema"))
  )
})

test_that("verify_genus_exists() returns FALSE when WoRMS finds no genus rank", {
  fake_result <- data.frame(
    rank            = "Species",
    scientificname  = "Fakeus invented",
    AphiaID         = 999999L,
    stringsAsFactors = FALSE
  )
  with_mocked_bindings(
    .wm_records_name_safe = function(genus_name, timeout_sec = 10) fake_result,
    .package = "phytaxr",
    expect_false(phytaxr:::verify_genus_exists("Fakeus"))
  )
})

test_that("verify_genus_exists() returns FALSE when WoRMS returns NULL", {
  with_mocked_bindings(
    .wm_records_name_safe = function(genus_name, timeout_sec = 10) NULL,
    .package = "phytaxr",
    expect_false(phytaxr:::verify_genus_exists("Anythingus"))
  )
})

test_that("verify_genus_exists() returns FALSE when WoRMS returns empty df", {
  empty_df <- data.frame(
    rank           = character(0),
    scientificname = character(0),
    AphiaID        = integer(0),
    stringsAsFactors = FALSE
  )
  with_mocked_bindings(
    .wm_records_name_safe = function(genus_name, timeout_sec = 10) empty_df,
    .package = "phytaxr",
    expect_false(phytaxr:::verify_genus_exists("Emptius"))
  )
})

# tests/testthat/test-vernacular.R

test_that("vernacular_list() returns a named character vector", {
  d <- vernacular_list()
  expect_type(d, "character")
  expect_true(length(d) > 0)
  expect_false(is.null(names(d)))
})

test_that("vernacular_lookup() finds an existing entry case-insensitively", {
  d    <- vernacular_list()
  key  <- names(d)[1]
  taxon <- unname(d[1])

  expect_equal(vernacular_lookup(key),          taxon)
  expect_equal(vernacular_lookup(toupper(key)),  taxon)
  expect_equal(vernacular_lookup(paste0(" ", key, " ")), taxon) # trimws
})

test_that("vernacular_lookup() returns NA for unknown name", {
  expect_true(is.na(vernacular_lookup("__no_such_name_xyz__")))
})

test_that("vernacular_add() inserts a new entry", {
  vernacular_add("test diatom", "Testus diatomus")
  expect_equal(vernacular_lookup("test diatom"), "Testus diatomus")
  # cleanup
  vernacular_remove("test diatom")
})

test_that("vernacular_add() refuses to overwrite an existing key", {
  vernacular_add("test diatom2", "Testus diatomus")
  expect_message(
    vernacular_add("test diatom2", "Other taxon"),
    "already exists"
  )
  # value unchanged
  expect_equal(vernacular_lookup("test diatom2"), "Testus diatomus")
  vernacular_remove("test diatom2")
})

test_that("vernacular_update() changes an existing entry", {
  vernacular_add("test alga", "Oldus algus")
  vernacular_update("test alga", "Newus algus")
  expect_equal(vernacular_lookup("test alga"), "Newus algus")
  vernacular_remove("test alga")
})

test_that("vernacular_update() warns when key does not exist", {
  expect_message(
    vernacular_update("__nonexistent__", "Foo bar"),
    "not found"
  )
})

test_that("vernacular_remove() deletes an entry", {
  vernacular_add("temp entry", "Tempus entrius")
  vernacular_remove("temp entry")
  expect_true(is.na(vernacular_lookup("temp entry")))
})

test_that("vernacular_remove() warns when key does not exist", {
  expect_message(
    vernacular_remove("__nonexistent__"),
    "not found"
  )
})

test_that("vernacular_search() matches in keys", {
  vernacular_add("search test alpha", "Alphus testus")
  res <- vernacular_search("search test", where = "keys")
  expect_true("search test alpha" %in% names(res))
  vernacular_remove("search test alpha")
})

test_that("vernacular_search() matches in values", {
  vernacular_add("search val test", "Uniqueus valueus")
  res <- vernacular_search("Uniqueus", where = "values")
  expect_true("search val test" %in% names(res))
  vernacular_remove("search val test")
})

test_that("vernacular_search() matches in both keys and values", {
  vernacular_add("combo key test", "Combous taxon")
  vernacular_add("another entry", "Testus combous")
  res <- vernacular_search("combo", where = "both")
  expect_true("combo key test" %in% names(res))
  expect_true("another entry" %in% names(res))
  vernacular_remove("combo key test")
  vernacular_remove("another entry")
})

test_that("vernacular_search() returns empty vector when no match", {
  res <- vernacular_search("__zzz_no_match_zzz__")
  expect_length(res, 0)
})

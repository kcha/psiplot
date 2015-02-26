context("Test getting sample names from PSI table")

test_that("Get sample names", {
  target <- get_psi_samples(psi)
  expected <- c("Sample1", "Sample2", "Sample3", "Sample4")
  expect_equal(target, expected)
})

test_that("Get sample names as indices", {
  target <- get_psi_samples(psi, value = FALSE)
  expected <- c(7,9,11,13)
  expect_equal(target, expected)
})

context("Test getting sample names with config")

test_that("Get sample names with config", {
  target <- get_psi_samples(psi, configdf = config[1:3,], value = TRUE)
  expected <- config$SampleName[1:3]
  expect_equal(target, expected)
})

test_that("Get sample names with config as indices", {
  target <- get_psi_samples(psi, configdf = config[1:3,], value = FALSE)
  expected <- c(13,11,9)
  expect_equal(target, expected)
})

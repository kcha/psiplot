context("Test getting sample names from PSI table")

test_that("Get sample names", {
  target <- get_psi_samples(psi)
  expected <- c("Sample1",
                "Sample2",
                "Sample3",
                "Sample4",
                "Sample5",
                "Sample6",
                "Sample7",
                "Sample8")
  expect_equal(target, expected)
})

test_that("Get sample names as indices", {
  target <- get_psi_samples(psi, value = FALSE)
  expected <- c(7,9,11,13,15,17,19,21)
  expect_equal(target, expected)
})

test_that("Get samples if quality columns are out of order", {
  psi2 <- psi
  set.seed(123)
  colnames(psi2)[7:ncol(psi2)] <- sample(colnames(psi2)[7:ncol(psi2)])
  target <- get_psi_samples(psi, value = TRUE)
  expected <- c("Sample1",
                "Sample2",
                "Sample3",
                "Sample4",
                "Sample5",
                "Sample6",
                "Sample7",
                "Sample8")
  expect_true(all(target %in% expected))
})

context("Test getting sample names with config")

test_that("Get sample names with config", {
  target <- get_psi_samples(psi, config = config[1:3,], value = TRUE)
  expected <- config$SampleName[1:3]
  expect_equal(target, expected)
})

test_that("Get sample names with config as indices", {
  target <- get_psi_samples(psi, config = config[1:3,], value = FALSE)
  expected <- c(21,7,9)
  expect_equal(target, expected)
})

context("Test re-ordering of event data using configuration file")

config <- data.frame(Order=c(1,2),
                     SampleName=c("Sample1", "Sample4"),
                     GroupName=c("ESC", "Neural"),
                     RColorCode=c("red", "blue")
                     )
formatted_psi <- format_table(psi)
formatted_crpkm <- format_table(crpkm, expr = TRUE)

test_that("Only samples in config are retained", {
  r <- preprocess_sample_colors(formatted_psi, config = config)
  expect_equal(ncol(r$data), nrow(config))
})

context("Test re-ordering of expression data using configuration file")

test_that("Quality scores is NULL", {
  r <- preprocess_sample_colors(formatted_crpkm, config = config, expr = TRUE)
  expect_true(is.null(r$qual))
  expect_equal(ncol(r$data), nrow(config))
  expect_equal(nrow(r$data), nrow(formatted_crpkm))
})

test_that("Quality scores is NULL when config is not used", {
  r <- preprocess_sample_colors(formatted_crpkm, config = NULL, expr = TRUE)
  expect_true(is.null(r$qual))
  expect_equal(ncol(r$data), 8)
  expect_equal(nrow(r$data), nrow(formatted_crpkm))
})

context("Test absence of optional columns")

test_that("Natural order is used with no Order column is specified", {
  config2 <- config
  config2$Order <- NULL
  r <- preprocess_sample_colors(formatted_psi, config2)
  expect_true(all(r$sample_order$SampleOrder== 1:nrow(config2)), "Using natural order")
  expect_true("Order" %in% colnames(r$config))
})

test_that("Default colors are used if RColorCode is missing", {
  config2 <- config
  config2$RColorCode <- NULL
  r <- preprocess_sample_colors(formatted_psi, config2)
  expect_true("RColorCode" %in% colnames(r$config))
})

context("Test that input PSI table is formattted correctly")

test_that("Error is returned if first column is not ID", {
  expect_error(preprocess_sample_colors(psi, config))
})

context("Test re-ordering of event data using configuration file")

config <- data.frame(Order=c(1,2),
                     SampleName=c("Sample1", "Sample4"),
                     GroupName=c("ESC", "Neural"),
                     RColorCode=c("red", "blue")
                     )

test_that("Only samples in config are retained", {
  r <- preprocess_sample_colors(psi, config = config)
  expect_equal(ncol(r$data), nrow(config))
})

context("Test re-ordering of expression data using configuration file")

test_that("Quality scores is NULL", {
  r <- preprocess_sample_colors(crpkm[,c(1,3:ncol(crpkm))], config = config, expr = TRUE)
  expect_true(is.null(r$qual))
  expect_equal(ncol(r$data), nrow(config))
  expect_equal(nrow(r$data), nrow(crpkm))
})

test_that("Quality scores is NULL when config is not used", {
  r <- preprocess_sample_colors(crpkm[,c(1,3:ncol(crpkm))], config = NULL, expr = TRUE)
  expect_true(is.null(r$qual))
  expect_equal(ncol(r$data), 8)
  expect_equal(nrow(r$data), nrow(crpkm))
})

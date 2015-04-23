context("Test plot title creation")

xname <- "S|TEST|chrX:100-101|1"

test_that("Rowname converts to human-friendly title", {
  target <- make_title(xname)
  expected <- "TEST\n(position = chrX:100-101, length = 1, type = S)"
  expect_equal(target, expected)
})

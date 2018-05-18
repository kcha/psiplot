context("Test formatting of PSI input table")

x <- data.frame(GENE = "TEST",
                EVENT = "HsaEX00XXXXXXX",
                COORD = "chrX:100-101",
                LENGTH = 1,
                FullCO = "chrX:90,100-101,110",
                COMPLEX = "S",
                SampleA = 0,
                SampleA.Q = "N,N,N,Bn,S@0,0",
                SampleB = 100,
                SampleB.Q = "SOK,SOK,SOK,OK,S@90,0",
                SampleC = 100,
                SampleC.Q = "N,N,N,Bn,S@1,0",
                stringsAsFactors=FALSE
)

test_that("PSI values with N quality scores are converted to NA", {
  target <- format_table(x)
  expected <- data.frame(ID = "S|TEST|HsaEX00XXXXXXX|chrX:100-101|1",
                         SampleA = as.numeric(NA),
                         SampleA.Q = "N,N,N,Bn,S@0,0",
                         SampleB = 100,
                         SampleB.Q = "SOK,SOK,SOK,OK,S@90,0",
                         SampleC = as.numeric(NA),
                         SampleC.Q = "N,N,N,Bn,S@1,0",
                         stringsAsFactors=FALSE)
  expect_equal(target, expected)
})

test_that("An error is thrown if input is missing the correct first column: GENE", {
  expect_error(format_table(x[,2:ncol(x)]))
  expect_error(format_table(x[,1:7]))
})

context("Test formatting of cRPKM input table")

z <- data.frame(ID = "ENSG000000000001",
                NAME = "TEST",
                SampleA = 0,
                SampleB = 1,
                SampleC = 100,
                stringsAsFactors=FALSE
)

test_that("Row names for cRPKM", {
  target <- format_table(z, expr = TRUE)
  expected <- data.frame(
                  ID = "TEST|ENSG000000000001",
                  SampleA = 0,
                  SampleB = 1,
                  SampleC = 100,
                  stringsAsFactors=FALSE
  )
  expect_equal(target, expected)
})

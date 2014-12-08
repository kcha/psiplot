context("Test conversion of low/bad quality PSI values")

x <- data.frame(SampleA = 0,
                SampleA.Q = "N,N,N,Bn,S@0,0",
                SampleB = 100,
                SampleB.Q = "SOK,SOK,SOK,OK,S@90,0",
                SampleC = 100,
                SampleC.Q = "N,N,N,Bn,S@1,0",
                stringsAsFactors=FALSE
)

test_that("PSI values with N quality scores are converted to NA", {
  target <- convert_psi(x)
  expected <- data.frame(SampleA = as.numeric(NA),
                         SampleA.Q = "N,N,N,Bn,S@0,0",
                         SampleB = 100,
                         SampleB.Q = "SOK,SOK,SOK,OK,S@90,0",
                         SampleC = as.numeric(NA),
                         SampleC.Q = "N,N,N,Bn,S@1,0",
                         stringsAsFactors=FALSE)
  expect_equal(target, expected)
})

x$SampleA.Q = "OK,OK,OK,OK,S@0,0"
x$SampleC.Q = "OK,OK,OK,OK,S@1,0"
test_that("PSI values with good quality scores are not converted to NA", {
  target <- convert_psi(x)
  expect_equal(target, x)
})

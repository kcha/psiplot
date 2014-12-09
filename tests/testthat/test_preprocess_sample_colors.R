context("Test re-ordering of event data using configuration file")

config <- data.frame(Order=c(1,2),
                     SampleName=c("Sample1", "Sample4"),
                     GroupName=c("Muscle", "Neural"),
                     RColorCode=c("red", "blue")
                     )

test_that("Only samples in config are retained", {
    r <- preprocess_sample_colors(psi, config = config)
    expect_equal(ncol(r$data), nrow(config))
    })

context("Test re-ordering of cRPKM data using configuration file")

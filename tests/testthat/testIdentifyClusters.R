context("Identify clusters")

testData = tibble::as.tibble(
  dget("extdata/randomData.txt", keep.source = FALSE)) # Tibble containing filter test data

test_that("A vector with sample names and cluster IDs",
  expect_equal(identifyClusters(testData,20000,
                                positionHeader = "start",
                                chromHeader = "chrom",
                                sampleIdHeader = "sampleIDs"),
               dget("extdata/expectedResultsIdentify.txt",keep.source = F))
)

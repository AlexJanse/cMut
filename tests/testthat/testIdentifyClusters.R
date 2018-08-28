context("Identify clusters")

testData = as.tibble(
  dget("extdata/expectedResultsRandom.txt", keep.source = FALSE)) # Tibble containing filter test data

test_that("A vector with sample names and cluster IDs",
  expect_equal(identifyClusters(testData,200,
                                positionHeader = "start",
                                chromHeader = "chrom",
                                sampleIdHeader = "sampleIDs"),
               dget("extdata/expectedResultsIdentify.txt",keep.source = F))
)

context("Identify and annotate clusters")

testData = tibble::as.tibble(
  dget("extdata/randomData.txt", keep.source = FALSE)) # Tibble containing filter test data

test_that("Check if the results are the same as the original file",
          expect_equal(
            dplyr::all_equal(
              identifyAndAnnotateClusters(testData,20000,
                                          positionHeader = "start",
                                          chromHeader = "chrom",
                                          sampleIdHeader = "sampleIDs"),
              tibble::as.tibble(dget("extdata/expectedResultsAnnotate.txt",keep.source = F))),
            TRUE)
          )

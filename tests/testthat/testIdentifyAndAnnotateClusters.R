context("Identify and annotate clusters")

testData = as.tibble(
  dget("extdata/expectedResultsRandom.txt", keep.source = FALSE)) # Tibble containing filter test data

test_that("Check if the results are the same as the original file",
          expect_equal(
            isTRUE(all_equal(
              identifyAndAnnotateClusters(testData,200,
                                          positionHeader = "start",
                                          chromHeader = "chrom",
                                          sampleIdHeader = "sampleIDs"),
              as.tibble(dget("extdata/expectedResultsAnnotate.txt",keep.source = F)))),
            TRUE)
          )

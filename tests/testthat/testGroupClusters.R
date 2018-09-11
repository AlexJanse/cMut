context("Group clusters and show mutation symbols")

testData = tibble::as.tibble(
  dget("extdata/expectedResultsAnnotate.txt", keep.source = FALSE)) # Tibble containing the data after identifyAndAnnotateClusters function

test_that("Check if the results are the same as the original file",
          expect_equal(
            dplyr::all_equal(
              groupClusters(testData,
                            altHeader = "alt",
                            refHeader = "ref",
                            clusterIdHeader = "clusterId")[,c(1,5,6,7)],
              tibble::as.tibble(dget("extdata/expectedResultsGroupClusters.txt",keep.source = F))[,c(1,5,6,7)]),
            TRUE)
)


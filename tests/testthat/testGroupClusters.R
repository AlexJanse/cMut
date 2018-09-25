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

testGroupClusters <- groupClusters(identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T),patternIntersect = T)

test_that("Check if the patterns match with the expected results",{
          expect_equal(
            testGroupClusters[testGroupClusters$has.intersect == T,][1,8][[1]],
            list(c("AID","AID(hotspot)")),
            TRUE)
          expect_equal(
            testGroupClusters[testGroupClusters$has.intersect == T,][2,8][[1]],
            list(c("A3F")),
            TRUE)
          expect_equal(
            testGroupClusters[testGroupClusters$has.intersect == T,][3,8][[1]],
            list(c("MMR(2)","MMR(1)[Rev.Com.]")),
            TRUE)
         })

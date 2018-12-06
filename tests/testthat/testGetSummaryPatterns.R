context("Check the getSummaryPatterns function")

testGroupClusters <- groupClusters(identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T),
                                   patternIntersect = T,searchClusterPatterns = T)
testSummary <- getSummaryPatterns(testGroupClusters)

test_that("Check if the results match with the expectation",{
  expect_equal(ncol(testSummary),
               3)
  expect_equal(any(testSummary$process == "Unidentified"),
               TRUE)
  expect_equal(testSummary[testSummary$process == "AID","frequency"][[1]],
               28)
  expect_equal(testSummary[testSummary$process == "A3F","frequency"][[1]],
               2)
  expect_equal(testSummary[testSummary$process == "MMR","frequency"][[1]],
               2)
  expect_equal(testSummary[testSummary$process == "PolZeta","frequency"][[1]],
               6)
  expect_equal(testSummary[testSummary$process == "PolZeta.endOnly","frequency"][[1]],
               12)
  expect_equal(nrow(testSummary),
               10)
})


#' convertFactor
#' @description A function to convert factor columns to character
convertFactor <- function(x){
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], as.character)
  return(x)
}

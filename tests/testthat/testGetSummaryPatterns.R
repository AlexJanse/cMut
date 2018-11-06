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
               2)
  expect_equal(testSummary[testSummary$process == "A3F","frequency"][[1]],
               2)
  expect_equal(testSummary[testSummary$process == "MMR","frequency"][[1]],
               2)
  expect_equal(testSummary[testSummary$process == "PolZetaPrimair","frequency"][[1]],
               2)
  expect_equal(testSummary[testSummary$process == "PolZetaSecundair","frequency"][[1]],
               2)
  expect_equal(nrow(testSummary),
               11)
})

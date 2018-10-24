context("Identify and annotate clusters")

testResults <- identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T)
validationResults <- identifyAndAnnotateClusters(validationTable,20000, sampleIdHeader = "id" ,linkPatterns = 2, patternsAsList = F)
test_that("Check if the linked patterns are as expected",{
  expect_equal(
    nrow(testResults[testResults$is.linked == T, ]) == 28,
    TRUE,
    TRUE)
  expect_equal(
    testResults[testResults$sampleIDs == "TEST","linkedPatterns"][1,][[1]][[1]],
    c("AID"),
    TRUE)
  expect_equal(
    nrow(validationResults[grepl("AID",validationResults$linkedPatterns),]),
    156
  )
  expect_equal(
    nrow(validationResults[grepl("MMR",validationResults$linkedPatterns),]),
    8
  )
  expect_equal(
    nrow(validationResults[grepl("A1/A3G",validationResults$linkedPatterns),]),
    44
  )
  expect_equal(
    nrow(validationResults[grepl("A3F",validationResults$linkedPatterns),]),
    12
  )
  expect_equal(
    nrow(validationResults[grepl("A3A",validationResults$linkedPatterns),]),
    8
  )
  expect_equal(
    nrow(validationResults[grepl("A3B",validationResults$linkedPatterns),]),
    8
  )
  expect_equal(
    nrow(validationResults[grepl("A3[^GFAB]",validationResults$linkedPatterns),]),
    32
  )
})


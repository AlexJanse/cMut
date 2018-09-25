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

testResults <- identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T)
test_that("Check if the linked patterns are as expected",{
  expect_equal(
    nrow(testResults[testResults$is.linked == T, ]) == 40,
    TRUE,
    TRUE)
  expect_equal(
    testResults[testResults$sampleIDs == "TEST","linkedPatterns"][1,][[1]][[1]],
    c("AID","AID(hotspot)"),
    TRUE)
  expect_equal(
    testResults[testResults$sampleIDs == "TEST","linkedPatterns"][2,][[1]][[1]],
    c("AID","AID(hotspot)","AID(pref.hotspot)"),
    TRUE)
})

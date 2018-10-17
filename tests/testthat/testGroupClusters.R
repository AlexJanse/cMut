context("Group clusters and show mutation symbols")

testGroupClusters <- groupClusters(identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T),patternIntersect = T)
validationGroups <- groupClusters(identifyAndAnnotateClusters(validationTable,20000, sampleIdHeader = "id" ,linkPatterns = 2, patternsAsList = F),patternIntersect = T)
test_that("Check if the patterns match with the expected results",{
          expect_equal(
            testGroupClusters[grepl("TEST",testGroupClusters$clusterId),][1,9][[1]],
            list(c("AID")),
            TRUE)
          expect_equal(
            nrow(validationGroups[grepl("AID",validationGroups$foundPatterns),]),
            156/2
          )
          expect_equal(
            nrow(validationGroups[grepl("MMR",validationGroups$foundPatterns),]),
            48/2
          )
          expect_equal(
            nrow(validationGroups[grepl("A1/A3G",validationGroups$foundPatterns),]),
            44/2
          )
          expect_equal(
            nrow(validationGroups[grepl("A3F",validationGroups$foundPatterns),]),
            12/2
          )
          expect_equal(
            nrow(validationGroups[grepl("A3A",validationGroups$foundPatterns),]),
            8/2
          )
          expect_equal(
            nrow(validationGroups[grepl("A3B",validationGroups$foundPatterns),]),
            8/2
          )
          expect_equal(
            nrow(validationGroups[grepl("A3[^GFAB]",validationGroups$foundPatterns),]),
            32/2
          )
         })

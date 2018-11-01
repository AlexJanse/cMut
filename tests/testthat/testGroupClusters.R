context("Group clusters and show mutation symbols")

testGroupClusters <- groupClusters(identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T),
                                   patternIntersect = T)
validationGroups <- groupClusters(identifyAndAnnotateClusters(validationTable,20000,
                                                              sampleIdHeader = "id",
                                                              linkPatterns = 2),
                                  patternIntersect = T)
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
            8/2
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

x <- testDataSet
names(x) <- c("chr","pos","endPos","reference","variant","id","sur")
x$sur <- gsub("\\.","-",x$sur)
testPatterns <- mutationPatterns
names(testPatterns) <- c("id","refer","alternative","sur","dist","reference")
data <- identifyAndAnnotateClusters(dataTable = x,
                                    chromHeader = "chr",
                                    positionHeader = "pos",
                                    refHeader = "reference",
                                    altHeader = "variant",
                                    sampleIdHeader = "id",
                                    contextHeader = "sur",
                                    mutationSymbol = "-",
                                    maxDistance = 20000,
                                    asTibble = F,
                                    linkPatterns = T,
                                    reverseComplement = F,
                                    searchPatterns = testPatterns,
                                    searchAltHeader = "alternative",
                                    searchContextHeader = "sur",
                                    searchRefHeader = "refer",
                                    searchIdHeader = "id",
                                    searchDistanceHeader = "dist",
                                    searchMutationSymbol = ".",
                                    searchReverseComplement = T,
                                    linkClustersOnly = T)
names(data)[8] <- "clId"
testGroupClusters2 <- groupClusters(data,
                                   clusterIdHeader = "clId",
                                   refHeader = "reference",
                                   altHeader = "variant",
                                   asTibble = F,
                                   searchClusterPatterns = T,
                                   searchPatterns = testPatterns,
                                   searchAltHeader = "alternative",
                                   searchRefHeader = "refer",
                                   searchIdHeader = "id",
                                   searchDistanceHeader = "dist",
                                   searchReverseComplement = T,
                                   patternIntersect = T)
testGroupClusters <- groupClusters(identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T),
                                   patternIntersect = T,
                                   searchClusterPatterns = T,
                                   asTibble = F)

test_that("check if changing the default parameters gives the same result",{
          expect_equal(all(testGroupClusters2[,-c(2,8,9)] == testGroupClusters[,-c(2,8,9)]),
                       TRUE)
          expect_equal(class(testGroupClusters2),
                       class(data.frame()))
          expect_equal(testGroupClusters[,9][[1]][1],
                       "PolZetaPrimair")
})

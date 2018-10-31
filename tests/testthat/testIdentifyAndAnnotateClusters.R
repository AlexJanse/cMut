context("Identify and annotate clusters")

testResults <- identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T)
validationResults <- identifyAndAnnotateClusters(validationTable,20000, sampleIdHeader = "id" ,linkPatterns = 2, patternsAsList = F)
test_that("Check if the linked patterns are as expected",{
  expect_equal(
    nrow(testResults[testResults$is.linked == T, ]) == 26,
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

testPatterns <- data.frame(id = c("testPat1","testPat2"),reference = c("A","G"),alternative = c("T","C"), dist = c(NA,NA),sur = c("-","-"))
x <- testDataSet
names(x) <- c("chr","pos","endPos","reference","variant","id","sur")
x$sur <- gsub("\\.","-",x$sur)
testResults2 <- identifyAndAnnotateClusters(dataTable = as.data.frame(x),
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
                                           reverseComplement = T,
                                           searchPatterns = testPatterns,
                                           searchAltHeader = "alternative",
                                           searchContextHeader = "sur",
                                           searchRefHeader = "reference",
                                           searchIdHeader = "id",
                                           searchDistanceHeader = "dist",
                                           searchReverseComplement = F,
                                           linkClustersOnly = F)
test_that("Check if the parameters works",{
          expect_equal(nrow(testResults),
                       1873)
          expect_equal(nrow(testResults[testResults$is.linked,]),
                       326)
          expect_equal(all(testResults[grepl("testPat1",testResults$linkedPatterns),]$ref == "T"),
                       TRUE)
          expect_equal(class(testResults),
                       class(data.frame()))
          expect_equal(all(testResults[testResults$is.clustered,"distance"] <= 20000),
                       TRUE)
          })

testPatterns <- mutationPatterns
names(testPatterns) <- c("id","refer","alternative","sur","dist","reference")
testResults2 <- identifyAndAnnotateClusters(dataTable = x,
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
names(testResults2) <- names(testResults)
test_that("Check if using the default gives the same results as setting",
          expect_equal(all.equal(testResults[,-11],testResults2[,-11]))
          )

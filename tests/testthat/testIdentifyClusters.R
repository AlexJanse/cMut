context("Identify and annotate clusters")

testResults <- identifyClusters(testDataSet, 20000, linkPatterns = T)
validationResults <- identifyClusters(validationTable, 20000, sampleIdHeader = "id" ,linkPatterns = T)
test_that("Check if the linked patterns are as expected",{
  expect_equal(
    nrow(testResults[testResults$is.linked == T, ]),
    61)
  expect_equal(
    testResults[testResults$sampleIDs == "TEST","linkedPatterns"][1,][[1]][[1]],
    "AID")
  expect_equal(
    nrow(validationResults[grepl("AID",validationResults$linkedPatterns),]),
    156
  )
  expect_equal(
    nrow(validationResults[grepl("PolEta",validationResults$linkedPatterns),]),
    8
  )
  expect_equal(
    nrow(validationResults[grepl("A1/A3F",validationResults$linkedPatterns),]),
    44
  )
  expect_equal(
    nrow(validationResults[grepl("A3G",validationResults$linkedPatterns),]),
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
testResults2 <- identifyClusters(dataTable = as.data.frame(x),
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
                                 searchMutationSymbol = "-",
                                 searchDistanceHeader = "dist",
                                 searchReverseComplement = F,
                                 linkClustersOnly = F)
testResults3 <- identifyClusters(dataTable = as.data.frame(x),
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
                                searchReverseComplement = T,
                                renameReverse = T)
test_that("Check if the parameters works",{
          expect_equal(nrow(testResults2),
                       1979)
          expect_equal(nrow(testResults2[testResults2$is.linked,]),
                       342)
          expect_equal(all(testResults2[grepl("testPat1",testResults2$linkedPatterns),]$reference == "T"),
                       TRUE)
          expect_equal(class(testResults2),
                       class(data.frame()))
          expect_equal(all(testResults2[testResults2$is.clustered,"distance"] <= 20000),
                       TRUE)
          expect_equal(nrow(testResults3[grepl("\\[Rev\\.Com\\.\\]",testResults3$linkedPatterns),]),
                       21)
          })

testPatterns <- mutationPatterns
names(testPatterns) <- c("id","refer","alternative","sur","dist","reference")
testResults2 <- identifyClusters(dataTable = x,
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
test_that("Check if using the default gives the same results as setting own parameters.",
          expect_equal(all(testResults[,-c(7,11)] == testResults2[,-c(7,11)]),
                       TRUE)
          )

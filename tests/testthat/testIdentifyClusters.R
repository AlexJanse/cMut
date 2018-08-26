context("Identify clusters")
library(tidyverse)

wd <- substr(getwd(),1,nchar(getwd())-15) # Get main work directory path

testData = as.tibble(
  dget(
    paste(wd,"/Data/testData.txt",sep = ""),
    keep.source = FALSE)
  ) # Tibble containing filter test data
testSample = testData[which(
  testData$sampleID == "102-00001-03" &
    testData$Chr == "chr1"),] # A tibble with one chromosome data from one sample

test_that("A vector with sample names and cluster IDs",
  expect_equal(identifyClusters(testSample,10^10),
               dget(paste(wd,"/Data/expectedResult1.txt",sep = "")))
)

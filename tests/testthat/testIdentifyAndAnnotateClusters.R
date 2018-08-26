context("Identify and annotate clusters")
library(tidyverse)

wd <- substr(getwd(),1,nchar(getwd())-15) # Get main work directory path

testData = as.tibble(
  dget(
    paste(wd,"/Data/testData.txt",sep = ""),
    keep.source = FALSE)
) # Tibble containing filter test data

test_that("Check if the results are the same as the original file",
          expect_equal(
            isTRUE(all_equal(
              identifyAndAnnotateClusters(testData,19),
              as.tibble(dget(paste(wd,"/Data/expectedResult2.txt",sep = ""))))),
            TRUE)
          )

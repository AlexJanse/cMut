context("Test the shuffle function")

set.seed(1)
identResults <- identifyAndAnnotateClusters(testDataSet,20000,linkPatterns = T)
test_that("Check if the shuffle results are the same",
  expect_equal(
    all.equal(
      shuffleMutations(identResults[identResults$is.clustered == T,],nBootstrap = 2, tibble = F),
      as.data.frame(dget("extdata/expectedResultsShuffleMutations.txt",keep.source = F))),
    TRUE)
)

context("Generate random data")

set.seed(1)
test_that("To get the same results as the test data",
          expect_equal(createRandomMutations(),
                       dget("extdata/expectedResultsRandom.txt", keep.source = F)))
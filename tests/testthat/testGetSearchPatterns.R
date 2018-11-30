context("Test getSearchPatterns")

test1 <- getSearchPatterns(asTibble = F)
test2 <- getSearchPatterns(reverse = F,asTibble = F)
test3 <- getSearchPatterns(reverse = T, renameReverse = T,asTibble = F)
test4 <- getSearchPatterns(reverse = F, renameReverse = T,asTibble = F)

test_that("See if the function returns the expected tables",{
  expect_equal(dplyr::all_equal(test2,
               test4),
               dplyr::all_equal(test4, mutationPatterns))
  expect_equal(nrow(test1),
               nrow(test3))
  expect_equal(nrow(test2)+11,
               nrow(test3))

})

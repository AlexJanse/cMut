context("Check shuffleMutation")

test <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 2)
while(test[test$process == "Unidentified","percentage"][[1]] == 100){
  test <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 2)
}
test2 <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = F, no.cores = 2)
while(test2[test2$process == "Unidentified","percentage"][[1]] == 100){
  test2 <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = F, no.cores = 2)
}
test_that("Check if percentage is correct",{
  expect_equal(sum(test2$percentage) >= 100,
               TRUE)
  expect_equal(nrow(test2) == 8,
               TRUE)
  expect_equal(sum(test$percentage) >= 100,
               TRUE)
  expect_equal(nrow(test) == 11,
               TRUE)
})

test_that("Check if parameters checked correctly",{
  expect_error(shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 20))
  expect_error(shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = -1))
  expect_error(shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 2, altHeader = "A"))
  expect_error(shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 2,searchAltHeader = "e"))
})

test <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T, no.cores = 2, returnEachBootstrap = T, renameReverse = T)
test_that("Check if returnEachBootstrap and renameReverse are working properly",{
  expect_equal(length(test),
               5)
  expect_equal(class(test),
               class(list()))
  expect_equal(nrow(test[[1]]),
               22)
})

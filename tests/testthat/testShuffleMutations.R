context("Check shuffleMutation")

test <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T)
while(test[test$process == "Unidentified","percentage"][[1]] == 100){
  test <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = T)
}
test2 <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = F)
while(test2[test2$process == "Unidentified","percentage"][[1]] == 100){
  test2 <- shuffleMutations(testDataSet,nBootstrap = 5, searchClusterPatterns = F)
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

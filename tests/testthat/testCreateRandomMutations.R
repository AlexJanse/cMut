context("Generate random data")

set.seed(1)
test_that("To get the same results as the test data",{
          expect_equal(createRandomMutations(),
                       dget("extdata/expectedResultsRandom.txt", keep.source = F))
  })

test <- createRandomMutations(nMut = 2,
                              sizeSur = 4,
                              asTibble = F,
                              refGenomeHg19 = F,
                              useChrom = c("chr1","chr2"))
test_that("Check if parameters are correctly processed",{
  expect_equal(nrow(createRandomMutations(nMut = 10)),
               10)
  expect_equal(nchar(test[1,"surrounding"][[1]]),
               9)
  expect_equal(class(test),
               class(data.frame()))
  expect_equal(test[1,"ref"][[1]],
               as.character(BSgenome::getSeq(
                 BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                 GenomicRanges::GRanges(test[1,"chrom"][[1]],
                                        IRanges::IRanges(test[1,"start"][[1]],
                                                         test[1,"start"][[1]])
                                        )
                 ))
               )
  expect_error(createRandomMutations(nMut = 2,
                                     sizeSur = 4,
                                     asTibble = F,
                                     refGenomeHg19 = F,
                                     useChrom = c("chr1","ch2")),
               "chromosome names")
})

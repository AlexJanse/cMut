#' GenrateRandomDNMs
#' @description Function to generate random data in a tibble with DNM information.
#' For explanation of the columns use cat(comment(<returned object>))
#' @param nDNM number of mutations that needed to be generated
#' @param sampleName A name for the test sample
#' @return A tibble with DNM information
#' @export
generateRandomDNMs <- function(nDNM = 100, sampleName = "testSample"){

  # Check arguments ---------------------------------------------------------
  stopifnot(nDNM > 1)

  # Human chromosomes information -------------------------------------------
  nucleotides <- c("A","C","G","T")
  nameChrom <- c("chr1","chr2","chr3","chr4",
                 "chr5","chr6","chr7","chr8",
                 "chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16",
                 "chr17","chr18","chr19","chr20",
                 "chr21","chr22","chrX","chrY")
  lenChrom <- c(248956422,	242193527,	198295559,	190214555,
                181538259,	170805979,	159345973,	145138636,
                138394717,	133797422,	135086622,	133275309,
                114364328,	107043718,	101991189,	90338345,
                83257441, 80373285,	58617616,	 64444167,
                46709983,	50818468, 156040895,	57227415)
  names(lenChrom) <- nameChrom
  probability <- lenChrom/sum(lenChrom)*100

  # Create random data ------------------------------------------------------
  randomTable <- data.frame(chrom = sample(nameChrom, nDNM, replace = T, prob = probability),
             stringsAsFactors = F) %>%
    mutate(chromLen = lenChrom[chrom]) %>%
    mutate(start = map_int(chromLen, ~sample(., 1))) %>%
    mutate(stop = start) %>%
    mutate(ref = sample(nucleotides, nDNM, replace = T)) %>%
    mutate(alt = map_chr(ref, ~sample(setdiff(nucleotides, c(.)),1))) %>%
    mutate(sampleIDs = sampleName) %>%
    mutate(surrounding = map_chr(sampleIDs,
                                 ~paste(sample(nucleotides,1),
                                        sample(nucleotides,1),
                                        sep = ".")))

  randomTable$chromLen <- NULL


  randomTable <- as.tibble(randomTable)
  comment(randomTable) <-
    " A random generated tibble with the following information
    chrom : Name of the chromosome
    start : start position of the mutation,
    (random generated but with the borders of the human reference genome)
    stop : stop position of the mutation
    ref : the nucleotide on the reference (also random generated)
    alt : the nucleotife that the sample has
    sampleName : name of the random sample
    surrounding : the direct linked nucleotides surrounding the mutation"
  return(randomTable)
}

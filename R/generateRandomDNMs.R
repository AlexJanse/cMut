#' GenrateRandomDNMs
#' @description Function to generate random tibble with DNM information.
#' For explanation of the columns use cat(comment(<returned object>))
#' @param nSamples number of samples
#' @param nChrom number of chromosomes per sample
#' @param nDNM number of mutations per chromosome
#' @return A tibble with DNM information
#' @export
generateRandomDNMs <- function(nSamples = 2, nChrom = 2, nDNM = 10){

  # Check arguments ---------------------------------------------------------
  stopifnot(nSamples > 1 && nChrom > 1 && nDNM > 1)

  # Human chromosomes information -------------------------------------------
  nucleotides <- c("A","C","G","T")
  nameChrom <- c("1","2","3","4",
                 "5","6","7","8",
                 "9","10","11","12",
                 "13","14","15","16",
                 "17","18","19","20",
                 "21","22","X","Y")
  lenChrom <- c(224999719,	237712649,	194704827,	187297063,
                177702766,	167273993,	154952424,	142612826,
                120312298,	131624737,	131130853,	130303534,
                95559980,	88290585,	81341915,	78884754,	77800220,
                74656155,	55785651,	59505254,	34171998,	34893953,
                151058754,	25121652)
  names(lenChrom) <- nameChrom

  # Create random data ------------------------------------------------------
  randomTable <- data.frame(chrom = character(),
                            start = integer(),
                            end = integer(),
                            ref = character(),
                            alt = character(),
                            sampleName = character(),
                            surrounding = character(),
                            stringsAsFactors = FALSE
                            )

  for (nrSample in seq.int(1, round(nSamples))) {
    sampleName = paste("randomSample-", nrSample, sep = "")

    for(nrChrom in seq.int(1, round(nChrom))) {
      chrom = sample(nameChrom, 1)
      size = lenChrom[chrom]

      for(nrDNM in seq.int(1, round(nDNM))) {
        start = sample(seq.int(1, size), 1)

        if(start == 1) {
          start <- start+1
        } else if(start == size) {
          start <- start-1
        } # To make sure that with the surrounding nucleotides we will stay within that size of the chromosome

        ref = sample(nucleotides, 1)
        alt = sample(setdiff(nucleotides, c(ref)), 1)
        surrounding = paste(sample(nucleotides, 1),sample(nucleotides, 1), sep = ".")

        randomRow <- data.frame(chrom = paste("chr",chrom, sep = ""),
                                  start = start,
                                  end = start,
                                  ref = ref,
                                  alt = alt,
                                  sampleName = sampleName,
                                  surrounding = surrounding,
                                  stringsAsFactors = FALSE)

        randomTable <- rbind.data.frame(randomTable, randomRow)
      }
    }
  }

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

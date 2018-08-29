#' GenerateRandomMutations
#' @description Function to generate random data in a tibble with DNA mutation information.
#' For explanation of the columns use cat(comment(<returned object>))
#' @param nMut number of mutations that needed to be generated
#' @param sampleName A name for the test sample
#' @return A tibble with mutation information
#' @export
#' @import magrittr
generateRandomMutations <- function(nMut = 500, sampleName = "testSample"){

  # Check arguments ---------------------------------------------------------
  stopifnot(nMut > 1)

  # Human chromosomes information -------------------------------------------
  nucleotides <- c("A","C","G","T")
  nameChrom <- c("chr1","chr2","chr3","chr4",
                 "chr5","chr6","chr7","chr8",
                 "chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16",
                 "chr17","chr18","chr19","chr20",
                 "chr21","chr22","chrX","chrY")
  lenChrom <- c(249250621, 243199373, 198022430, 191154276,
                180915260, 171115067, 159138663, 146364022,
                141213431, 135534747, 135006516, 133851895,
                115169878, 107349540, 102531392, 90354753,
                81195210, 78077248, 59128983, 63025520,
                48129895, 51304566, 155270560, 59373566)
  names(lenChrom) <- nameChrom
  probability <- lenChrom/sum(lenChrom)*100

  # Create random data ------------------------------------------------------
  randomTable <- data.frame(chrom = sample(nameChrom, nMut, replace = T, prob = probability),
                            stringsAsFactors = F) %>%
    plyr::mutate(chromLen = lenChrom[chrom]) %>%
    plyr::mutate(start = purrr::map_int(chromLen, ~sample(., 1))) %>%
    plyr::mutate(stop = start) %>%
    plyr::mutate(irange = paste(chrom,start,stop,sep = "-")) %>%
    plyr::mutate(ref = purrr::map_chr(irange, ~getRef(.))) %>%
    plyr::mutate(alt = purrr::map_chr(ref, ~sample(setdiff(nucleotides, c(.)),1))) %>%
    plyr::mutate(sampleIDs = sampleName) %>%
    plyr::mutate(surrounding = paste(purrr::map2_chr(chrom,start-1,getSurrounding, table = .),
                                     purrr::map2_chr(chrom,start+1,getSurrounding, table = .),
                                     sep = "."))

  randomTable$irange <- NULL
  randomTable$chromLen <- NULL

  randomTable <- tibble::as.tibble(randomTable)
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


#' getRef
#' @description A function to get the reference nucleotide
#' @param x string that contains the chromosome name, start and stop location and is seperated by "-"
getRef <- function(x){
  data <- unlist(strsplit(x,"\\-"))
  chr <- as.character(data[1])
  start <- as.numeric(data[2])
  stop <- as.numeric(data[3])
  range <- GenomicRanges::GRanges(chr,IRanges::IRanges(start,stop))
  return(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,range)))
}

getSurrounding <- function(chr,loc,table){
  sameLoc <- table[table$start == loc,]
  if(nrow(sameLoc) > 0){
    altSurr <- sameLoc[sameLoc$chrom == chr,]

    if(nrow(altSurr) > 0){
      stopifnot(nrow(altSurr) == 1)
      return(altSurr[1,"alt"])
    }
  }
  return(getRef(paste(chr,loc,loc,sep = "-")))
}

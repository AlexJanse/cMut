#' createRandomMutations
#' @description Function to generate random data in a tibble with DNA mutation
#'   information. For explanation of the columns use cat(comment(<returned
#'   object>))
#' @param nMut number of mutations that needed to be generated
#' @param sampleName A name for the test sample
#' @param tibble A boolean if the table has to be a tibble
#' @param sizeSur A number with the ammount of nucleotides left and right of the
#'   mutation. (e.g. sizeSur = 2; CC.GT)
#' @return A tibble with mutation information
#' @export
#' @import magrittr
#' @examples
#' x <- createRandomMutations(nMut = 100,
#'                            sampleName = "test",
#'                            sizeSur = 3)
#'
#' # See explanation of table columns
#' cat(comment(x))
createRandomMutations <- function(nMut = 500,
                                  sampleName = "testSample",
                                  tibble = TRUE,
                                  sizeSur = 2){

  stopifnot(sizeSur > 1)

  # Human chromosomes information (GRCh37/hg19) ------------------------------
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

  # Check argument ----------------------------------------------------------
  stopifnot(nMut > 1)

  # Create random data ------------------------------------------------------
  probability <- lenChrom/sum(lenChrom)*100
  lenChromEnquo <- dplyr::enquo(lenChrom)
  sizeSur <- dplyr::enquo(sizeSur)

  randomTable <- data.table::data.table(chrom = sample(nameChrom, nMut, replace = T, prob = probability),
                            stringsAsFactors = F)
  randomTable <- dplyr::mutate(randomTable, chromLen = lenChrom[chrom])
  randomTable <- dplyr::mutate(randomTable, start = purrr::map_int(chromLen, ~sample(., 1)))
  randomTable <- dplyr::mutate(randomTable, stop = start)
  randomTable <- dplyr::mutate(randomTable, irange = paste(chrom,start,stop,sep = "-"))
  randomTable <- dplyr::mutate(randomTable, sampleIDs = sampleName)
  randomTable <- dplyr::mutate(randomTable, refdata = purrr::map2_chr(chrom,start,function(x,y){getRefData(x,y,sizeSur, lenChrom = !!lenChromEnquo)}))
  randomTable <- dplyr::mutate(randomTable, surrounding = purrr::map_chr(refdata,function(x){strsplit(x,"-")[[1]][1]}))
  randomTable <- dplyr::mutate(randomTable, ref = purrr::map_chr(refdata,function(x){strsplit(x,"-")[[1]][2]}))
  randomTable <- dplyr::mutate(randomTable, alt = purrr::map_chr(ref, ~sample(setdiff(nucleotides, c(.)),1)))

  randomTable$refdata <- NULL
  randomTable$irange <- NULL
  randomTable$chromLen <- NULL

  randomTable <- randomTable[c("chrom","start","stop","ref","alt","sampleIDs","surrounding")]

  if(tibble){
    randomTable <- tibble::as.tibble(randomTable)
  }

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
#' @param x A string that contains the chromosome name, start and stop location
#'   and is seperated by "-"
#' @param lenChrom A vector with chromosomes length
getRef <- function(x,lenChrom){
  data <- unlist(strsplit(x,"\\-"))
  chr <- as.character(data[1])
  start <- as.numeric(data[2])
  stop <- as.numeric(data[3])
  if(start < 1 | stop > lenChrom[chr]){
    return("")
  }
  range <- GenomicRanges::GRanges(chr,IRanges::IRanges(start,stop))
  return(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,range)))
}

#' getRefData
#' @description Function to get the surrounding of a mutation
#' @param chr A string with the chromosome name
#' @param pos A number with the mutation location
#' @param sizeSur A number with the amount of nucleotides there has to be
#'   arround the mutation.
#' @param table A table containing the mutation information
#' @inherit getRef
getRefData <- function(chr,pos,sizeSur,lenChrom){
  lenChrom <- rlang::get_expr(lenChrom)
  sizeSur <- rlang::get_expr(sizeSur)
  maxPos <- lenChrom[chr]
  start <- pos-sizeSur
  stop <- pos+sizeSur

  if(start < 1){
    start <- 1
  }
  if(stop > maxPos){
    stop <- maxPos
  }
  context <- getRef(paste(chr,start,stop,sep = "-"),lenChrom)
  if(pos == start){
    mutPos <- 1
  } else if(pos == maxPos){
    mutpos <- nchar(context)
  } else{
    mutpos <- pos-start+1
  }
  refData <- paste(paste(substr(context,1,mutpos-1),
                       substr(context,mutpos+1,nchar(context)),
                       sep = "."),substr(context,mutpos,mutpos),sep = "-")
  return(refData)
}

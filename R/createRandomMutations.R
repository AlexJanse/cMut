#' createRandomMutations
#' @description Function to generate random mutation data.
#' @param nMut The number of mutations that needed to be generated.
#' @param sampleName A name for the test sample.
#' @param asTibble A Boolean to tell if the returned table needs to be a tibble.
#'   Returns a data.frame when FALSE.
#' @param sizeSur A number with the amount of nucleotides left and right of the
#'   mutation. (e.g. sizeSur = 2 --> "CC.GT")
#' @param refGenomeHg19 A Boolean to tell if the preferred reference genome has
#'   to be Hg19. Hg38 will be used when \code{refGenome = FALSE}.
#' @param useChrom A vector with the chromosomes that are used. Make sure that
#'   the names are notated as e.g. "chr1".
#' @return A tibble or data.frame with random mutation data.
#' @export
#' @import magrittr
#' @examples
#' x <- createRandomMutations(nMut       = 100,
#'                            sampleName = "test",
#'                            sizeSur    = 3)
#'
#' # Proper way of using all the parameters:
#' createRandomMutations(nMut          = 10,
#'                       sampleName    = "name",
#'                       asTibble      = FALSE,
#'                       sizeSur       = 4,
#'                       refGenomeHg19 = FALSE,
#'                       useChrom      = c("chr1", "chr2"))
#'
#' # See explanation of table columns
#' cat(comment(x))
createRandomMutations <- function(nMut,
                                  sampleName = "testSample",
                                  asTibble = TRUE,
                                  sizeSur = 2,
                                  refGenomeHg19 = TRUE,
                                  useChrom = c("chr1","chr2","chr3","chr4",
                                               "chr5","chr6","chr7","chr8",
                                               "chr9","chr10","chr11","chr12",
                                               "chr13","chr14","chr15","chr16",
                                               "chr17","chr18","chr19","chr20",
                                               "chr21","chr22","chrX","chrY")){

  # Human chromosomes information (hg19 or hg38) ------------------------------
  nameChrom <- c("chr1",  "chr2",  "chr3",  "chr4",
                 "chr5",  "chr6",  "chr7",  "chr8",
                 "chr9",  "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16",
                 "chr17", "chr18", "chr19", "chr20",
                 "chr21", "chr22", "chrX",  "chrY")
  if (refGenomeHg19) {
    lenChrom <- c(249250621, 243199373, 198022430, 191154276,
                  180915260, 171115067, 159138663, 146364022,
                  141213431, 135534747, 135006516, 133851895,
                  115169878, 107349540, 102531392, 90354753,
                  81195210,  78077248,  59128983,  63025520,
                  48129895,  51304566,  155270560, 59373566)
  } else {
    lenChrom <- c(248956422, 242193529, 198295559, 190214555,
                  181538259, 170805979, 159345973, 145138636,
                  138394717, 133797422, 135086622, 133275309,
                  114364328, 107043718, 101991189, 90338345,
                  83257441,  80373285,  58617616,  64444167,
                  46709983,  50818468,  156040895, 57227415)
  }
  names(lenChrom) <- nameChrom

  # Check argument ----------------------------------------------------------
  stopifnot(nMut > 0)
  stopifnot(sizeSur > 0)
  if (!all(startsWith(useChrom, "chr"))){
    stop ("Make sure that the chromosome names starts with chr.
          See the default value of useChrom in the createRandomMutations help
          page for an example.")
  }

  # Create randomTable ------------------------------------------------------
  randomTable <- createRandomTable(nMut          = nMut,
                                   sampleName    = sampleName,
                                   sizeSur       = sizeSur,
                                   refGenomeHg19 = refGenomeHg19,
                                   useChrom      = useChrom,
                                   lenChrom      = lenChrom)

  # Repeat the rows with "N" as reference -----------------------------------
  while (nrow(randomTable[randomTable$ref == "N",]) > 0) {
    randomTable <- dplyr::bind_rows(randomTable[randomTable$ref != "N",],
                                    createRandomTable(nMut          = nrow(randomTable[randomTable$ref == "N",]),
                                                      sizeSur       = sizeSur,
                                                      sampleName    = sampleName,
                                                      refGenomeHg19 = refGenomeHg19,
                                                      lenChrom      = lenChrom,
                                                      useChrom      = useChrom))
  }

  # Make sure there are no repeats in the table ------------------------------
  randomTable <- unique(randomTable)
  while (nrow(randomTable) != nMut) {
    randomTable <- dplyr::bind_rows(randomTable[randomTable$ref != "N",],
                                    createRandomTable(nMut          = nMut-nrow(randomTable),
                                                      sizeSur       = sizeSur,
                                                      sampleName    = sampleName,
                                                      refGenomeHg19 = refGenomeHg19,
                                                      lenChrom      = lenChrom,
                                                      useChrom      = useChrom))
    randomTable <- unique(randomTable)
  }

  # Convert if wished
  if (asTibble) {
    randomTable <- tibble::as.tibble(randomTable)
  } else {
    randomTable <- as.data.frame(randomTable)
  }

  comment(randomTable) <-
    paste0(" A random generated ",
ifelse(asTibble,"tibble","data.table")," with the following
                information
  chrom       : Name of the chromosome
  start       : Start position of the mutation,
                (random generated but with the borders of
                the human reference genome)
  stop        : Stop position of the mutation. Always the
                same as start.
  ref         : The nucleotide on the reference genome ",
ifelse(refGenomeHg19,"Hg19","Hg38"),".","
  alt         : The nucleotife that the sample has
  sampleName  : Name of the random sample
  surrounding : The direct linked nucleotides surrounding
                the mutation")

  return(randomTable)
}


#' getRef
#' @description A function to get the reference nucleotide.
#' @param x A string that contains the chromosome name, start and stop location
#'   and is separated by "-".
#' @param lenChrom A vector with chromosomes length
#' @inheritParams createRandomMutations
#' @return A string with the reference nucleotides.
getRef <- function(x,lenChrom, refGenomeHg19) {


  # Prepare parameters ------------------------------------------------------
  data  <- unlist(strsplit(x,"\\-"))
  chr   <- as.character(data[1])
  start <- as.numeric(data[2])
  stop  <- as.numeric(data[3])


  # Check if the start and stop positions are not outside the range ---------
  if (start < 1 | stop > lenChrom[chr]) {
    return("")
  }


  # Create Range object and return the reference nucleotide ------------------
  range <- GenomicRanges::GRanges(chr,IRanges::IRanges(start,stop))
  if (refGenomeHg19) {
    return(as.character(BSgenome::getSeq(
                          BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                          range)))
  } else {
    return(as.character(BSgenome::getSeq(
                          BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                          range)))
  }
}


#' getRefData
#' @description Function to get the reference nucleotides of a mutation.
#' @param chr A string with the chromosome name.
#' @param pos A number with the mutation location.
#' @param sizeSur A number with the amount of nucleotides there has to be
#'   around the mutation.
#' @inheritParams createRandomMutations
#' @return A string with the surrounding and mutated position reference
#'   nucleotides of a mutation. Each part separated with "-"
#' @inheritParams getRef
getRefData <- function(chr,pos,sizeSur,lenChrom, refGenomeHg19) {

  # Prepare parameters ------------------------------------------------------
  lenChrom <- rlang::get_expr(lenChrom)
  sizeSur  <- rlang::get_expr(sizeSur)
  maxPos   <- lenChrom[chr][[1]]
  start    <- pos-sizeSur
  stop     <- pos+sizeSur


  # Adjust the start stop parameters if it outside the range ----------------
  if (start < 1) {
    start <- 1
  }
  if (stop > maxPos) {
    stop <- maxPos
  }


  # Get the reference nucleotides -------------------------------------------
  context <- getRef(paste(chr,start,stop,sep = "-"),lenChrom,refGenomeHg19)


  # Determine the location of the mutated nucleotide ------------------------
  if(pos == start){
    mutPos <- 1
  } else if(pos == maxPos){
    mutpos <- nchar(context)
  } else{
    mutpos <- pos-start+1
  }


  # Create reference symbol and return it ----------------------------------
  refData <- paste(paste(substr(x     = context,
                                start = 1,
                                stop  = mutpos-1),
                         substr(x     = context,
                                start = mutpos+1,
                                stop  = nchar(context)),
                         sep = "."),
                   substr(x     = context,
                          start = mutpos,
                          stop  = mutpos),
                   sep = "-")
  return(refData)
}


#' createRandomTable
#' @description Function to build the random mutation table and call supporting
#'   functions.
#' @param lenChrom vector with the length of the chromosomes
#' @inheritParams createRandomMutations
#' @return A table with random mutations
#' @importFrom rlang .data
createRandomTable <- function(nMut,
                              sampleName,
                              sizeSur,
                              refGenomeHg19,
                              useChrom,
                              lenChrom) {

  # Prepare parameters ------------------------------------------------------
  probability   <- lenChrom[useChrom]/sum(lenChrom[useChrom])*100
  lenChromEnquo <- dplyr::enquo(lenChrom)
  sizeSur       <- as.numeric(sizeSur)
  sizeSur       <- dplyr::enquo(sizeSur)


  # Nucleotides to choose from ----------------------------------------------
  nucleotides <- c("A", "C", "G", "T")


  # Build the table (no pipe to save time) ----------------------------------

  # Start with creating column with random chromosomes:
  randomTable <- data.table::data.table(chrom = sample(useChrom,
                                                       nMut,
                                                       replace = TRUE,
                                                       prob    = probability),
                                        stringsAsFactors = FALSE)

  # Add information about the chromosome maximum size:
  randomTable <- dplyr::mutate(randomTable, chromLen = lenChrom[.data$chrom])

  # Take a random number with the range of the chromosome:
  randomTable <- dplyr::mutate(randomTable, start = purrr::map_int(.data$chromLen,
                                                                   ~sample(., 1)))
  # Stop is the same as start because thease are single mutations only:
  randomTable <- dplyr::mutate(randomTable, stop = .data$start)

  # Create a column with the necessary data for the getRefData function:
  randomTable <- dplyr::mutate(randomTable, irange = paste(.data$chrom, .data$start,
                                                           stop, sep = "-"))

  # Create a column with the sample name:
  randomTable <- dplyr::mutate(randomTable, sampleIDs = sampleName)

  # Create a column with the information needed from the reference genome:
  randomTable <- dplyr::mutate(randomTable,
                               refdata = purrr::map2_chr(.data$chrom,.data$start,
                                                         function(x,y){
                                                           getRefData(chr = x,
                                                                      pos = y,
                                                                      sizeSur = sizeSur,
                                                                      lenChrom = !!lenChromEnquo,
                                                                      refGenomeHg19 = refGenomeHg19)
                                                         }))

  # Create a column with the surrounding nucleotides from the refData column:
  randomTable <- dplyr::mutate(randomTable,
                               surrounding = purrr::map_chr(.data$refdata,
                                                            function(x){
                                                              strsplit(x,
                                                                       "-")[[1]][1]
                                                            }))

  # Create a column with the reference nucleotide from the refData column:
  randomTable <- dplyr::mutate(randomTable, ref = purrr::map_chr(.data$refdata,
                                                                 function(x){
                                                                   strsplit(x,
                                                                            "-")[[1]][2]
                                                                 }))

  # Create a column with the variant nucleotide:
  randomTable <- dplyr::mutate(randomTable, alt = purrr::map_chr(.data$ref,
                                                                 ~sample(setdiff(nucleotides,
                                                                                 c(.)),1)))


  # Return the table with only the columns that are needed ------------------
  return(randomTable[c("chrom", "start", "stop",
                       "ref",   "alt",   "sampleIDs",
                       "surrounding")])
}

#' groupClusters
#' @description A function that will group the clusters and give the mutation consensus back
#' @param table A table with columns containing cluster ID's, reference and alternative nucleotide
#' @param clusterIdHeader Contains the name of the column with the cluster IDs
#' @param refHeader Contains the name of the column with the reference nucleotides
#' @param altHeader Contains the name of the column with the alternative nucleotides
#' @param groupPatterns Boolean if the table contains patterns and these needed to be processed aswell
#' @export
#' @import magrittr
#' @import foreach
#' @examples
#' # Example of a table containing the right columns and data
#' test <- createRandomMutations()
#' test <- identifyAndAnnotateClusters(test, 20000, chromHeader = "chrom",
#'                                     sampleIdHeader = "sampleIDs",positionHeader = "start")
groupClusters <- function(table, clusterIdHeader = "clusterId",
                          refHeader = "ref", altHeader = "alt",
                          patterns = FALSE, patternHeader = "linkedPatterns"){

  # Check data --------------------------------------------------------------------------------------------
  stopifnot(nrow(table) > 0)
  stopifnot(!any(is.na(dplyr::select(table,clusterIdHeader,refHeader, altHeader))))

  # Build table -------------------------------------------------------------------------------------------
  table <- table %>%
    dplyr::group_by_(clusterIdHeader) %>%
    tidyr::nest(.key = "cMuts") %>%
    dplyr::filter(clusterId!="") %>%
    dplyr::mutate(refs = purrr::map(cMuts, ~as.character(dplyr::pull(., refHeader))),
           alts = purrr::map(cMuts, ~as.character(dplyr::pull(., altHeader)))) %>%
    dplyr::mutate(plusStrand = purrr::map2_chr(refs, alts, formatClusterMutations),
           minusStrand = purrr::map2_chr(refs, alts, formatClusterMutations, convert=TRUE)) %>%
    dplyr::mutate(clusterType = purrr::map2_chr(plusStrand, minusStrand, getClusterType))

  if(patterns){
    patternHeader <- dplyr::enquo(patternHeader)
    table <- table %>%
      dplyr::mutate(patternInternsect = purrr::map(cMuts,groupPatterns,!!patternHeader))
  }

  if(nrow(table) == 0){
    warning("No rows found. Please make sure the cluster IDs are present and try again.")
  }

  return(table)

}

#' convertLetter
#' @description Function to get the reverse complement nucleotides
#' @param chars List with nucleotides
convertLetter <- function(chars){

  res <- foreach::foreach(char = chars,
                 .combine = c) %do% {
                   as.character(Biostrings::reverseComplement(Biostrings::DNAString(char)))
                 }
  return(res)
}

#' formatClusterMutations
#' @description Function to make a mutation description symbol (e.g. G>C)
#' @param refs A list containing the reference nucleotides
#' @param alts A list containing the alternative nucleotides
#' @param convert A boolean that tells if the nucleotides need to be converted to the reverse complement
formatClusterMutations <- function(refs, alts, convert=FALSE) {
  # Check parameters -----------------------------------------------------------------
  stopifnot(length(refs) == length(alts))
  stopifnot(is.logical(convert))

  # Create the description symbols per mutation-------------------------------------
  combinedResult <-
    foreach::foreach(i = 1:length(refs),
            .combine = c) %do% {
              if(convert){
                res <- paste0(convertLetter(refs[i]), ">", convertLetter(alts[i]))
              } else {
                res <- paste0(refs[i], ">", alts[i])
              }
            }

  # Combine the mutation pairs -----------------------------------------------------
  finalResult <- dplyr::if_else(convert,
                         do.call(paste, as.list(rev(combinedResult))),
                         do.call(paste, as.list(combinedResult)))
  return(finalResult)
}

#' getClusterType
#' @description Function to sort and create a cluster type (e.g. C>G C>G / G>C G>C)
#' @param plusstrand A string with the plusStrand mutation description
#' @param minusStrand A string with the minusStrand mutation description
getClusterType <- function(plusStrand, minusStrand) {
  s <- sort(c(plusStrand, minusStrand))
  stringr::str_glue("{s[1]} / {s[2]}")
}

# TODO: controleer of het de juiste gegevens uithaalt en wat er gebeurt bij lege patronen resultaten
groupPatterns <- function(clusterList,patternHeader){
  patternHeader <- rlang::get_expr(patternHeader)

  allPatterns <-
    foreach::foreach(index = 1:nrow(clusterList),
                     .combine = c) %do% {
                       patterns <- clusterList[index,patternHeader][[1]]

                     }
  print(allPatterns)
  readline(prompt = "")
}

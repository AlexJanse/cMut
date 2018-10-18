#' groupClusters
#' @description A function that will group the clusters and give the mutation
#'   consensus back.
#' @param table A table with columns containing cluster ID's, reference and
#'   alternative nucleotide.
#' @param clusterIdHeader Contains the name of the column with the cluster IDs.
#' @param refHeader Contains the name of the column with the reference
#'   nucleotides..
#' @param altHeader Contains the name of the column with the alternative
#'   nucleotides.
#' @param patternIntersect Boolean if the table contains patterns and these
#'   needed to be processed aswell.
#' @param patternHeader A string with the column name of the patterns. Only in
#'   use when patternIntersect is TRUE.
#' @param showWarning A boolean if there need to be a warning if nrow is 0
#' @param searchClusterPatterns A boolean if it's needed to search to cluster
#'   patterns (e.g. GA > TT).
#' @inheritParams identifyAndAnnotateClusters
#' @export
#' @import magrittr
#' @import foreach
#' @examples
#' # Example of a table containing the right columns and data
#' test <- testDataSet
#'
#' # Example of using this function without linking it with patterns
#' mutations <- identifyAndAnnotateClusters(x = test,
#'                                          maxDistance = 20000,
#'                                          chromHeader = "chrom",
#'                                          sampleIdHeader = "sampleIDs",
#'                                          positionHeader = "start",
#'                                          linkPatterns = FALSE)
#' clusters <- groupClusters(table = mutations,
#'                           patternIntersect = FALSE)
#'
#' # Example of using this function with linking it with patterns
#' mutations <- identifyAndAnnotateClusters(x = test,
#'                                          maxDistance = 20000,
#'                                          chromHeader = "chrom",
#'                                          sampleIdHeader = "sampleIDs",
#'                                          positionHeader = "start",
#'                                          linkPatterns = TRUE)
#' clusters <- groupClusters(table = mutations,
#'                           patternIntersect = TRUE)
groupClusters <- function(table,
                          clusterIdHeader = "clusterId",
                          refHeader = "ref",
                          altHeader = "alt",
                          tibble = TRUE,
                          patternIntersect = FALSE,
                          searchClusterPatterns = FALSE,
                          patternHeader = "linkedPatterns",
                          showWarning = TRUE,
                          searchPatterns = NULL, searchRefHeader = "ref",
                          searchAltHeader = "alt",
                          searchIdHeader = "process", searchDistanceHeader = "maxDistance",
                          searchReverseComplement = TRUE){

  # Check data --------------------------------------------------------------------------------------------
  stopifnot(nrow(table) > 0)
  stopifnot(!any(is.na(dplyr::select(table,clusterIdHeader,refHeader, altHeader))))

  # Build table -------------------------------------------------------------------------------------------
  table <- convertFactor(table)
  table <- dplyr::group_by_(table, clusterIdHeader)
  table <- tidyr::nest(table, .key = "cMuts")
  table <- dplyr::filter(table, clusterId!="")
  table <- dplyr::mutate(table, refs = purrr::map(cMuts, ~as.character(dplyr::pull(., refHeader))),
           alts = purrr::map(cMuts, ~as.character(dplyr::pull(., altHeader))))
  table <- dplyr::mutate(table, refs = purrr::map_chr(refs, function(x){paste0(x,collapse = "")}),
                         alts = purrr::map_chr(alts, function(x){paste0(x,collapse = "")}))
  table <- dplyr::mutate(table, plusStrand = purrr::map2_chr(refs, alts, formatClusterMutations),
           minusStrand = purrr::map2_chr(refs, alts, formatClusterMutations, convert=TRUE))
  table <- dplyr::mutate(table, clusterType = purrr::map2_chr(plusStrand, minusStrand, getClusterType))
  table <- dplyr::mutate(table, distance = purrr::map_int(cMuts,function(x){x$distance[1]}))

  # Find the pattern intersect if needed --------------------------------------------------------
  if(patternIntersect){
    patternHeader <- dplyr::enquo(patternHeader)
    table <- dplyr::mutate(table, foundPatterns = purrr::map(cMuts,getPatternIntersect,!!patternHeader))
    table <- dplyr::mutate(table, has.intersect = purrr::map_lgl(foundPatterns, function(x){ifelse(length(x) > 0 && x[[1]] != "",TRUE,FALSE)}))
  }

  if(searchClusterPatterns){
    if(is.null(searchPatterns)){
      searchPatterns <- getSearchPatterns(searchReverseComplement)
      searchPatterns <- searchPatterns[nchar(searchPatterns$ref) > 1,]
    }
    table <- searchClusterPatterns(table,
                                    searchPatterns,
                                    searchRefHeader,
                                    searchAltHeader,
                                    searchDistanceHeader,
                                    searchIdHeader )
  }

  # Let know if no rows are found ---------------------------------------------------------------
  if(showWarning){
    if(nrow(table) == 0){
      warning("No rows found. Please make sure the cluster IDs are present and try again.")
    }
  }

  if(tibble){
    return(tibble::as.tibble(table))
  } else {
    return(as.data.frame(table))
  }


}

#' formatClusterMutations
#' @description Function to make a mutation description symbol (e.g. G>C).
#' @param refs A list containing the reference nucleotides.
#' @param alts A list containing the alternative nucleotides.
#' @param convert A boolean that tells if the nucleotides need to be converted to the reverse complement.
formatClusterMutations <- function(refs, alts, convert=FALSE) {

  # Check parameters -----------------------------------------------------------------
  stopifnot(length(refs) == length(alts))
  stopifnot(is.logical(convert))

  # Create the description symbols per mutation-------------------------------------

  if(convert){
    seq <- paste0(alts,".",refs)
    rev <- strsplit(as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq))),"\\.")[[1]]
    res <- paste(rev[1],rev[2],sep = ">")
  } else {
    res <- paste0(refs, ">", alts)
  }


  # Combine the mutation pairs -----------------------------------------------------
  finalResult <- dplyr::if_else(convert,
                         do.call(paste, as.list(rev(res))),
                         do.call(paste, as.list(res)))
  return(finalResult)
}

#' getClusterType
#' @description Function to sort and create a cluster type (e.g. C>G C>G / G>C G>C).
#' @param plusstrand A string with the plusStrand mutation description.
#' @param minusStrand A string with the minusStrand mutation description.
getClusterType <- function(plusStrand, minusStrand) {
  s <- sort(c(plusStrand, minusStrand))
  stringr::str_glue("{s[1]} / {s[2]}")
}

#' getPatternIntersect
#' @description A function to find the intersect of patterns between mutations.
#' @param clusterList A tibble with the mutations information.
#' @param patternHeader A string with the column header of the patterns.
getPatternIntersect <- function(clusterList,patternHeader){
  clusterList <- as.data.frame(clusterList)
  patternHeader <- rlang::get_expr(patternHeader)
  patterns <- c()
  counter <- 0
  for(index in 1:nrow(clusterList)){
   if(counter == 0){
     patterns <- clusterList[index,patternHeader][[1]]
   }
   patterns <- intersect(patterns, clusterList[index,patternHeader][[1]])
   counter <- counter+1
  }

  if(length(patterns) == 0){
    return(c(""))
  }
  return(patterns)
}

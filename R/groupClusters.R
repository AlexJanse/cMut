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
groupClusters <- function(table, clusterIdHeader = "clusterId",
                          refHeader = "ref", altHeader = "alt",
                          tibble = TRUE, patternIntersect = FALSE,
                          patternHeader = "linkedPatterns", showWarning = TRUE){

  # Check data --------------------------------------------------------------------------------------------
  stopifnot(nrow(table) > 0)
  stopifnot(!any(is.na(dplyr::select(table,clusterIdHeader,refHeader, altHeader))))

  # Build table -------------------------------------------------------------------------------------------
  table <- convertFactor(table)
  table <- data.table::as.data.table(table)
  table <- dplyr::group_by_(table, clusterIdHeader)
  table <- tidyr::nest(table, .key = "cMuts")
  table <- dplyr::filter(table, clusterId!="")
  table <- dplyr::mutate(table, refs = purrr::map(cMuts, ~as.character(dplyr::pull(., refHeader))),
           alts = purrr::map(cMuts, ~as.character(dplyr::pull(., altHeader))))
  table <- dplyr::mutate(table, plusStrand = purrr::map2_chr(refs, alts, formatClusterMutations),
           minusStrand = purrr::map2_chr(refs, alts, formatClusterMutations, convert=TRUE))
  table <- dplyr::mutate(table, clusterType = purrr::map2_chr(plusStrand, minusStrand, getClusterType))

  # Find the pattern intersect if needed --------------------------------------------------------
  if(patternIntersect){
    patternHeader <- dplyr::enquo(patternHeader)
    table <- dplyr::mutate(table, patternIntersect = purrr::map(cMuts,getPatternIntersect,!!patternHeader))
    table <- dplyr::mutate(table, has.intersect = purrr::map_lgl(patternIntersect, function(x){ifelse(length(x) > 0 && x[[1]] != "",TRUE,FALSE)}))
  }

  # Let know if no rows are found ---------------------------------------------------------------
  if(nrow(table) == 0){
    warning("No rows found. Please make sure the cluster IDs are present and try again.")
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
  combinedResult <-
    foreach::foreach(i = 1:length(refs),
            .combine = c) %do% {
              if(convert){
                res <- paste0(revNuc[refs[i]][[1]], ">", revNuc[alts[i]][[1]])
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
  patternHeader <- rlang::get_expr(patternHeader)
  patterns <- c()
  counter <- 0
  allPatterns <-
    foreach::foreach(index = 1:nrow(clusterList),
                     .combine = c) %do% {
                       if(counter == 0){
                         patterns <- clusterList[index,patternHeader][[1]][[1]]
                       }
                       patterns <- intersect(patterns, clusterList[index,patternHeader][[1]][[1]])
                       counter <- counter+1
                     }
  if(length(patterns) == 0){
    return(c(""))
  }
  return(patterns)
}

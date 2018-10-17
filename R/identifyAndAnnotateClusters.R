#' identifyAndAnnotateClusters
#' @description A function that finds and annotate clusters in a genomic data
#'   tibble.
#' @param x A tibble that contains at least chromosome nr., sampleID and
#'   position information. The data cannot contain any NA. For an example use
#'   createRandomMutations function.
#' @param maxDistance A number; The maximum distance between DNA mutations that
#'   count as clustered.
#' @param tibble A boolean if the table has to be a tibble
#' @param chromHeader A string; The name of the column with the chromosome nr.
#' @param sampleIdHeader A string; The name of the column with the sample ID.
#' @param positionHeader A string; The name of the column with the position nr.
#' @param linkPatterns A boolean to tell if patterns are needed to be found. If
#'   FALSE then the search... parameters are irrelevant.
#' @inheritParams linkPatterns
#' @inheritParams groupClusters
#' @param contextHeader A String; The name of the column with the context.
#' @return The tibble that was sent as an argument for this fuction with extra
#'   columns: clusterId, is.clustered and distance till nearest mutation below
#'   the maximum distance.
#' @export
#' @import magrittr
#' @examples
#' data <- testDataSet
#' results <- identifyAndAnnotateClusters(x = data,
#'                                        maxDistance = 20000)
#' resultsWithPatterns <- identifyAndAnnotateClusters(x = data,
#'                                                    maxDistance = 20000,
#'                                                    linkPatterns = TRUE)
identifyAndAnnotateClusters <- function(x, maxDistance, tibble = TRUE,
                                        chromHeader = "chrom", sampleIdHeader = "sampleIDs",
                                        positionHeader = "start", refHeader = "ref",
                                        altHeader = "alt", contextHeader = "surrounding",
                                        mutationSymbol = ".", linkPatterns = FALSE,
                                        reverseComplement = FALSE, searchPatterns = NULL,
                                        searchRefHeader = "ref", searchAltHeader = "alt", searchContextHeader = "surrounding",
                                        searchIdHeader = "process", searchReverseComplement = TRUE, patternsAsList = TRUE) {
  # Check if arguments are correct ------------------------------------------
  stopifnot(!any(is.na(dplyr::select(x,chromHeader,sampleIdHeader, positionHeader))))
  stopifnot(is.numeric(maxDistance))

  # Sort data ---------------------------------------------------------------
  x <- convertFactor(x)
  x <- data.table::as.data.table(x)
  x <- dplyr::arrange(x,
          dplyr::pull(x, chromHeader),
          dplyr::pull(x, sampleIdHeader),
          dplyr::pull(x, positionHeader))  # Used pulled function because we want the
                                           # string inside the variable and not the
                                           # variable name itself as column name.

  # Create GRange object ----------------------------------------------------
  ranges <- x %>%
    with(GenomicRanges::GRanges(seqnames = paste(
                   dplyr::pull(x, chromHeader),
                   dplyr::pull(x, sampleIdHeader)),
                 ranges = IRanges::IRanges(
                   dplyr::pull(x,positionHeader),
                   dplyr::pull(x,positionHeader))
                 )
         )


  # Add distance to nearest mutation information to the GRange object -------
  ranges <- addDistance(ranges,maxDistance)

  # Get clusterIDs ----------------------------------------------------------
  tempIds <- paste(dplyr::pull(x, chromHeader),
                   dplyr::pull(x, sampleIdHeader))
  clusterIds <- by(x,
                   factor(tempIds,
                     levels = unique(tempIds),
                          ordered = TRUE),
                   identifyClusters,
                   round(maxDistance), # Just to be sure that the max is a rounded number
                   positionHeader = positionHeader,
                   chromHeader = chromHeader,
                   sampleIdHeader = sampleIdHeader)

  # Add the cluster information to the x tibble ------------------------------
  ranges$clusterId <- unlist(clusterIds)
  x <- dplyr::arrange(x, dplyr::pull(x, chromHeader),
                  dplyr::pull(x, sampleIdHeader),
                  dplyr::pull(x, positionHeader))
  x <- dplyr::mutate(x,clusterId = ranges$clusterId,
                       is.clusteredTemp = ranges$is.clustered)
  x <- dplyr::mutate(x, is.clustered = purrr::map2_lgl(is.clusteredTemp, !!rlang::sym(refHeader), function(x,y){ifelse(x & y != "N",return(TRUE),return(FALSE))}),
                        distance = ranges$distance)

  x$is.clusteredTemp <- NULL

  # Add information about if the mutation can be linked to a certain pattern ------
  if(linkPatterns){
    if(is.null(searchPatterns)){
      searchPatterns <- getSearchPatterns(searchReverseComplement)
      searchReverseComplement <- FALSE
    }
    x <- addLinkPatterns(x, refHeader, altHeader, contextHeader,
                         mutationSymbol, reverseComplement,
                         searchPatterns, searchRefHeader,
                         searchAltHeader,  searchContextHeader,
                         searchIdHeader, searchReverseComplement)
  }
  if(tibble){
    x <- tibble::as.tibble(x)
  } else {
    x <- as.data.frame(x)
  }
  return(x)

}

#-----------------------------------------------------------------------------------------------
#' addDistance
#' @description A function that adds distances to nearest mutation to the tibble of the identifyAndAnnotateClusters function
#' @inheritParams identifyAndAnnotateClusters
#' @param ranges A GRange object which were created during identifyAndAnnotateClusters function
#' @param maxDistance A number with the maximum distance
#' @return A GRange with added distance and a logical column as metadata
addDistance <- function(ranges, maxDistance) {

  hits <- GenomicRanges::distanceToNearest(ranges)
  ranges$distance <- NA
  ranges$distance[S4Vectors::queryHits(hits)] <- S4Vectors::elementMetadata(hits)$distance

  ranges$is.clustered <- (ranges$distance <= maxDistance & ! is.na(ranges$distance))

  return(ranges)
}

#' callLinkPatterns
#' @description A function to correctly call the linkPattern function
callLinkPatterns <- function(x,linkedVariables){
  # linkedVariables <- rlang::get_expr(linkedVariables)
  mutation <- strsplit(x,"!")
  return(linkPatterns(ref = mutation[[1]][1], alt = mutation[[1]][2],
                      context = mutation[[1]][3], distance = as.numeric(mutation[[1]][4]),
                      mutationSymbol = linkedVariables[[1]], reverseComplement = linkedVariables[[2]],
                      searchPatterns = linkedVariables[[3]], searchRefHeader = linkedVariables[[4]],
                      searchAltHeader = linkedVariables[[5]], searchContextHeader = linkedVariables[[6]],
                      searchIdHeader = linkedVariables[[7]], searchReverseComplement = linkedVariables[[8]]))
}

#' addLinkPatterns
#' @inheritParams identifyAndAnnotateClusters
#' @export
#' @import foreach
addLinkPatterns <- function(x, refHeader = "ref",
                            altHeader = "alt",
                            contextHeader = "surrounding",
                            mutationSymbol = ".",
                            reverseComplement = FALSE,
                            searchPatterns = NULL,
                            searchRefHeader = "ref",
                            searchAltHeader = "alt",
                            searchContextHeader = "surrounding",
                            searchIdHeader = "process",
                            searchReverseComplement = TRUE,
                            checkHeader = "is.clustered"){

  linkVariables <- list(mutationSymbol, reverseComplement,
                        searchPatterns, searchRefHeader,
                        searchAltHeader,  searchContextHeader,
                        searchIdHeader, searchReverseComplement) # Variables are put in a list to reduce the amount of code

  x <- dplyr::mutate(x, tempMutColumn = paste(!!rlang::sym(refHeader),
                                        !!rlang::sym(altHeader),
                                        !!rlang::sym(contextHeader),
                                        distance,
                                        sep = "!"))
  x <- dplyr::mutate(x, linkedPatterns = purrr::map2(tempMutColumn,
                                               !!rlang::sym(checkHeader),
                                               function(x,y){
                                                 ifelse(y,callLinkPatterns(x,linkVariables),list(""))
                                               }))
  x <- dplyr::mutate(x, linkedPatterns = purrr::map(linkedPatterns,function(x){x[[1]]}))
  x <- dplyr::mutate(x, is.linked = purrr::map_lgl(linkedPatterns,
                                             function(x){
                                               dplyr::if_else(x[[1]][[1]] != "" && x[[1]][[1]] != "NA",
                                                              TRUE, FALSE)
                                             }))


  x$tempMutColumn <- NULL
  return(x)
}

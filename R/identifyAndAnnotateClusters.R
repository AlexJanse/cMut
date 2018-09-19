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
#' @param linkPatterns A boolean to tell if patterns are needed to be found.
#' @inheritParams linkPatterns
#' @inheritParams groupClusters
#' @param contextHeader A String; The name of the column with the context.
#' @return The tibble that was sent as an argument for this fuction with extra
#'   columns: clusterId, is.clustered and distance till nearest mutation below
#'   the maximum distance.
#' @export
#' @import magrittr
#' @examples
#' data <- createRandomMutations(1000)
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
                                        searchIdHeader = "proces", searchReverseComplement = TRUE) {
  # Check if arguments are correct ------------------------------------------
  stopifnot(!any(is.na(dplyr::select(x,chromHeader,sampleIdHeader, positionHeader))))
  stopifnot(is.numeric(maxDistance))

  # Sort data ---------------------------------------------------------------
  x <- dplyr::arrange(x,
          dplyr::pull(x, chromHeader),
          dplyr::pull(x, sampleIdHeader),
          dplyr::pull(x,positionHeader))  # Used pulled function because we want the
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
  clusterIds <- by(x,
                   factor(
                     paste(dplyr::pull(x, chromHeader),
                           dplyr::pull(x, sampleIdHeader)),
                     levels = unique(
                                paste(
                                  dplyr::pull(x, chromHeader),
                                  dplyr::pull(x, sampleIdHeader))),
                                ordered = TRUE),
                   identifyClusters,
                   round(maxDistance), # Just to be sure that the max is a rounded number
                   positionHeader = positionHeader,
                   chromHeader = chromHeader,
                   sampleIdHeader = sampleIdHeader)

  # Add the cluster information to the x tibble ------------------------------
  ranges$clusterId <- unlist(clusterIds)
  x <- x %>%
    dplyr::arrange(dplyr::pull(x, chromHeader),
                  dplyr::pull(x, sampleIdHeader),
                  dplyr::pull(x, positionHeader)) %>%
    dplyr::mutate(clusterId = ranges$clusterId,
           is.clusteredTemp = ranges$is.clustered) %>%
    dplyr::mutate(is.clustered = purrr::map2_lgl(is.clusteredTemp, !!rlang::sym(refHeader), function(x,y){ifelse(x & y != "N",return(TRUE),return(FALSE))})) %>%
    dplyr::mutate(distance = ranges$distance)

  x$is.clusteredTemp <- NULL

  # Add information about if the mutation can be linked to a certain pattern ------
  if(linkPatterns){
    linkVariables <- list(mutationSymbol, reverseComplement,
                       searchPatterns, searchRefHeader,
                       searchAltHeader,  searchContextHeader,
                       searchIdHeader, searchReverseComplement) # Variables are put in a list to reduce the amount of code

    x <- x %>%
      dplyr::mutate(tempMutColumn = paste(!!rlang::sym(refHeader),
                                          !!rlang::sym(altHeader),
                                          !!rlang::sym(contextHeader),
                                          sep = "!")) %>%
      dplyr::mutate(linkedPatterns = purrr::map2(tempMutColumn,
                                                 is.clustered,
                                                 function(x,y){
                                                   ifelse(y,callLinkPatterns(x,linkVariables),list(""))
                                                   })) %>%
      dplyr::mutate(is.linked = purrr::map_lgl(linkedPatterns,
                                               function(x){
                                                 dplyr::if_else(x[[1]] != "" && x[[1]] != "NA",
                                                                TRUE, FALSE)
                                                 }))


    x$tempMutColumn <- NULL

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
  return(linkPatterns(mutation[[1]][1],mutation[[1]][2],mutation[[1]][3], linkedVariables[[1]],
                      linkedVariables[[2]], linkedVariables[[3]], linkedVariables[[4]],
                      linkedVariables[[5]], linkedVariables[[6]],
                      linkedVariables[[7]], linkedVariables[[8]]))
}

#' identifyAndAnnotateClusters
#' @description A function that finds and annotate clusters in a genomic data tibble.
#' @param x A tibble that contains at least chromosome nr., sampleID and position information.
#' The data cannot contain any NA.
#' @param maxDistance A number; The maximum distance between DNA mutations that count as clustered.
#' @param chromHeader A string; The name of the column with the chromosome nr.
#' @param sampleIdHeader A string; The name of the column with the sample ID.
#' @param positionHeader A string; The name of the column with the position nr.
#' @inheritParams linkPatterns
#' @inheritParams groupClusters
#' @param contextHeader A String; The name of the column with the context.
#' @return The tibble that was sent as an argument for this fuction with extra columns:
#' clusterId, is.clustered and distance till nearest mutation below the maximum distance.
#' @export
#' @import magrittr
identifyAndAnnotateClusters <- function(x, maxDistance,
                                        chromHeader = "chrom", sampleIdHeader = "sampleIDs",
                                        positionHeader = "start", linkPatterns = FALSE,
                                        reverseComplement = FALSE, searchPatterns = NULL,
                                        searchRefHeader = "ref", searchAltHeader = "alt", searchContextHeader = "surrounding",
                                        searchIdHeader = "proces", searchReverseComplement = TRUE,
                                        refHeader = "ref", altHeader = "alt", contextHeader = "surrounding", mutationSymbol = ".") {

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
         ) # TODO Make the function able to handle CNVs


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
           is.clustered = ranges$is.clustered,
           distance = ranges$distance)

  if(linkPatterns){
    reverseComplement <- dplyr::enquo(reverseComplement)
    searchPatterns <- dplyr::enquo(searchPatterns)
    searchRefHeader <- dplyr::enquo(searchRefHeader)
    searchAltHeader <- dplyr::enquo(searchAltHeader)
    searchContextHeader <- dplyr::enquo(searchContextHeader)
    searchIdHeader <- dplyr::enquo(searchIdHeader)
    searchReverseComplement <- dplyr::enquo(searchReverseComplement)
    mutationSymbol <- dplyr::enquo(mutationSymbol)


    x <- x %>%
      dplyr::mutate(tempMutColumn = paste(!!rlang::sym(refHeader),
                                          !!rlang::sym(altHeader),
                                          !!rlang::sym(contextHeader),
                                          sep = "$")) %>%
      dplyr::mutate(linkedPatterns = purrr::map(tempMutColumn,
                                                callLinkPatterns,reverseComplement,
                                                searchPatterns, searchRefHeader,
                                                searchAltHeader,  searchContextHeader,
                                                searchIdHeader, searchReverseComplement,
                                                mutationSymbol))
    x$tempMutColumn <- NULL
    return(x)
    } else {
    return(x)
  }
}

#-----------------------------------------------------------------------------------------------
#' addDistance
#' @description A function that adds distances to nearest mutation to the tibble of the identifyAndAnnotateClusters function
#' @inheritParams identifyAndAnnotateClusters
#' @param ranges A GRange object which were created during identifyAndAnnotateClusters function
#' @return A GRange with added distance and a logical column as metadata
addDistance <- function(ranges, maxDistance) {

  hits <- GenomicRanges::distanceToNearest(ranges)
  ranges$distance <- NA
  ranges$distance[S4Vectors::queryHits(hits)] <- S4Vectors::elementMetadata(hits)$distance

  ranges$is.clustered <- (ranges$distance <= maxDistance & ! is.na(ranges$distance))

  return(ranges)
}

callLinkPatterns <- function(x,reverseComplement,
                             searchPatterns, searchRefHeader,
                             searchAltHeader,  searchContextHeader,
                             searchIdHeader, searchReverseComplement,
                             mutationSymbol){

  reverseComplement <- dplyr::get_expr(reverseComplement)
  searchPatterns <- dplyr::get_expr(searchPatterns)
  searchRefHeader <- dplyr::get_expr(searchRefHeader)
  searchAltHeader <- dplyr::get_expr(searchAltHeader)
  searchContextHeader <- dplyr::get_expr(searchContextHeader)
  searchIdHeader <- dplyr::get_expr(searchIdHeader)
  searchReverseComplement <- dplyr::get_expr(searchReverseComplement)
  mutationSymbol <- dplyr::get_expr(mutationSymbol)
  mutation <- strsplit(x,"$")

  return(linkPatterns(mutation[1],mutation[2],mutation[3], mutationSymbol,
                      reverseComplement, searchPatterns, searchRefHeader,
                      searchAltHeader,  searchContextHeader,
                      searchIdHeader, searchReverseComplement))
}

#' identifyAndAnnotateClusters
#' @description A function that finds and annotate clusters in a genomic data tibble.
#' @param x A tibble that contains at least chromosome nr., sampleID and position information.
#' The data cannot contain any NA.
#' @param maxDistance A number; The maximum distance between DNA mutations that count as clustered.
#' @param chromHeader A string; The name of the column with the chromosome nr.
#' @param sampleIdHeader A string; The name of the column with the sample ID.
#' @param positionHeader A string; The name of the column with the position nr.
#' @return The tibble that was sent as an argument for this fuction with extra columns:
#' clusterId, is.clustered and distance till nearest mutation below the maximum distance.
#' @export
#' @import magrittr
identifyAndAnnotateClusters <- function(x,
                                        maxDistance,
                                        chromHeader = "Chr",
                                        sampleIdHeader = "sampleID",
                                        positionHeader = "Pos") {

  # Check if arguments are correct ------------------------------------------
  stopifnot(!any(is.na(dplyr::select(x,chromHeader,sampleIdHeader, positionHeader))))
  stopifnot(is.numeric(maxDistance))

  # Sort data ---------------------------------------------------------------
  x <- dplyr::arrange(x,
          dplyr::pull(x, chromHeader),
          dplyr::pull(x, sampleIdHeader),
          dplyr::pull(x,positionHeader)) # Used pulled function because we want the
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
  x %>%
    dplyr::arrange(dplyr::pull(x, chromHeader),
                  dplyr::pull(x, sampleIdHeader),
                  dplyr::pull(x,positionHeader)) %>%
    dplyr::mutate(clusterId = ranges$clusterId,
           is.clustered = ranges$is.clustered,
           distance = ranges$distance) %>%
    return()
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

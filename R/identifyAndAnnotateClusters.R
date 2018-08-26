#' identifyAndAnnotateClusters
#' @description A function that finds and annotate clusters in a genomic data tibble.
#' @param x A tibble that contains at least chromosome nr., sampleID and position information.
#' The data cannot contain any NA.
#' @param maxDistance A integer; The maximum distance between DNMs that count as clustered.
#' @param chromHeader A string; The name of the column with the chromosome nr.
#' @param sampleIDHeader A string; The name of the column with the sample ID.
#' @param positionHeader A string; The name of the column with the position nr.
#' @return The tibble that was sent as an argument for this fuction with extra columns:
#' clusterId, is.clustered and distance till nearest mutation below the maximum distance.
#' @export
identifyAndAnnotateClusters <- function(x,
                                        maxDistance,
                                        chromHeader = "Chr",
                                        sampleIdHeader = "sampleID",
                                        positionHeader = "Pos") {
  # Sort data ---------------------------------------------------------------

  x <- arrange(x,
          pull(x, chromHeader),
          pull(x, sampleIdHeader),
          pull(x,positionHeader)) # Used pulled function because we want the
                                  # string inside the variable and not the
                                  # variable name it self as column name.

  # Create GRange object ----------------------------------------------------
  ranges <- x %>%
    with(GRanges(seqnames = paste(
                   pull(x, chromHeader),
                   pull(x, sampleIdHeader)),
                 ranges = IRanges(
                   pull(x,positionHeader),
                   pull(x,positionHeader))
                 )
         ) # TODO Make the function able to handle CNVs


  # Add distance to nearest mutation information to the GRange object -------
  ranges <- addDistance(ranges,maxDistance)

  # Get clusterIDs ----------------------------------------------------------
  clusterIds <- by(x,
                   factor(paste(x$Chr, x$sampleID),
                          levels = unique(paste(x$Chr, x$sampleID)),
                          ordered = TRUE),
                   identifyClusters,
                   maxDistance,
                   positionHeader="Pos",
                   chromHeader="Chr",
                   sampleIdHeader="sampleID")

  # Add the cluster information to the x tibble ------------------------------
  ranges$clusterId <- unlist(clusterIds)
  x %>%
    arrange(pull(x, chromHeader),
            pull(x, sampleIdHeader),
            pull(x,positionHeader)) %>%
    mutate(clusterId = ranges$clusterId,
           is.clustered = ranges$is.clustered,
           distance = ranges$distance) %>%
    return()
}

#-----------------------------------------------------------------------------------------------
#' addDistance
#' @description A function that adds distances to nearest DNM to the tibble of the identifyAndAnnotateClusters function
#' @inheritParams identifyAndAnnotateClusters
#' @param ranges A GRange object which were created during identifyAndAnnotateClusters function
#' @return A GRange with added distance and a logical column as metadata
addDistance <- function(ranges, maxDistance) {

  hits <- distanceToNearest(ranges)
  ranges$distance <- NA
  ranges$distance[queryHits(hits)] <- elementMetadata(hits)$distance

  ranges$is.clustered <- (ranges$distance <= maxDistance & ! is.na(ranges$distance))

  return(ranges)
}

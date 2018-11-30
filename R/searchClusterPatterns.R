#' searchClusterPatterns
#' @description A function to search cluster patterns based on reference,
#'   variant and distance of the nucleotides within the cluster.
#' @param groupedClusters Table gained by using the \code{\link{groupClusters}}
#'   function.
#' @inheritParams groupClusters
#' @import magrittr
searchClusterPatterns <- function(groupedClusters,
                                  searchPatterns       = mutationPatterns,
                                  searchRefHeader      = "ref",
                                  searchAltHeader      = "alt",
                                  searchDistanceHeader = "maxDistance",
                                  searchIdHeader       = "process") {

  # Check parameters -----------------------------------------------
  stopifnot(all(nchar(searchPatterns[ ,searchRefHeader]) >  1 |
                nchar(searchPatterns[ ,searchRefHeader]) == 0))

  stopifnot(all(nchar(searchPatterns[ ,searchAltHeader]) >  1))

  stopifnot(!any(is.na(dplyr::select(searchPatterns,  searchRefHeader,
                                     searchAltHeader, searchDistanceHeader))))

  stopifnot(!any(is.null(dplyr::select(searchPatterns,  searchRefHeader,
                                       searchAltHeader, searchDistanceHeader))))

  stopifnot("cMuts" %in% names(groupedClusters))

  # Convert search table to the correct class --------------------
  searchPatterns <- as.data.frame(searchPatterns)


  # Try to find cluster patterns in the table --------------------

  # Create list to fill in the results:
  clusterPatterns <- list()

  # Loop over the sent data:
  for (index in 1:nrow(groupedClusters)) {
    row <- groupedClusters[index, ]
    subclusterPatterns <- c()

    for (pattIndex in 1:nrow(searchPatterns)) {

      # Collect data from the search table:
      id          <- searchPatterns[pattIndex, searchIdHeader]
      refPat      <- searchPatterns[pattIndex, searchRefHeader]
      altPat      <- searchPatterns[pattIndex, searchAltHeader]
      maxDistance <- searchPatterns[pattIndex, searchDistanceHeader]
      if (is.na(refPat)) {
        refPat <- ""
      }

      # Check if the data from the search table match with the data
      #   in the sent data row:
      if ((refPat == row$refs | refPat == "") & altPat == row$alts) {
        clusterDistance <- max(row$distance[[1]])

        if (clusterDistance <= maxDistance) {
          subclusterPatterns[length(subclusterPatterns) + 1] <- id
        }
      }
    }

    # Adjust the subresults based on the content:
    if (is.null(subclusterPatterns)) {
      clusterPatterns[[length(clusterPatterns) + 1]] <- c("")
    } else {
      clusterPatterns[[length(clusterPatterns) + 1]] <- unique(subclusterPatterns)
    }

  }

  # List for the column that tells if results were found for that row:
  checkList <- clusterPatterns != ""

  # Add results to table -----------------------------------------------

  # Add temporary column with the cluster patterns results:
  groupedClusters <- dplyr::mutate(groupedClusters,
                                   clusterPatterns = clusterPatterns)

  # Fuse the found mutation patterns with the found cluster patterns:
  groupedClusters <- dplyr::mutate(groupedClusters,
                                   foundPatterns = purrr::map2(foundPatterns,
                                                               clusterPatterns,
                                                               function(x,y) {
                                                                 ifelse(x == "",
                                                                        list(c(y)),
                                                                        ifelse(y == "",
                                                                               list(c(x)),
                                                                               list(c(x, y))))
                                                                  }))

  # Unlist foundPatterns column:
  groupedClusters <- dplyr::mutate(groupedClusters,
                                   foundPatterns = purrr::map(foundPatterns,
                                                              function(x) {
                                                                x[[1]]
                                                              }))

  # Add check column:
  groupedClusters <- dplyr::mutate(groupedClusters,
                                   has.clusterPatterns = checkList)

  # Remove temporary column and return results:
  groupedClusters$clusterPatterns <- NULL
  return(groupedClusters)
}



#' identifyAndAnnotateClusters
#' @description A function that finds and annotate clusters in a genomic data
#'   tibble.
#' @param x A tibble that contains at least chromosome name, sample ID and
#'   position information. The data cannot contain any NA. For an example use
#'   \code{\link{createRandomMutations}} function.
#' @param maxDistance A number; The maximum distance between DNA mutations that
#'   count as clustered.
#' @param asTibble A boolean if the result table has to be a tibble
#' @param chromHeader A string with the name of the column with the chromosome
#'   name. (So the data in the column needs to be notated as e.g. "chr2")
#' @param sampleIdHeader A string with the name of the column with the sample
#'   ID.
#' @param positionHeader A string with the name of the column with the position
#'   of the mutation. The data in the column needs to be numberic.
#' @param linkPatterns A Boolean to tell if it's necessary to try and link the
#'   mutations to patterns. If FALSE then the search... parameters are
#'   irrelevant.
#' @inheritParams linkPatterns
#' @inheritParams groupClusters
#' @param contextHeader A string with the name of the column with the context.
#'   The data inside this column is e.g. "C.G" hereby stands the "." for the
#'   location of the mutation. What symbol is used to describe this location is
#'   irrelevant but be sure to adjust the \code{mutationSymbol} accordingly when
#'   searching for patterns. The \code{contextHeader} is irrelevant if
#'   \code{linkPatterns} is FALSE.
#' @return The tibble that was sent as an argument for this fuction with extra
#'   columns: clusterId, is.clustered and distance till nearest mutation below
#'   the maximum distance.
#' @export
#' @import magrittr
#' @examples
#' # Example data set:
#' data <- testDataSet
#'
#' # Example for just clustering:
#' results <- identifyAndAnnotateClusters(x = data,
#'                                        maxDistance = 20000)
#'
#' # Example for clustering and linking patterns with the default searchPattern table:
#' results <- identifyAndAnnotateClusters(x = data,
#'                                                    maxDistance = 20000,
#'                                                    linkPatterns = TRUE)
#' # See the getSearchPattern funtion for more information about it.
#'
#' # For more information about the added columns, use:
#' cat(comment(results))
#' @seealso
#' \itemize{
#'   \item \code{\link{createRandomMutations}} for an example
#'   of data as input for parameter \code{x}
#'   \item \code{\link{mutationPatterns}} for looking at the default pattern search
#'   table
#'   }
identifyAndAnnotateClusters <- function(x, maxDistance, asTibble = TRUE,
                                        chromHeader = "chrom", sampleIdHeader = "sampleIDs",
                                        positionHeader = "start", refHeader = "ref",
                                        altHeader = "alt", contextHeader = "surrounding",
                                        mutationSymbol = ".", linkPatterns = FALSE,
                                        reverseComplement = FALSE, searchPatterns = NULL,
                                        searchRefHeader = "ref", searchAltHeader = "alt", searchContextHeader = "surrounding",
                                        searchIdHeader = "process", searchDistanceHeader = "maxDistance",
                                        searchReverseComplement = TRUE, patternsAsList = TRUE) {
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
    x <- addLinkPatterns(x, refHeader = refHeader,
                         altHeader = altHeader,
                         contextHeader = contextHeader,
                         mutationSymbol = mutationSymbol,
                         reverseComplement = reverseComplement,
                         searchPatterns = searchPatterns,
                         searchRefHeader = searchRefHeader,
                         searchAltHeader = searchAltHeader,
                         searchContextHeader = searchContextHeader,
                         searchIdHeader = searchIdHeader,
                         searchReverseComplement = searchReverseComplement,
                         searchDistanceHeader = searchDistanceHeader)
  }
  if(asTibble){
    x <- tibble::as.tibble(x)
  } else {
    x <- as.data.frame(x)
  }

  comment(x) <-
  "Information about the added columns:
  distance        : Column with the distances to
                    the nearest mutation.
  clusterId       : Column with the ID of a column.
                    It consist of the chromosome name,
                    sampleID and the sample unique
                    cluster ID number. This is al
                    seperated with a space.
  is.clustered    : Column with Boolean if the mutation
                    is part of a cluster.

  (When using the linkPatterns parameter with TRUE):
  linkedPatterns  : Column with the names from the
                    searchPatterns table that matched
                    with the mutation. The names are
                    put in a vector.
  is.linked       : Column with Boolean if there are
                    found patterns for the mutation."
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
                      searchIdHeader = linkedVariables[[7]], searchReverseComplement = linkedVariables[[8]],
                      searchDistanceHeader = linkedVariables[[9]]))
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
                            checkHeader = "is.clustered",
                            searchDistanceHeader = "maxDistance"){

  linkVariables <- list(mutationSymbol, reverseComplement,
                        searchPatterns, searchRefHeader,
                        searchAltHeader,  searchContextHeader,
                        searchIdHeader, searchReverseComplement,
                        searchDistanceHeader) # Variables are put in a list to reduce the amount of code

  x <- dplyr::mutate(x, tempMutColumn = paste(!!rlang::sym(refHeader),
                                        !!rlang::sym(altHeader),
                                        !!rlang::sym(contextHeader),
                                        distance,
                                        sep = "!"))
  x <- dplyr::mutate(x, linkedPatterns = purrr::map2(tempMutColumn,
                                               !!rlang::sym(checkHeader),
                                               function(x,y){
                                                 ifelse(y,callLinkPatterns(x,linkVariables),"")
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

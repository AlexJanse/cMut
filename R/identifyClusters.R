#' identifyClusters
#' @description A function that finds and annotate clusters in a genomic data
#'   tibble.
#' @param dataTable A data.frame or tibble that contains at least chromosome
#'   name, sample ID and position information. The data cannot contain any NA.
#'   For an example use \code{\link{testDataSet}}.
#' @param maxDistance A number with the maximum distance between DNA mutations
#'   that are defined as being in a cluster.
#' @param chromHeader A string with the name of the column with the chromosome
#'   name. (So the data in the column needs to be notated as e.g. "chr2")
#' @param sampleIdHeader A string with the name of the column with the sample
#'   ID.
#' @param positionHeader A string with the name of the column with the position
#'   of the mutation. (The data in the column needs to be nummeric.)
#' @param linkPatterns A Boolean to tell if it's necessary to try and link the
#'   mutations to patterns. If FALSE then the \code{search...} parameters and
#'   the \code{linkClustersOnly} parameter are irrelevant. For more information
#'   see \code{\link{linkPatterns}}.
#' @inheritParams linkPatterns
#' @inheritParams groupClusters
#' @param contextHeader A string with the name of the column with the context.
#'   The data inside this column is e.g. "C.G" hereby stands the "." for the
#'   location of the mutation. What symbol is used to describe this location is
#'   arbitrary but be sure to adjust the \code{mutationSymbol} accordingly when
#'   searching for patterns. The \code{contextHeader} is irrelevant if
#'   \code{linkPatterns} is FALSE.
#' @param linkClustersOnly A boolean to tell if only the clustered mutations are
#'   needed to be linked with the patterns in the \code{searchPatterns} table.
#'   When it is FALSE all the mutations will be used.
#' @param asTibble A boolean to tell if the result table has to be a tibble.
#'   When it is FALSE it will return data.frame
#' @return The tibble that was sent as an argument for this function with extra
#'   columns: clusterId, is.clustered and distance till nearest mutation below
#'   the maximum distance.
#' @export
#' @import magrittr
#' @examples
#' # Example data set:
#' data <- testDataSet
#'
#' # Example for just clustering:
#' results <- identifyClusters(dataTable   = data,
#'                             maxDistance = 20000,
#'                             linkPatterns = FALSE)
#'
#' # Example for clustering and linking patterns with the default searchPattern table:
#' results <- identifyClusters(dataTable    = data,
#'                             maxDistance  = 20000,
#'                             linkPatterns = TRUE)
#'
#' # For more information about the added columns, use:
#' cat(comment(results))
#'
#' @seealso \itemize{ \item \code{\link{testDataSet}} for an example of data as
#'   input for parameter \code{dataTable} \item \code{\link{mutationPatterns}}
#'   for looking at the default pattern search table \item Use the following
#'   code to access the vignette with detailed examples of how to use the
#'   functions of cMut: vignette("analysis_of_clusterpattterns",package =
#'   "cMut") }
identifyClusters <- function(dataTable,                               maxDistance,
                             chromHeader             = "chrom",       sampleIdHeader       = "sampleIDs",
                             positionHeader          = "start",       refHeader            = "ref",
                             altHeader               = "alt",         contextHeader        = "surrounding",
                             mutationSymbol          = ".",           linkPatterns         = TRUE,
                             reverseComplement       = FALSE,         searchPatterns       = NULL,
                             searchRefHeader         = "ref",         searchAltHeader      = "alt",
                             searchContextHeader     = "surrounding", searchIdHeader       = "process",
                             searchDistanceHeader    = "maxDistance", searchMutationSymbol = ".",
                             searchReverseComplement = TRUE,          linkClustersOnly     = TRUE,
                             renameReverse           = FALSE,         asTibble             = TRUE) {

  # Check if arguments are correct ------------------------------------------
  if (length(setdiff(c(chromHeader,sampleIdHeader, positionHeader),
                      names(dataTable))) > 0) {
    stop ("Error: Please check if the header parameters match with
          the column names of the sent table.")
  }

  stopifnot(is.numeric(maxDistance))
  stopifnot(is.logical(c(reverseComplement,       renameReverse,
                         searchReverseComplement, asTibble,
                         linkClustersOnly,        linkPatterns)))

  dataTable <- createIdentTable(dataTable      = dataTable,      maxDistance    = maxDistance,
                                chromHeader    = chromHeader,    sampleIdHeader = sampleIdHeader,
                                positionHeader = positionHeader, refHeader      = refHeader,
                                altHeader      = altHeader,      contextHeader  = contextHeader)

  # Search for patterns if asked --------------------------------------------
  if (linkPatterns) {
    dataTable <- addLinkPatterns(table                = dataTable,            refHeader               = refHeader,
                                 altHeader            = altHeader,            contextHeader           = contextHeader,
                                 mutationSymbol       = mutationSymbol,       reverseComplement       = reverseComplement,
                                 searchPatterns       = searchPatterns,       searchRefHeader         = searchRefHeader,
                                 searchAltHeader      = searchAltHeader,      searchContextHeader     = searchContextHeader,
                                 searchIdHeader       = searchIdHeader,       searchMutationSymbol    = searchMutationSymbol,
                                 searchDistanceHeader = searchDistanceHeader, searchReverseComplement = searchReverseComplement,
                                 linkClustersOnly     = linkClustersOnly,     renameReverse           = renameReverse)
  }


  # Add comments to explain the added columns -------------------------------
  dataTable <- addIdentTableComment(dataTable, linkPatterns)


  # Return the results in the prefered class --------------------------------
  if(asTibble){
    return(tibble::as.tibble(dataTable))
  } else {
    return(as.data.frame(dataTable))
  }

}


#' addDistance
#' @description A function that adds distances to nearest mutation to the tibble
#'   of the identifyClusters function
#' @inheritParams identifyClusters
#' @param ranges A GRange object which were created during
#'   identifyClusters function
#' @return A GRange object with added distance and a logical column as metadata
addDistance <- function(ranges, maxDistance) {

  hits <- GenomicRanges::distanceToNearest(ranges)
  ranges$distance <- NA
  ranges$distance[S4Vectors::queryHits(hits)] <- S4Vectors::elementMetadata(hits)$distance

  ranges$is.clustered <- (ranges$distance <= maxDistance & !is.na(ranges$distance))

  return(ranges)
}


#' addLinkPatterns
#' @description A function to add columns to the sent table with information
#'   gathered from the \code{\link{linkPatterns}} function.
#' @inheritParams identifyClusters
#' @param table table with the clusters information
#' @param checkHeader name of the Boolean check column
#' @import foreach
#' @importFrom rlang .data
addLinkPatterns <- function(table,                refHeader,
                            altHeader,            contextHeader,
                            mutationSymbol,       reverseComplement,
                            searchPatterns,       searchRefHeader,
                            searchAltHeader,      searchContextHeader,
                            searchIdHeader,       searchReverseComplement,
                            searchMutationSymbol, checkHeader,
                            searchDistanceHeader, linkClustersOnly,
                            renameReverse) {

  # Get deafault table if nothing is sent -----------------------------------
  # We do this before calling linkPatterns to save time:
  if (is.null(searchPatterns)) {
    searchPatterns <- getSearchPatterns(reverse       = searchReverseComplement,
                                        renameReverse = renameReverse)
  } else {
    searchPatterns <- convertFactor(searchPatterns)

    if (searchReverseComplement) {
      searchPatterns <- dplyr::bind_rows(searchPatterns,
                                         getRevComTable(table         = searchPatterns,  refHeader     = searchRefHeader,
                                                        altHeader     = searchAltHeader, contextHeader = searchContextHeader,
                                                        renameReverse = renameReverse,   idHeader      = searchIdHeader))
    }
  }
  # set the reverse parameters on FALSE beause it has already been processed above:
  searchReverseComplement <- FALSE
  renameReverse           <- FALSE


  # Combine variables for easier transport ----------------------------------
  linkVariables <- list(mutationSymbol,       reverseComplement,
                        searchPatterns,       searchRefHeader,
                        searchAltHeader,      searchContextHeader,
                        searchIdHeader,       searchReverseComplement,
                        searchDistanceHeader, searchMutationSymbol,
                        renameReverse)


  # Build table -------------------------------------------------------------
  # Create a temporary column with all the necessary information for the
  #   linkPattern function:
  table <- dplyr::mutate(table,
                         tempMutColumn = paste(!!rlang::sym(refHeader),
                                               !!rlang::sym(altHeader),
                                               !!rlang::sym(contextHeader),
                                               .data$distance,
                                               sep = "!"))

  # Call the callLinkPatterns function:
  table <- dplyr::mutate(table,
                         linkedPatterns = purrr::map2(.data$tempMutColumn,
                                                      .data$is.clustered,
                                               function(x, y) {
                                                 ifelse(linkClustersOnly,
                                                        ifelse(y,
                                                          callLinkPatterns(x,
                                                                           linkVariables),
                                                          ""),
                                                        callLinkPatterns(x,
                                                                         linkVariables))
                                               }))

  # Unlist the linked patterns:
  table <- dplyr::mutate(table,
                         linkedPatterns = purrr::map(.data$linkedPatterns,
                                                     function(x) {
                                                       x[[1]]
                                                     }))

  # Add a Boolean column to tell if the row contains linked patterns:
  table <- dplyr::mutate(table,
                         is.linked = purrr::map_lgl(.data$linkedPatterns,
                                                    function(x){
                                                      x[[1]][[1]] != "" &
                                                      x[[1]][[1]] != "NA"
                                                    }))


  # remove the temporary column and return the table:
  table$tempMutColumn <- NULL
  return(table)
}


#' callLinkPatterns
#' @description A function to correctly call the linkPattern function
#' @param x contains mutation information
#' @param linkedVariables Contains all the parameters needed voor
#'   \code{\link{linkPatterns}}
callLinkPatterns <- function(x, linkedVariables) {
  # linkedVariables <- rlang::get_expr(linkedVariables)
  mutation <- strsplit(x,"!")
  return(linkPatterns(ref                  = mutation[[1]][1],     alt                     = mutation[[1]][2],
                      context              = mutation[[1]][3],     distance                = as.numeric(mutation[[1]][4]),
                      mutationSymbol       = linkedVariables[[1]], reverseComplement       = linkedVariables[[2]],
                      searchPatterns       = linkedVariables[[3]], searchRefHeader         = linkedVariables[[4]],
                      searchAltHeader      = linkedVariables[[5]], searchContextHeader     = linkedVariables[[6]],
                      searchIdHeader       = linkedVariables[[7]], searchReverseComplement = linkedVariables[[8]],
                      searchDistanceHeader = linkedVariables[[9]], searchMutationSymbol    = linkedVariables[[10]],
                      renameReverse        = linkedVariables[[11]]))
}


#' createIdentTable
#' @description Creates the table with extra annotation about the present
#'   clusters
#' @inheritParams identifyClusters
#' @inheritParams linkPatterns
createIdentTable <- function(dataTable,      maxDistance,
                             chromHeader,    sampleIdHeader,
                             positionHeader, refHeader,
                             altHeader,      contextHeader) {

  # Convert the sent data to the correct object -----------------------------
  dataTable <- convertFactor(dataTable)
  dataTable <- data.table::as.data.table(dataTable)

  # Sort data ---------------------------------------------------------------
  dataTable <- dplyr::arrange(dataTable,
                              dplyr::pull(dataTable, chromHeader),
                              dplyr::pull(dataTable, sampleIdHeader),
                              dplyr::pull(dataTable, positionHeader))

  # Create GRange object ----------------------------------------------------
  ranges <- dataTable %>%
    with(GenomicRanges::GRanges(seqnames = paste(
      dplyr::pull(dataTable, chromHeader),
      dplyr::pull(dataTable, sampleIdHeader)),
      ranges = IRanges::IRanges(
        dplyr::pull(dataTable, positionHeader),
        dplyr::pull(dataTable, positionHeader))
    )
    )

  # Add distance to nearest mutation information to the GRange object -------
  ranges <- addDistance(ranges, maxDistance)

  # Get clusterIDs ----------------------------------------------------------
  tempIds <- paste(dplyr::pull(dataTable, chromHeader),
                   dplyr::pull(dataTable, sampleIdHeader))
  clusterIds <- by(dataTable,
                   factor(tempIds,
                          levels  = unique(tempIds),
                          ordered = TRUE),
                   getClusterId,
                   round(maxDistance), # Just to be sure that the max is a rounded number
                   positionHeader = positionHeader,
                   chromHeader    = chromHeader,
                   sampleIdHeader = sampleIdHeader)

  # Add the cluster information to the dataTable tibble ---------------------
  ranges$clusterId <- unlist(clusterIds)

  # Sort again in the correct order:
  dataTable <- dplyr::arrange(dataTable,
                              dplyr::pull(dataTable, chromHeader),
                              dplyr::pull(dataTable, sampleIdHeader),
                              dplyr::pull(dataTable, positionHeader))
  dataTable <- dplyr::mutate(dataTable,
                             clusterId        = ranges$clusterId,
                             is.clustered     = ranges$is.clustered)
  dataTable <- dplyr::mutate(dataTable,
                             is.clustered     = purrr::map2_lgl(.data$is.clustered,
                                                                !!rlang::sym(refHeader),
                                                                function(x, y) {
                                                                  x & y != "N"
                                                                }),
                             distance         = ranges$distance)
}

#' addIdentTableComment
#' @description Adds commentary to the results of the
#'   \code{\link{identifyClusters}} function.
#' @inheritParams identifyClusters
#' @param linkPatterns Boolean if pattern columns are added
addIdentTableComment <- function(dataTable, linkPatterns) {
  comment(dataTable) <-
    paste0("Information about the added columns:
    distance        : Column with the distances to
                      the nearest mutation.
    clusterId       : Column with the ID of a column.
                      It consist of the chromosome name,
                      sampleID and the sample unique
                      cluster ID number. This is al
                      seperated with a space.
    is.clustered    : Column with Boolean if the mutation
                      is part of a cluster.",ifelse(linkPatterns,"
    linkedPatterns  : Column with the names from the
                      searchPatterns table that matched
                      with the mutation. The names are
                      put in a vector.
    is.linked       : Column with Boolean if there are
                      found patterns for the mutation.",""))
  return(dataTable)
}

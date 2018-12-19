#' getSummaryPatterns
#' @description A optional follow-up function to summarize the results of the
#'   \code{foundPatterns} column from \code{\link{groupClusters}}.
#' @param groupedClusters A table generated from the
#'   \code{\link{groupClusters}}.
#' @inheritParams linkPatterns
#' @param asTibble A Boolean to tell if the returned table needs to be a tibble
#'   or a data.frame.
#' @param renameReverse A Boolean to tell if the pattern IDs from the
#'   searchPatterns also needs to be used to find the patterns with "
#'   [Rev.Com.]". If FALSE and those patterns are available, then they will be
#'   counted in the Unidentified row.
#' @export
#' @note If the \code{groupedClusters} table contains patterns that are not
#'   present in the \code{searchPattern} table, then they will be marked as
#'   \code{Unidentified} together with clusters without patterns.
#' @examples
#' # Example dataset
#' data <- cMut::testDataSet
#'
#' # Use the following functions to get the necessary table
#' results <- identifyClusters(dataTable    = data,
#'                             maxDistance  = 20000,
#'                             linkPatterns = TRUE)
#' groupResults <- groupClusters(dataTable             = results,
#'                               searchClusterPatterns = TRUE,
#'                               patternIntersect      = TRUE)
#'
#' # Use the getSummaryPatterns function to see the summary of found patterns
#' summary <- getSummaryPatterns(groupResults)
#' summary
#'
#' # For more information about the columns use:
#' cat(comment(summary))
#' @importFrom rlang .data
getSummaryPatterns <- function(groupedClusters,
                               searchPatterns = NULL,
                               searchIdHeader = "process",
                               renameReverse  = FALSE,
                               asTibble       = TRUE) {

  # Get or check the searchPatterns table -----------------------------------
  if (is.null(searchPatterns)) {
    # Get the default pattern IDs if nothing is sent:
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(mutationPatterns[ ,"process"]),
                                                   stringsAsFactors = FALSE))
  } else {
    # Check if the assigned headers are present in the given table:
    stopifnot(any(grepl(searchIdHeader, names(searchPatterns))))

    # Get the pattern IDs:
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(searchPatterns[ ,searchIdHeader]),
                                                   stringsAsFactors = FALSE))

  }

  # Add the reverse complement ID's if asked --------------------------------
  if (renameReverse) {
    searchPatterns <- dplyr::bind_rows(searchPatterns,
                                       dplyr::mutate(searchPatterns,
                                                     process = paste0(.data$process,
                                                                      " [Rev.Com.]")))
  }

  # Add an unidentified column and a frequence column to fill up during bootstrapping
  searchPatterns[nrow(searchPatterns) + 1, searchIdHeader] <- "Unidentified"
  searchPatterns <- dplyr::mutate(searchPatterns,
                                  frequency = rep.int(x     = 0,
                                                      times = nrow(searchPatterns)))

  # Count the amount of mutations per pattern per cluster -------------------
  table <- createSummaryPatterns(clusterTable   = groupedClusters,
                                 searchPatterns = searchPatterns,
                                 searchIdHeader = searchIdHeader)

  # Determine the total clustered mutations ---------------------------------
  total <- 0
  for (cMut in groupedClusters$cMuts) {
    total <- total + nrow(cMut)
  }


  # Determine the percentage per pattern over the total ---------------------
  if (total == 0) {
    table <- dplyr::mutate(table,
                           percentage = rep.int(x     = 0,
                                                times = nrow(table)))
  } else {
    table <- dplyr::mutate(table,
                           percentage = purrr::map_dbl(.data$frequency,
                                                       function(x){
                                                         x / total * 100
                                                         }))
  }


  # Explanation about the table ---------------------------------------------
  comment(table) <-
paste0("Information about the summary table columns:
",searchIdHeader,"    : Column with the pattern IDs found
             in sent the searchPattern table.
frequency  : The number of clustered mutations
             that were linked to this pattern ID.
percentage : The frequency divided by the total
             number of clustered mutations times 100.
")

  # return the table in the desired class -----------------------------------
  if (asTibble) {
    return(tibble::as.tibble(table))
  } else {
    return(as.data.frame(table))
  }

}

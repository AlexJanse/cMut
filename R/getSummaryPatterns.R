#' getSummaryPatterns
#' @description A optional folowup function to summarize the results in the
#'   \code{foundPatterns} column from \code{\link{groupClusters}}.
#' @param groupedClusters A table generated from the
#'   \code{\link{groupClusters}}.
#' @inheritParams createSummaryPatterns
#' @param asTibble A Boolean if the returned table needs to be a tibble or a
#'   data.frame
#' @param renameReverse A Boolean if the pattern IDs from the searchPatterns
#'   also need to be used to find the patterns with " [Rev.Com.]". If FALSE and
#'   those patterns are available, then they will be counted in the Unidentified
#'   row. This parameter will be irrelevant if the \code{reverse} parameter is FALSE.
#' @export
#' @note If the \code{groupedClusters} table contains patterns that are not
#'   present in the \code{searchPattern} table, then it will be marked as
#'   unidentified together with clusters without patterns.
#' @examples
#' data <- testDataSet
#' results <- identifyAndAnnotateClusters(data,20000,linkPatterns = TRUE)
#' groupResults <- groupClusters(results,
#'                               searchClusterPatterns = TRUE,
#'                               patternIntersect = TRUE)
#' summary <- getSummaryPatterns(groupResults)
#' summary
getSummaryPatterns <- function(groupedClusters,
                               searchPatterns = NULL,
                               searchIdHeader = "process",
                               renameReverse = FALSE,
                               asTibble = T){

  # get or check the searchPatterns table -----------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(mutationPatterns[,"process"])
                                                ,stringsAsFactors = FALSE))
  } else {
    # check if the assigned headers are present in the given table
    stopifnot(any(grepl(searchIdHeader,names(searchPatterns))))
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(searchPatterns[,searchIdHeader])
                                                   ,stringsAsFactors = FALSE))

  }

  if(renameReverse){
    searchPatterns <- dplyr::bind_rows(searchPatterns,dplyr::mutate(searchPatterns, process = paste0(process," [Rev.Com.]")))
  }

  # Add a unidientified column and a frequence column to fill up during bootstrapping
  searchPatterns[nrow(searchPatterns)+1,searchIdHeader] <- "Unidentified"
  searchPatterns <- dplyr::mutate(searchPatterns, frequency = rep.int(0,nrow(searchPatterns)))

  # Count the amount of mutations per pattern per cluster ------------------------
  table <- createSummaryPatterns(groupedClusters,
                                 searchPatterns,
                                 searchIdHeader, random = FALSE)

  # Determine the total clustered mutations --------------------------------------
  total <- 0
  for(cMut in groupedClusters$cMuts){
    total <- total+nrow(cMut)
  }


  # Determine the percentage per pattern over the total --------------------------
  if(total == 0){
    table <- dplyr::mutate(table, percentage = rep.int(0,nrow(table)))
  } else {
    table <- dplyr::mutate(table, percentage = purrr::map_dbl(frequency,function(x){x/total*100}))
  }

  if(asTibble){
    return(tibble::as.tibble(table))
  } else {
    return(as.data.frame(table))
  }

}

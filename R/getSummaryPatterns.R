#' getSummaryPatterns
#' @description A function to summarise the found patterns
#' @inheritParams createSummaryPatterns
#' @export
#' @examples
#' data <- testDataSet
#' results <- identifyAndAnnotateclusters(data,20000,linkPatterns == T)
getSummaryPatterns <- function(clusterTable,
                             searchPatterns = NULL,
                             searchIdHeader = "process",
                             searchRefHeader = "ref",
                             searchAltHeader = "alt",
                             searchContextHeader = "surrounding",
                             grouped = FALSE,
                             searchReverseComplement = TRUE){

  # get or check the searchPatterns table -----------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(mutationPatterns[,"process"])
                                                ,stringsAsFactors = F))

  } else {
    # check if the assigned headers are present in the given table
    stopifnot(any(grepl(searchIdHeader,names(searchPatterns))))
    searchPatterns <- tibble::as.tibble(data.frame(process = unique(mutationPatterns[,searchIdHeader])
                                                   ,stringsAsFactors = F))

  }

  # Add a frequence column to fill up during bootstrapping
  searchPatterns[nrow(searchPatterns)+1,searchIdHeader] <- "Unidentified"
  searchPatterns <- dplyr::mutate(searchPatterns, frequency = rep.int(0,nrow(searchPatterns)))


  # Count the amount of mutations per pattern per cluster ------------------------
  table <- createSummaryPatterns(clusterTable,
                                 searchPatterns,
                                 searchIdHeader, random = F)

  # Determine the total clustered mutations --------------------------------------
  if(grouped){
    total <- 0
    foreach::foreach(cMut = clusterTable$cMuts) %do% {
      total <- total+nrow(cMut)
    }
  } else {
    total = nrow(clusterTable[clusterTable$is.clustered == T,])
  }

  # Determine the percentage per pattern over the total --------------------------
  if(total == 0){
    table <- dplyr::mutate(table, percentage = rep.int(0,nrow(table)))
  } else {
    table <- dplyr::mutate(table, percentage = purrr::map_dbl(frequency,function(x){x/total*100}))
  }

  return(table)

}

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
                             searchReverseComplement = TRUE,
                             locationBased = FALSE){

  # get or check the searchPatterns table -----------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    if(locationBased){
      searchPatterns <- tibble::as.tibble(data.frame(process = unique(locationPatterns[,"process"])
                                                     ,stringsAsFactors = F))
    } else {
      searchPatterns <- tibble::as.tibble(data.frame(process = unique(mutationPatterns[,"process"])
                                                ,stringsAsFactors = F))
    }

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
                                 searchIdHeader, random = F, locationBased = locationBased)

  # Determine the total clustered mutations --------------------------------------
  total <- 0
  for(cMut in clusterTable$cMuts){
    total <- total+nrow(cMut)
  }


  # Determine the percentage per pattern over the total --------------------------
  if(total == 0){
    table <- dplyr::mutate(table, percentage = rep.int(0,nrow(table)))
  } else {
    table <- dplyr::mutate(table, percentage = purrr::map_dbl(frequency,function(x){x/total*100}))
  }

  return(table)

}

#' createSummaryPatterns
#' @description A function to create a table with frequencies of patterns.
#' @param clusterTable Table with cluster information
#' @inheritParams shuffleMutations
#' @import magrittr
#' @import foreach
createSummaryPatterns <- function(clusterTable,
                                  searchPatterns,
                                  searchIdHeader){


  # Check if there is data is correct ---------------------------------------
  if (nrow(clusterTable) == 0) {
    return(searchPatterns)
  }

  # Peprare the data for the loop -------------------------------------------
  clusterTable <- data.table::as.data.table(clusterTable)
  condition    <- clusterTable[,"has.intersect"]
  if (any(grepl("has.clusterPatterns",
                names(clusterTable)))) {
    condition <- condition | clusterTable[ ,"has.clusterPatterns"]
  }

  # Counter for the unidentified patterns:
  nonIntersectFreq <- 0

  # Start counting the patterns ---------------------------------------------
  for(index in 1:nrow(clusterTable)){

    # Look for the patterns if the row meets the condition
    if(condition[index]){

      nrowCounter  <- 0
      nothingFound <- TRUE
      patterns     <- clusterTable[index,"foundPatterns"][[1]]
      # Loop over the found patterns of the specific row:
      for(pattern in patterns) {
        nrowCounter <- nrowCounter + 1

        addFreq <- nrow(clusterTable[index,"cMuts"][[1]])

        if(any(pattern == dplyr::pull(searchPatterns,
                                      searchIdHeader))){
          nothingFound <- FALSE
          # Get the current frequency in the result pattern table:
          frequency <- searchPatterns[dplyr::pull(searchPatterns,
                                                  searchIdHeader) == pattern,
                                      "frequency"][[1]]
          # Add the new frequency to the result table:
          searchPatterns[dplyr::pull(searchPatterns,
                                     searchIdHeader) == pattern, "frequency"] <-
                                                             frequency + addFreq

        } else {
          if(nrowCounter == length(patterns) & nothingFound){
            # Add the frequency to the unidentified if the found pattern
            #   didn't match with the patterns result table and the cluster
            #   didn't had any other found patterns.
            nonIntersectFreq <- nonIntersectFreq + nrow(clusterTable[index,
                                                                     "cMuts"][[1]])
          }
        }
      }

      # If the row doesn't meet the condition than add the number
      #   of mutations to the nonIntersectFreq:
    } else {
      nonIntersectFreq <- nonIntersectFreq + nrow(clusterTable[index,"cMuts"][[1]])
    }

  }

  # Get the old unidentified frequency
  oldNonIntFreq <- searchPatterns[dplyr::pull(searchPatterns,
                                              searchIdHeader) == "Unidentified",
                                 "frequency"]
  # Save the new unidentified frequency
  searchPatterns[dplyr::pull(searchPatterns,
                             searchIdHeader) == "Unidentified", "frequency"] <-
                                                nonIntersectFreq + oldNonIntFreq

  return(searchPatterns)
}

#' convertFactor
#' @description A function to convert factor columns to character
#' @param x data frame were the columns with string need to converted to
#'   character
convertFactor <- function(x){
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)], as.character)
  return(x)
}

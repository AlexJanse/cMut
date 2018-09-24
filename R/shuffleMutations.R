#' shuffleMutations
#' @description A function to shuffle the reference, alternative and surrounding
#'   nucleotides. Usefull to determine the chance that the results of your data
#'   might be due randomness.
#' @param x A table with the reference, alternative and surrounding
#'   nucleotides.
#' @param refHeader A string with the column header of the reference nucleotide.
#' @param altHeader A string with the column header of the alternative nucleotide.
#' @param surroundingHeader A string with the column header of the surrounding nucleotides.
#' @param nBootstrap A number with the ammount of bootstraps there need to be excecuted.
#' @inheritParams linkPatterns
#' @inheritParams identifyAndAnnotateClusters
#' @import magrittr
#' @import foreach
#' @export
shuffleMutations <- function(x,chromHeader = "chrom",
                             positionHeader = "start",
                             refHeader = "ref",
                             altHeader = "alt",
                             surroundingHeader = "surrounding",
                             sampleIdHeader = "sampleIDs",
                             nBootstrap = 1000,
                             maxDistance = 20000,
                             reverseComplement = FALSE,
                             searchPatterns = NULL,
                             searchRefHeader = "ref",
                             searchAltHeader = "alt",
                             searchContextHeader = "surrounding",
                             searchIdHeader = "proces",
                             searchReverseComplement = TRUE,
                             tibble = TRUE){
  # Get the search table with known mutation patterns -------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    tempWd <- setwd(paste0(.libPaths(),"/cMut")) # To avoid connection issues with the required files
    on.exit(setwd(tempWd), add = T) # Set the wd back to the orignal one
    resultTable <- tibble::as.tibble(readRDS("data/mutationPatterns.rds")) %>%
      dplyr::mutate_all(as.character)
  } else {
    resultTable <- searchPatterns
  }

  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)],as.character)

  # Add a frequence column to fill up during bootstrapping
  resultTable <- dplyr::mutate(resultTable, frequency = rep.int(0,nrow(resultTable)))

  # Preform bootstrap ------------------------------------------------------------
  for(bootstrap in 1:nBootstrap){
    # Create a table with shuffled mutations and contexts ------------------------
    shuffleTable <- createShuffleTable(x,chromHeader,
                                       positionHeader,
                                       refHeader,
                                       altHeader,
                                       surroundingHeader,
                                       sampleIdHeader)

    # Add the frequencies of patterns to the resultTable -----------------------
    resultTable <- addToResultTable(chromHeader = "chrom",
                                     positionHeader = "pos",
                                     refHeader = "ref",
                                     altHeader = "alt",
                                     surroundingHeader = "surrounding",
                                     sampleIdHeader = "sampleIDs",
                                     maxDistance,
                                     shuffleTable,
                                     reverseComplement,
                                     searchPatterns = resultTable,
                                     searchRefHeader,
                                     searchAltHeader,
                                     searchContextHeader,
                                     searchIdHeader,
                                     searchReverseComplement)
  }

  # Calculate the percentage of the frequecies -----------------------------------
  total = sum(resultTable$frequency)
  results <- resultTable %>% dplyr::mutate(percentage = purrr::map_dbl(frequency,function(x){x/total*100}))

  if(tibble){
    return(tibble::as.tibble(results))
  } else {
    return(as.data.frame(results))
  }

}

#' createShuffleTable
#' @description A function to bootstrap the mutation table
#' @inheritParams shuffleMutations
#' @import magrittr
#' @import foreach
createShuffleTable <-  function(x,chromHeader,
                                positionHeader,
                                refHeader,
                                altHeader,
                                surroundingHeader,
                                sampleIdHeader){

  # Extract the muation data from the table --------------------------------------
  ref <- dplyr::pull(x,refHeader)
  alt <- dplyr::pull(x,altHeader)
  surrounding <- dplyr::pull(x,surroundingHeader)

  # Create a table with random picked reference and surrounding nucleotides -------
  shuffleTable <- data.frame(ref = sample(ref,length(ref)),stringsAsFactors = F) %>%
    dplyr::mutate(surrounding = sample(surrounding,length(surrounding)))

  again <- TRUE
  while(again){
    again <- FALSE

    shuffleAlt <- c()
    alt <- dplyr::pull(x,altHeader)

    # A loop to chosse a random alternative nucleotide while making sure that it's not the same as the reference ---------
    foreach::foreach(index = 1:length(ref)) %do% {
      grabbag <- alt

      refNuc = shuffleTable[index,1][[1]]
      grabbag <- setdiff(grabbag,c(refNuc))
      if(length(grabbag) == 0){
        again = TRUE
      } else{
        random <- as.character(sample(grabbag,1)[[1]])
      }

      alt <- alt[-match(random,alt)]
      shuffleAlt[length(shuffleAlt)+1] <- random
    }
  }

  shuffleTable <- shuffleTable %>%
    dplyr::mutate(sampleIDs = dplyr::pull(x,sampleIdHeader),
                  chrom = dplyr::pull(x,chromHeader),
                  pos = dplyr::pull(x,positionHeader),
                  alt = shuffleAlt,
                  check = !alt == ref)

  nTable <- dplyr::count(dplyr::group_by(shuffleTable[c("sampleIDs","chrom","pos","ref","alt","surrounding","check")],
                                         sampleIDs,
                                         chrom,
                                         pos,
                                         ref,
                                         alt,
                                         surrounding,
                                         check))

  tempDF <- data.frame(stringsAsFactors = f)
  shuffleTable <- dplyr::bind_rows(tempDF,nTable) # Needed to avoid a grouped table which will lead to errors further down the path

  return(shuffleTable)
}

#' addToResultTable
#' @description A function to create a table with frequencies of patterns.
#' @inheritParams shuffleMutations
#' @import magrittr
#' @import foreach
addToResultTable <- function(chromHeader,
                              positionHeader,
                              refHeader,
                              altHeader,
                              surroundingHeader,
                              sampleIdHeader,
                              maxDistance,
                              shuffleTable,
                              reverseComplement,
                              searchPatterns,
                              searchRefHeader,
                              searchAltHeader,
                              searchContextHeader,
                              searchIdHeader,
                              searchReverseComplement){

  # Add patterns that match with the random mutations in the shuffleTable ----------------
  clusterTablePerMut <- identifyAndAnnotateClusters(x = shuffleTable,maxDistance = maxDistance,positionHeader = "pos",linkPatterns = T)
  clusterTable <- groupClusters(clusterTablePerMut,patternIntersect = T)

  for(index in 1:nrow(clusterTable)){
    if(clusterTable[index,"has.intersect"] == TRUE){
      for(patroon in clusterTable[index,"patternIntersect"][[1]][[1]]){
        frequency <- searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == patroon,"frequency"][[1]]
        addFreq <- sum(clusterTable[index,"cMuts"][[1]][[1]][,"n"])
        searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == patroon,"frequency"] <- frequency + addFreq
      }
    }
  }
  return(searchPatterns)
}
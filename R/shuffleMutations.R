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
#' @param showProgress A boolean whether or not to show the progress.
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
                             showProgress = TRUE,
                             linkPatterns = FALSE,
                             maxDistance = 20000,
                             reverseComplement = FALSE,
                             searchPatterns = NULL,
                             searchRefHeader = "ref",
                             searchAltHeader = "alt",
                             searchContextHeader = "surrounding",
                             searchIdHeader = "proces",
                             searchReverseComplement = TRUE,
                             tibble = TRUE){

  shuffleTable <- createShuffleTable(x,chromHeader,
                                     positionHeader,
                                     refHeader,
                                     altHeader,
                                     surroundingHeader,
                                     sampleIdHeader,
                                     nBootstrap,
                                     showProgress)

  if(linkPatterns == TRUE){
    # Create a table with the patterns and their frequency ----------------------
    patternFreq <- createPatternFreq(chromHeader = "chrom",
                                     positionHeader = "pos",
                                     refHeader = "ref",
                                     altHeader = "alt",
                                     surroundingHeader = "surrounding",
                                     sampleIdHeader = "sampleIDs",
                                     maxDistance,
                                     shuffleTable,
                                     linkPatterns,
                                     reverseComplement,
                                     searchPatterns,
                                     searchRefHeader,
                                     searchAltHeader,
                                     searchContextHeader,
                                     searchIdHeader,
                                     searchReverseComplement)

    total = sum(patternFreq$frequency)
    results <- dplyr::mutate(patternFreq, percentage = frequency/total*100)
  } else {
    # Add the percentage to the shuffleTable ------------------------------------
    total = sum(shuffleTable$n)
    results <- dplyr::mutate(shuffleTable,percentage = n/total*100)
  }


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
                                sampleIdHeader,
                                nBootstrap,
                                showProgress){
  shuffleTable <- data.frame(chom = c(), pos = c(),sampleIDs = c(), stringsAsFactors = F)

  foreach::foreach(i = 1:nBootstrap) %do% {
    if(showProgress){
      svMisc::progress(i/nBootstrap*100)
    }
    ref <- dplyr::pull(x,refHeader)

    surrounding <- dplyr::pull(x,surroundingHeader)

    tempShuffleTable <- data.frame(ref = sample(ref,length(ref)),stringsAsFactors = F) %>%
      dplyr::mutate(surrounding = sample(surrounding,length(surrounding)))

    again <- TRUE
    while(again){
      again <- FALSE
      shuffleAlt <- c()
      alt <- dplyr::pull(x,altHeader)

      foreach::foreach(index = 1:length(ref)) %do% {
        grabbag <- alt
        random <- sample(grabbag,1)[[1]]
        while(length(grabbag) != 0 & random == tempShuffleTable[index,1][[1]]){
          grabbag <- grabbag[-match(random,grabbag)]
          if(length(grabbag) > 0){
            random <- sample(grabbag,1)[[1]]
          }
        }

        if(length(grabbag) == 0 & random == tempShuffleTable[index,1][[1]]){
          again = TRUE
        }

        alt <- alt[-match(random,alt)]
        shuffleAlt[length(shuffleAlt)+1] <- random
      }
    }

    tempShuffleTable <- tempShuffleTable %>%
      dplyr::mutate(sampleIDs = dplyr::pull(x,sampleIdHeader),
                    chrom = dplyr::pull(x,chromHeader),
                    pos = dplyr::pull(x,positionHeader),
                    alt = shuffleAlt,
                    check = !alt == ref)

    shuffleTable <- dplyr::bind_rows(shuffleTable,tempShuffleTable)
  }
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
  if(showProgress){
    print("Bootstrap done!")
  }

  return(shuffleTable)
}

#' createPatternFreq
#' @description A function to create a table with frequencies of patterns.
#' @inheritParams shuffleMutations
#' @import magrittr
#' @import foreach
createPatternFreq <- function(chromHeader,
                              positionHeader,
                              refHeader,
                              altHeader,
                              surroundingHeader,
                              sampleIdHeader,
                              maxDistance,
                              shuffleTable,
                              linkPatterns,
                              reverseComplement,
                              searchPatterns,
                              searchRefHeader,
                              searchAltHeader,
                              searchContextHeader,
                              searchIdHeader,
                              searchReverseComplement){

  # Add patterns that match with the random mutations in the shuffleTable ----------------
  clusterTablePerMut <- identifyAndAnnotateClusters(x = shuffleTable,maxDistance = maxDistance,positionHeader = "pos",linkPatterns = T)
  print("ident")
  clusterTable <- groupClusters(clusterTablePerMut,patternIntersect = T)
  print("clusterd")
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent ----------------------------------------
    tempWd <- setwd(paste0(.libPaths(),"/cMut")) # To avoid connection issues with the required files
    on.exit(setwd(tempWd), add = T) # Set the wd back to the orignal one
    searchPatterns <- tibble::as.tibble(readRDS("data/mutationPatterns.rds")) %>%
      dplyr::mutate_all(as.character)
  }

  searchPatterns$frequency <- rep.int(0,nrow(searchPatterns)) # Add a column to the pattern table for the frequencies

  foreach::foreach(index = 1:nrow(clusterTable)) %do% {
    if(clusterTable[index,"has.intersect"] == TRUE){
      foreach::foreach(patroon = clusterTable[index,"patternIntersect"][[1]][[1]]) %do% {
        frequency <- searchPatterns[grep(searchIdHeader, names(searchPatterns)) == patroon,"frequency"]
        addFreq <- sum(clusterTable[index,"cMuts"][[1]][[1]][,"n"])
        searchPatterns[grep(searchIdHeader, names(searchPatterns)) == patroon,"frequency"] <- frequency + addFreq
      }
    }
  }
  return(searchPatterns)
}

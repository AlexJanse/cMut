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
                             searchIdHeader = "process",
                             searchReverseComplement = TRUE,
                             tibble = TRUE){

  # Get the search table with known mutation patterns -------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    resultTable <- mutationPatterns
  } else {
    resultTable <- searchPatterns
  }

  # Add the reverse complement of the known table to the search table -----------------------------------------------
  if(searchReverseComplement){
    i <- getRevComTable(resultTable,searchRefHeader,searchAltHeader,searchContextHeader,searchIdHeader)
    resultTable <- rbind(resultTable,i)
  }


  x <- convertFactor(x)

  # Add a frequence column to fill up during bootstrapping
  resultTable[nrow(resultTable)+1,searchIdHeader] <- "Unidentified"
  resultTable <- tibble::tibble(!!rlang::sym(searchIdHeader) := unique(resultTable[,searchIdHeader]))
  resultTable <- dplyr::mutate(resultTable, frequency = rep.int(0,nrow(resultTable)))

  total <- 0 # Variable to keep up the total nrows
  # Preform bootstrap ------------------------------------------------------------
  for(bootstrap in 1:nBootstrap){
    # Create a table with shuffled mutations and contexts ------------------------
    shuffleTable <- createShuffleTable(x,chromHeader,
                                       positionHeader,
                                       refHeader,
                                       altHeader,
                                       surroundingHeader,
                                       sampleIdHeader)

    # Identify, annotate and group clustered mutations --------------------------
    clusterTablePerMut <- identifyAndAnnotateClusters(x = shuffleTable,
                                                      maxDistance = maxDistance,
                                                      positionHeader = "pos",
                                                      linkPatterns = T)
    clusterTable <- groupClusters(clusterTablePerMut,patternIntersect = T, showWarning = F)


    total <- total+nrow(clusterTablePerMut[clusterTablePerMut$is.clustered == T,])


    # Add the frequencies of patterns to the resultTable -----------------------
    resultTable <- createSummaryPatterns(clusterTable,
                                     searchPatterns = resultTable,
                                     searchIdHeader,
                                     random = T)
  }

  # Calculate the percentage of the frequecies -----------------------------------



  if(total != 0){
    results <- dplyr::mutate(resultTable, percentage = purrr::map_dbl(frequency,function(x){x/total*100}))
  } else {
    results <- dplyr::mutate(resultTable, percentage = rep.int(0, nrow(resultTable)))
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
                                sampleIdHeader){

  # Extract the muation data from the table --------------------------------------
  ref <- dplyr::pull(x,refHeader)
  alt <- dplyr::pull(x,altHeader)
  surrounding <- dplyr::pull(x,surroundingHeader)

  # Create a table with random picked reference and surrounding nucleotides -------
  shuffleTable <- data.table::data.table(ref = sample(ref,length(ref)),stringsAsFactors = F) %>%
    dplyr::mutate(surrounding = sample(surrounding,length(surrounding)))

  again <- TRUE
  while(again){
    again <- FALSE

    shuffleAlt <- c()
    alt <- dplyr::pull(x,altHeader)

    # A loop to chosse a random alternative nucleotide while making sure that it's not the same as the reference ---------
    for(index in 1:length(ref)) {
      grabbag <- alt
      refNuc = shuffleTable[index,1]
      grabbag <- setdiff(grabbag,c(refNuc))
      if(length(grabbag) == 0){
        indexPrev <- index-1
        failed <- TRUE
        while(indexPrev != 0 & failed){
          refNucPrev <- shuffleTable[indexPrev,1]
          altPrev <- shuffleAlt[indexPrev]
          if(refNuc != refNucPrev &
             altPrev != refNuc){
            shuffleAlt[indexPrev] <- refNuc
            shuffleAlt[index] <- altPrev
            alt <- alt[-match(refNuc,alt)]
            failed <- FALSE
          }
          indexPrev <- indexPrev-1
        }
        if(failed){
          again <- TRUE        }
      } else{
        random <- as.character(sample(grabbag,1)[[1]])
        alt <- alt[-match(random,alt)]
        shuffleAlt[length(shuffleAlt)+1] <- random
      }
    }
  }

  # Fill the table with mutation information ----------------------------------
  shuffleTable <- shuffleTable %>%
    dplyr::mutate(sampleIDs = dplyr::pull(x,sampleIdHeader),
                  chrom = dplyr::pull(x,chromHeader),
                  pos = dplyr::pull(x,positionHeader),
                  alt = shuffleAlt,
                  check = !alt == ref)

  # Create a table with the frequency of each mutation ------------------------
  nTable <- dplyr::count(dplyr::group_by(shuffleTable[c("sampleIDs","chrom","pos","ref","alt","surrounding","check")],
                                         sampleIDs,
                                         chrom,
                                         pos,
                                         ref,
                                         alt,
                                         surrounding,
                                         check))

  tempDF <- data.table::data.table(stringsAsFactors = f)
  shuffleTable <- dplyr::bind_rows(tempDF,nTable) # Needed to avoid a grouped table which will lead to errors further down the path

  return(shuffleTable)
}

#' createSummaryPatterns
#' @description A function to create a table with frequencies of patterns.
#' @param clusterTable Table with cluster information
#' @param grouped A Boolean to tell if the data is grouped or not
#' @inheritParams shuffleMutations
#' @import magrittr
#' @import foreach
createSummaryPatterns <- function(clusterTable,
                             searchPatterns,
                             searchIdHeader,
                             random = FALSE){

  if(nrow(clusterTable) == 0){
    return(searchPatterns)
  }
  nonIntersectFreq <- 0
  for(index in 1:nrow(clusterTable)){
    if(clusterTable[index,"has.intersect"][[1]] == TRUE){
      for(pattern in clusterTable[index,"patternIntersect"][[1]][[1]]) {
        if(random){
          addFreq <- sum(clusterTable[index,"cMuts"][[1]][[1]][,"n"])
        } else {
          addFreq <- nrow(clusterTable[index,"cMuts"][[1]][[1]])
        }
        frequency <- searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"][[1]]
        searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"] <- frequency + addFreq
      }
    } else {
      nonIntersectFreq <- nonIntersectFreq + nrow(clusterTable[index,"cMuts"][[1]][[1]])
    }
  }
  oldNonIntFreq <- searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == "Unidentified","frequency"]
  searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == "Unidentified","frequency"] <- nonIntersectFreq + oldNonIntFreq
  return(searchPatterns)
}

#' convertFactor
#' @description A function to convert factor columns to character
convertFactor <- function(x){
  x[sapply(x, is.factor)] <- lapply(x[sapply(x, is.factor)],as.character)
  return(x)
}



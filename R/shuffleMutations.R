#' shuffleMutations
#' @description A function to shuffle the reference, alternative and surrounding
#'   nucleotides. Then it will use the \code{\link{identifyAndAnnotateClusters}}
#'   and \code{\link{groupClusters}} functions and returns a summary of the
#'   frequency of found patterns.
#' @param x A table with the reference, alternative and surrounding nucleotides.
#'   The best data to use is the output of the
#'   \code{\link{identifyAndAnnotateClusters}} where \code{is.clustered} is
#'   TRUE.
#' @param refHeader A string with the column header of the reference nucleotide.
#' @param altHeader A string with the column header of the alternative
#'   nucleotide.
#' @param surroundingHeader A string with the column header of the surrounding
#'   nucleotides.
#' @param nBootstrap A number with the ammount of bootstraps there need to be
#'   excecuted.
#' @inheritParams linkPatterns
#' @inheritParams identifyAndAnnotateClusters
#' @param saveEachBootstrap A Boolean if the summaries per bootstrap are needed
#'   to be saved
#' @param saveFileName A string with the name and location of the file that
#'   contains the idivudual summaries. Only needed if saveEachBootstrap is TRUE.
#' @param no.cores A number with the amount of clusters that is allowed to
#'   use during shuffle. Default is maximum amount of cores present on the
#'   system.
#' @import magrittr
#' @import foreach
#' @import doParallel
#' @import compiler
#' @export
#' @examples
#' identResults <- identifyAndAnnotateClusters(testDataSet,20000, linkPatterns = T)
#'
#' # If only the mutation patterns are needed searched:
#' shuffleResults <- shuffleMutations(identResults,nBootstrap = 5)
#'
#' # If also the cluster patterns are needed to be added:
#' shuffleResults <- shuffleMutations(identResults,
#'                                    nBootstrap = 5,
#'                                    searchClusterPatterns = T)
#'
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
                             searchDistanceHeader = "maxDistance",
                             searchReverseComplement = TRUE,
                             tibble = TRUE,
                             saveEachBootstrap = FALSE,
                             saveFileName = NULL,
                             searchClusterPatterns = FALSE,
                             no.cores = parallel::detectCores()){

  stopifnot(no.cores > 0)

  # Get the search table with known mutation patterns -------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    resultTable <- getSearchPatterns(reverse = F, asTibble = F)
    searchPatterns <- getSearchPatterns(reverse = F, asTibble = F)
  } else {
    resultTable <- convertFactor(as.data.frame(searchPatterns))
    searchPatterns <- convertFactor(as.data.frame(searchPatterns))
  }

  if(!searchClusterPatterns){
    searchPatterns <- searchPatterns[nchar(dplyr::pull(searchPatterns, searchRefHeader)) == 1,]
    resultTable <- resultTable[nchar(dplyr::pull(resultTable, searchRefHeader)) == 1,]
  }

  # Add the reverse complement of the known table to the search table -----------------------------------------------
  if(searchReverseComplement){
    i <- getRevComTable(resultTable,searchRefHeader,searchAltHeader,searchContextHeader,searchIdHeader)
    resultTable <- rbind(resultTable,i)
    searchPatterns <- dplyr::bind_rows(searchPatterns, getRevComTable(searchPatterns,searchRefHeader,searchAltHeader,searchContextHeader,searchIdHeader))
  }


  x <- convertFactor(x)

  # Add a frequence column to fill up during bootstrapping
  resultTable[nrow(resultTable)+1,searchIdHeader] <- "Unidentified"
  resultTable <- tibble::tibble(!!rlang::sym(searchIdHeader) := unique(resultTable[,searchIdHeader]))
  resultTable <- dplyr::mutate(resultTable, frequency = rep.int(0,nrow(resultTable)))

  # Prepare for parallel loop
  if(no.cores > parallel::detectCores()){
   stop("The no.cores parameter is higher than amount of clusters present.")
  }
  clusters <- parallel::makeCluster(no.cores)
  doParallel::registerDoParallel(clusters)

  # Preform bootstrap ------------------------------------------------------------
  resultTables <- foreach::foreach(iterators::icount(nBootstrap)) %dopar% {
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
                                                      linkPatterns = TRUE,
                                                      searchPatterns = searchPatterns,
                                                      searchRefHeader = searchRefHeader,
                                                      searchAltHeader = searchAltHeader,
                                                      searchContextHeader = searchContextHeader,
                                                      searchIdHeader = searchIdHeader,
                                                      searchDistanceHeader = searchDistanceHeader,
                                                      searchReverseComplement = FALSE)

    clusterTable <- groupClusters(clusterTablePerMut,
                                  patternIntersect = TRUE,
                                  showWarning = F,
                                  searchClusterPatterns = searchClusterPatterns,
                                  searchPatterns = searchPatterns[nchar(dplyr::pull(searchPatterns,searchRefHeader)) > 1,],
                                  searchRefHeader = searchRefHeader,
                                  searchAltHeader = searchAltHeader,
                                  searchIdHeader = searchIdHeader,
                                  searchDistanceHeader = searchDistanceHeader,
                                  searchReverseComplement = FALSE)


    # Add the frequencies of patterns to the resultTable -----------------------
    subResultTable <- resultTable
    subResultTable <- createSummaryPatterns(clusterTable,
                                            searchPatterns = subResultTable,
                                            searchIdHeader,
                                            random = T)

    total <- list("total",nrow(clusterTablePerMut[clusterTablePerMut$is.clustered == T,]))
    subResultTable <- rbind(subResultTable,total)
  }
  parallel::stopCluster(clusters)

  if(saveEachBootstrap){
    dput(resultTables,saveFileName)
  }

  # Calculate the percentage of the frequecies -----------------------------------
  total <- list("total",0)
  resultTable <- rbind(resultTable,total)
  for(table in resultTables){
    resultTable <- data.table::as.data.table(tibble::tibble(!!rlang::sym(searchIdHeader) := dplyr::pull(table,!!rlang::sym(searchIdHeader)),
                                          frequency = resultTable$frequency + table$frequency))
  }

  total <- resultTable[dplyr::pull(resultTable,!!rlang::sym(searchIdHeader)) == "total",2][[1]]
  resultTable <- resultTable[dplyr::pull(resultTable,!!rlang::sym(searchIdHeader)) != "total",]


  if(total != 0){
    results <- dplyr::mutate(resultTable, percentage = purrr::map_dbl(frequency,function(x){x/total*100}))
  } else {
    results <- dplyr::mutate(resultTable, percentage = rep.int(0, nrow(resultTable)))
  }

  comment(results) <-
    "
    Information about the summary table columns:
    *searchIdHeader*     : Column with the pattern IDs from
                           the searchPatterns table.
    frequency            : Column with the number of mutations
                           that had this pattern ID as intersect
                           of their cluster.
    percentage           : Column with the percentage of the
                           frequency. Please notice that it is
                           possible to have a total percentage
                           above 100% since clusters may be
                           linked with multiple patterns. However
                           the total percentage should never be
                           below 100%.
  "

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
                                sampleIdHeader
                                ){

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
  clusterTable <- data.table::as.data.table(clusterTable)
  condition <- clusterTable[,"has.intersect"]
  if(length(which(clusterTable == "has.clusterPatterns")) != 0){
    condition <- condition | clusterTable[,"has.clusterPatterns"]
  }
  nonIntersectFreq <- 0

  for(index in 1:nrow(clusterTable)){
    if(condition[index]){
      for(pattern in clusterTable[index,"foundPatterns"][[1]]) {
        if(random){
          addFreq <- sum(clusterTable[index,"cMuts"][[1]][,"n"])
        } else {
          addFreq <- nrow(clusterTable[index,"cMuts"][[1]])
        }
        frequency <- searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"][[1]]
        searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"] <- frequency + addFreq
      }

    } else {
      nonIntersectFreq <- nonIntersectFreq + nrow(clusterTable[index,"cMuts"][[1]])
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



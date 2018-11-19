#' shuffleMutations
#' @description A function to shuffle the reference, alternative and surrounding
#'   nucleotides. Then it will use the \code{\link{identifyAndAnnotateClusters}}
#'   and \code{\link{groupClusters}} functions and returns a summary of the
#'   frequency of found patterns.
#' @param dataTable A table with the reference, alternative and surrounding
#'   nucleotides. The best data to use is the output of the
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
#' @param returnEachBootstrap A Boolean if the summaries per bootstrap are
#'   needed to be returned. If TRUE then the return value will be a list with a
#'   summary of each bootstrap.
#' @param no.cores A number with the amount of clusters that is allowed to use
#'   during shuffle. Default is maximum amount of cores present on the system.
#' @import magrittr
#' @import foreach
#' @import doParallel
#' @import compiler
#' @import doSNOW
#' @export
#' @examples
#' identResults <- identifyAndAnnotateClusters(testDataSet,20000, linkPatterns = TRUE)
#'
#' # If only the mutation patterns are needed searched:
#' shuffleResults <- shuffleMutations(identResults[identResults$is.clustered,],
#'                                    nBootstrap = 5,
#'                                    no.cores = 2)
#'
#' # If also the cluster patterns are needed to be added:
#' shuffleResults <- shuffleMutations(identResults[identResults$is.clustered,],
#'                                    nBootstrap = 5,
#'                                    searchClusterPatterns = TRUE,
#'                                    no.cores = 2)
#' # The no.cores is set to 2 because of CRAN limits when testing the examples.
shuffleMutations <- function(dataTable,chromHeader = "chrom",
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
                             asTibble = TRUE,
                             returnEachBootstrap = FALSE,
                             searchClusterPatterns = FALSE,
                             no.cores = parallel::detectCores(),
                             renameReverse = FALSE){


  # Get the search table with known mutation patterns -------------------------------
  if(is.null(searchPatterns)){
    # Get default table if nothing is sent
    resultTable <- getSearchPatterns(reverse = FALSE, asTibble = FALSE)
    searchPatterns <- getSearchPatterns(reverse = FALSE, asTibble = FALSE)
  } else {
    resultTable <- convertFactor(as.data.frame(searchPatterns))
    searchPatterns <- convertFactor(as.data.frame(searchPatterns))
  }

  checkParametersShuffleMutations(dataTable,chromHeader, positionHeader,
                                  refHeader,altHeader,surroundingHeader,
                                  sampleIdHeader,nBootstrap,maxDistance,
                                  reverseComplement,searchPatterns,searchRefHeader,
                                  searchAltHeader,searchContextHeader,searchIdHeader,
                                  searchDistanceHeader,searchReverseComplement,
                                  asTibble,returnEachBootstrap,
                                  searchClusterPatterns,no.cores)

  if(!searchClusterPatterns){
    searchPatterns <- searchPatterns[nchar(dplyr::pull(searchPatterns, searchRefHeader)) == 1,]
    resultTable <- resultTable[nchar(dplyr::pull(resultTable, searchRefHeader)) == 1,]
  }

  # Add the reverse complement of the known table to the search table -----------------------------------------------
  if(searchReverseComplement){
    i <- getRevComTable(resultTable,searchRefHeader,searchAltHeader,
                        searchContextHeader,searchIdHeader,
                        renameReverse = renameReverse)
    resultTable <- rbind(resultTable,i)
    searchPatterns <- dplyr::bind_rows(searchPatterns,
                                       getRevComTable(searchPatterns,searchRefHeader,
                                                      searchAltHeader,searchContextHeader,
                                                      searchIdHeader, renameReverse = renameReverse))
  }


  dataTable <- convertFactor(dataTable)

  # Add a frequence column to fill up during bootstrapping
  resultTable[nrow(resultTable)+1,searchIdHeader] <- "Unidentified"
  resultTable <- tibble::tibble(!!rlang::sym(searchIdHeader) := unique(resultTable[,searchIdHeader]))
  resultTable <- dplyr::mutate(resultTable, frequency = rep.int(0,nrow(resultTable)))

  resultTables <- shuffleParallel(dataTable,chromHeader, positionHeader,
                                  refHeader,altHeader,surroundingHeader,
                                  sampleIdHeader,nBootstrap,maxDistance,
                                  reverseComplement,searchPatterns,searchRefHeader,
                                  searchAltHeader,searchContextHeader,searchIdHeader,
                                  searchDistanceHeader,searchReverseComplement,
                                  searchClusterPatterns,no.cores,resultTable)

  if(returnEachBootstrap){
    return(resultTables)
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

  if(asTibble){
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
  if(any(grepl("has.clusterPatterns",names(clusterTable)))){
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
        if(any(grepl(pattern,dplyr::pull(searchPatterns,searchIdHeader),fixed = TRUE))){
          frequency <- searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"][[1]]
          searchPatterns[searchPatterns[,grep(searchIdHeader, names(searchPatterns))] == pattern,"frequency"] <- frequency + addFreq

        } else {
          nonIntersectFreq <- nonIntersectFreq + nrow(clusterTable[index,"cMuts"][[1]])
        }
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

#' checkParametersShuffleMutations
#' @description A function to check if the parameters are correct of the
#'   \code{\link{shuffleMuations}} function
checkParametersShuffleMutations <- function(dataTable,chromHeader, positionHeader,
                                            refHeader,altHeader,surroundingHeader,
                                            sampleIdHeader,nBootstrap,maxDistance,
                                            reverseComplement,searchPatterns,searchRefHeader,
                                            searchAltHeader,searchContextHeader,searchIdHeader,
                                            searchDistanceHeader,searchReverseComplement,
                                            asTibble,returnEachBootstrap,
                                            searchClusterPatterns,no.cores){
  # Check headers ------------------------------------------
  sapply(c(chromHeader,
           positionHeader,
           refHeader,
           altHeader,
           surroundingHeader),
         function(x){
           if(!any(x == names(dataTable))){
             stop("Please check if the column names in the parameters match with the dataTable columns.")
           }
         })

  sapply(c(searchAltHeader,
           searchRefHeader,
           searchContextHeader,
           searchIdHeader),
         function(x){
           if(!any(x == names(searchPatterns))){
             stop("Please check if the column names in the parameters match with the searchPatterns table columns.")
           }
         })

  # Check numeric values ------------------------------------
  stopifnot(is.numeric(maxDistance))
  stopifnot(no.cores > 0)

  if(no.cores > parallel::detectCores()){
    stop("The no.cores parameter is higher than amount of clusters present.")
  }

  # Check Boolean values --------------------------------------
  stopifnot(is.logical(searchClusterPatterns) & is.logical(reverseComplement) & is.logical(searchReverseComplement))
}


#' shuffleParallel
#' @description Function to preform the shuffling and annotating in parallel
shuffleParallel <- function(dataTable,chromHeader, positionHeader,
                            refHeader,altHeader,surroundingHeader,
                            sampleIdHeader,nBootstrap,maxDistance,
                            reverseComplement,searchPatterns,searchRefHeader,
                            searchAltHeader,searchContextHeader,searchIdHeader,
                            searchDistanceHeader,searchReverseComplement,
                            searchClusterPatterns,no.cores,resultTable){

  # Prepare for parallel loop
  clusters <- parallel::makeCluster(no.cores)
  doSNOW::registerDoSNOW(clusters)
  bar <- txtProgressBar(max = nBootstrap, style = 3)
  progress <- function(x){setTxtProgressBar(bar,x)}
  opts <- list(progress = progress)

  # Preform bootstrap ------------------------------------------------------------
  resultTables <- foreach::foreach(iterators::icount(nBootstrap),
                                   .verbose = TRUE,
                                   .packages = c("magrittr","cMut"),
                                   .export = c("createShuffleTable"),
                                   .options.snow = opts) %dopar% {

                                     # Create a table with shuffled mutations and contexts ------------------------
                                     shuffleTable <- createShuffleTable(dataTable,chromHeader,
                                                                        positionHeader,
                                                                        refHeader,
                                                                        altHeader,
                                                                        surroundingHeader,
                                                                        sampleIdHeader)

                                     # Identify, annotate and group clustered mutations --------------------------
                                     clusterTablePerMut <- identifyAndAnnotateClusters(dataTable = shuffleTable,
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
                                                                   showWarning = FALSE,
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
                                                                             random = TRUE)

                                     total <- list("total",nrow(clusterTablePerMut[clusterTablePerMut$is.clustered == TRUE,]))
                                     subResultTable <- rbind(subResultTable,total)
                                   }
  parallel::stopCluster(clusters)
  return(resultTables)
}


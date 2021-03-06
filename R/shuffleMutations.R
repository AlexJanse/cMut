#' shuffleMutations
#' @description A function to shuffle the reference, alternative and surrounding
#'   nucleotides. Then it will use the \code{\link{identifyClusters}}
#'   and \code{\link{groupClusters}} functions and returns a summary of the
#'   frequency of found patterns.
#' @param dataTable A table with the reference, alternative and surrounding
#'   nucleotides. The best data to use is the output of the
#'   \code{\link{identifyClusters}} where \code{is.clustered} is
#'   TRUE.
#' @param refHeader A string with the column header of the reference nucleotide.
#' @param altHeader A string with the column header of the alternative
#'   nucleotide.
#' @param contextHeader A string with the column header of the context
#'   nucleotides.
#' @param nBootstrap A number with the amount of bootstraps there need to be
#'   executed.
#' @inheritParams linkPatterns
#' @inheritParams identifyClusters
#' @inheritParams groupClusters
#' @param asTibble A Boolean if the returned results needs to be a tibble. It
#'   will return a data.frame otherwise. Irrelevant if \code{returnEachBootstrap
#'   == TRUE}
#' @param returnEachBootstrap A Boolean if the summaries per bootstrap are
#'   needed to be returned. If TRUE then the return value will be a list with a
#'   summary of each bootstrap.
#' @param no.cores A number with the amount of clusters that is allowed to use
#'   during shuffle. Default is maximum amount of cores present on the system.
#' @importFrom magrittr %>%
#' @import foreach
#' @import doParallel
#' @import compiler
#' @importFrom rlang :=
#' @export
#' @examples
#' identResults <- identifyClusters(dataTable    = cMut::testDataSet,
#'                                  maxDistance  = 20000,
#'                                  linkPatterns = TRUE)
#' clusteredMutations <- identResults[identResults$is.clustered, ]
#'
#' # If only the mutation patterns are needed searched:
#' \donttest{shuffleResults <- shuffleMutations(dataTable             = clusteredMutations,
#'                                    searchClusterPatterns = FALSE,
#'                                    no.cores              = 2)}
#'
#' # If also the cluster patterns are needed to be added:
#' \donttest{shuffleResults <- shuffleMutations(dataTable  = identResults[identResults$is.clustered,],
#'                                     no.cores   = 2)}
shuffleMutations <- function(dataTable,                             chromHeader             = "chrom",
                             positionHeader        = "start",       refHeader               = "ref",
                             altHeader             = "alt",         contextHeader           = "surrounding",
                             sampleIdHeader        = "sampleIDs",   nBootstrap              = 1000,
                             maxDistance           = 20000,         linkPatterns            = TRUE,
                             reverseComplement     = FALSE,         searchPatterns          = NULL,
                             searchRefHeader       = "ref",         searchAltHeader         = "alt",
                             searchContextHeader   = "surrounding", searchIdHeader          = "process",
                             searchDistanceHeader  = "maxDistance", searchReverseComplement = TRUE,
                             asTibble              = TRUE,          returnEachBootstrap     = FALSE,
                             searchClusterPatterns = TRUE,          renameReverse           = FALSE,
                             no.cores              = parallel::detectCores()) {




  # Prepare the parameters for the functions table --------------------------
  # Get and/or adjust the search pattern table:
  if (is.null(searchPatterns)) {
    resultTable    <- getSearchPatterns(reverse = FALSE, asTibble = FALSE)
    searchPatterns <- getSearchPatterns(reverse = FALSE, asTibble = FALSE)
  } else {
    resultTable    <- convertFactor(as.data.frame(searchPatterns))
    searchPatterns <- convertFactor(as.data.frame(searchPatterns))
  }

  checkParametersShuffleMutations(dataTable               = dataTable,               chromHeader           = chromHeader,
                                  positionHeader          = positionHeader,          refHeader             = refHeader,
                                  altHeader               = altHeader,               contextHeader         = contextHeader,
                                  sampleIdHeader          = sampleIdHeader,          nBootstrap            = nBootstrap,
                                  maxDistance             = maxDistance,             reverseComplement     = reverseComplement,
                                  searchPatterns          = searchPatterns,          searchRefHeader       = searchRefHeader,
                                  searchAltHeader         = searchAltHeader,         searchContextHeader   = searchContextHeader,
                                  searchIdHeader          = searchIdHeader,          searchDistanceHeader  = searchDistanceHeader,
                                  searchReverseComplement = searchReverseComplement, asTibble              = asTibble,
                                  returnEachBootstrap     = returnEachBootstrap,     searchClusterPatterns = searchClusterPatterns,
                                  no.cores                = no.cores,                linkPatterns          = linkPatterns)

  if (!searchClusterPatterns) {
    searchPatterns <- searchPatterns[nchar(dplyr::pull(searchPatterns,
                                                       searchRefHeader)) == 1, ]
    resultTable    <- resultTable[nchar(dplyr::pull(resultTable,
                                                    searchRefHeader)) == 1, ]
  }

  # Add the reverse complement to the tables if asked:
  if (searchReverseComplement) {
    revResultTable <- getRevComTable(table     = resultTable,     refHeader     = searchRefHeader,
                                     altHeader = searchAltHeader, contextHeader = searchContextHeader,
                                     idHeader  = searchIdHeader,  renameReverse = renameReverse)

    resultTable <- rbind(resultTable, revResultTable)

    searchPatterns <- dplyr::bind_rows(searchPatterns,
                                       getRevComTable(table     = searchPatterns,  refHeader     = searchRefHeader,
                                                      altHeader = searchAltHeader, contextHeader = searchContextHeader,
                                                      idHeader  = searchIdHeader,  renameReverse = renameReverse))
  }

  # Convert the sent data to the correct class and make sure no factors columns are present:
  dataTable <- convertFactor(as.data.frame(dataTable))

  # Adjust resultTable for the createSummaryPatterns function
  resultTable[nrow(resultTable) + 1, searchIdHeader] <- "Unidentified"
  resultTable <- tibble::tibble(!!rlang::sym(searchIdHeader) := unique(resultTable[ ,searchIdHeader]))
  # Add a frequence column to fill up during bootstrapping
  resultTable <- dplyr::mutate(resultTable,
                               frequency = rep.int(0, nrow(resultTable)))

  # Perform the bootstrap ---------------------------------------------------
  resultTables <- shuffleParallel(dataTable       = dataTable,       chromHeader             = chromHeader,
                                  positionHeader  = positionHeader,  refHeader               = refHeader,
                                  altHeader       = altHeader,       contextHeader           = contextHeader,
                                  sampleIdHeader  = sampleIdHeader,  nBootstrap              = nBootstrap,
                                  maxDistance     = maxDistance,     reverseComplement       = reverseComplement,
                                  searchPatterns  = searchPatterns,  searchRefHeader         = searchRefHeader,
                                  searchAltHeader = searchAltHeader, searchContextHeader     = searchContextHeader,
                                  searchIdHeader  = searchIdHeader,  searchDistanceHeader    = searchDistanceHeader,
                                  resultTable     = resultTable,     searchClusterPatterns   = searchClusterPatterns,
                                  no.cores        = no.cores,        searchReverseComplement = searchReverseComplement,
                                  linkPatterns    = linkPatterns)


  # Prepare to return the results -------------------------------------------
  # Return the raw results if asked
  if (returnEachBootstrap) {
    return(resultTables)
  }

  # Summarize the results and add comment about the summary:
  results <- summarizeBootstrap(resultTable    = resultTable,
                                resultTables   = resultTables,
                                searchIdHeader = searchIdHeader)
  results <- addShuffleMutationComment(results, searchIdHeader)

  # Return the results in the desired class:
  if (asTibble) {
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
createShuffleTable <-  function(dataTable,      chromHeader,
                                positionHeader, refHeader,
                                altHeader,      contextHeader,
                                sampleIdHeader
                                ){

  # Extract the muation data from the table --------------------------------------
  ref         <- dplyr::pull(dataTable, refHeader)
  alt         <- dplyr::pull(dataTable, altHeader)
  surrounding <- dplyr::pull(dataTable, contextHeader)

  # Create a table with random picked reference and surrounding nucleotides -------
  shuffleTable <- data.table::data.table(ref = sample(ref,
                                                      length(ref)),
                                         stringsAsFactors = F) %>%
    dplyr::mutate(surrounding = sample(surrounding,
                                       length(surrounding)))

  # A Boolean variable to enter the while loop:
  again <- TRUE
  while (again) {
    # Stays FALSE until to correct table is produced
    again <- FALSE

    # Vector to fill with random alternative nucleotides
    shuffleAlt <- c()

    # The alternative nucleotides to choose from
    alt <- dplyr::pull(dataTable, altHeader)

    # A loop to choose a random alternative nucleotide
    #  while making sure that it's not the same as the reference:
    for (index in 1:length(ref)) {

      # Temporary copy of the alternative nucleotides to choose from
      grabbag <- alt

      # The reference nucleotide that should not match with the picked
      #   alternative nucleotide:
      refNuc  <- shuffleTable[index, 1]

      # Remove all nucleotides from the grabbag that match with the reference
      grabbag <- setdiff(grabbag, c(refNuc))

      # If grabbag is empty (which should only be possible if the
      #   remaining nucleotides are all a match with the reference
      #   nucleotide) look at previous rows to see if it is
      #   possible to switch the picked alternative nucleotides:
      if (length(grabbag) == 0) {

        # Variable with the previous index nr.
        indexPrev <- index - 1

        # Loop over the build table till a candidate
        #   is found or till the end of the table:
        failed <- TRUE
        while (indexPrev != 0 & failed) {
          # Data of the candidate row:
          refNucPrev <- shuffleTable[indexPrev,1]
          altPrev    <- shuffleAlt[indexPrev]

          # Check if the candidate can be used for the switch:
          if (refNuc  != refNucPrev &
              altPrev != refNuc) {
            shuffleAlt[indexPrev] <- refNuc
            shuffleAlt[index]     <- altPrev
            alt                   <- alt[-match(refNuc,alt)]
            failed                <- FALSE
          }
          indexPrev <- indexPrev - 1

        }

        # If no candidate is found, remake the whole table.
        if(failed){
          again <- TRUE
        }

    # When there are nucleotides left in the grabbag:
    #   Grab one and remove this from the alternative
    #   nucleotide list and put them in the shuffle list.
      } else{
        random <- as.character(sample(grabbag, 1)[[1]])
        alt    <- alt[-match(random, alt)]
        shuffleAlt[length(shuffleAlt) + 1] <- random
      }
    }
  }

  # Fill the table with mutation information ----------------------------------
  shuffleTable <- shuffleTable %>%
    dplyr::mutate(sampleIDs = dplyr::pull(dataTable, sampleIdHeader),
                  chrom     = dplyr::pull(dataTable, chromHeader),
                  pos       = dplyr::pull(dataTable, positionHeader),
                  alt       = shuffleAlt,
                  check     = !alt == ref)

  # Return the shuffled table:
  return(shuffleTable)
}


#' checkParametersShuffleMutations
#' @description A function to check if the parameters are correct of the
#'   \code{\link{shuffleMutations}} function
#' @inheritParams shuffleMutations
#' @import doRNG
checkParametersShuffleMutations <- function(dataTable,               chromHeader,
                                            positionHeader,          refHeader,
                                            altHeader,               contextHeader,
                                            sampleIdHeader,          nBootstrap,
                                            maxDistance,             reverseComplement,
                                            searchPatterns,          searchRefHeader,
                                            searchAltHeader,         searchContextHeader,
                                            searchIdHeader,          searchDistanceHeader,
                                            searchReverseComplement, asTibble,
                                            returnEachBootstrap,     searchClusterPatterns,
                                            no.cores,                linkPatterns){


  # Check headers -----------------------------------------------------------
  sapply(c(chromHeader,
           positionHeader,
           refHeader,
           altHeader,
           contextHeader),
         function(x) {
           if (!any(x == names(dataTable))) {
             stop ("Please check if the column names in the parameters match with the dataTable columns.")
           }
         })
  if(nrow(dataTable[is.na(dplyr::pull(dataTable,chromHeader)) &
                    is.na(dplyr::pull(dataTable,refHeader)) &
                    is.na(dplyr::pull(dataTable,altHeader)) &
                    is.na(dplyr::pull(dataTable,contextHeader)) &
                    is.na(dplyr::pull(dataTable,positionHeader)) &
                    is.na(dplyr::pull(dataTable,sampleIdHeader)) ,]) > 0){
    stop("Error: NA's were found in the sent data.")
  }
  sapply(c(searchAltHeader,
           searchRefHeader,
           searchContextHeader,
           searchIdHeader),
         function(x) {
           if (!any(x == names(searchPatterns))) {
             stop ("Error: Please check if the column names in the parameters match with the searchPatterns table columns.")
           }
         })

  # Check numeric values ------------------------------------
  stopifnot(is.numeric(maxDistance))
  stopifnot(no.cores > 0)

  if (no.cores > parallel::detectCores()) {
    stop ("Error: The no.cores parameter is higher than amount of clusters present.")
  }

  # Check Boolean values --------------------------------------
  stopifnot(is.logical(searchClusterPatterns) &
            is.logical(reverseComplement) &
            is.logical(searchReverseComplement))
  if(!searchClusterPatterns & !linkPatterns){
    stop ("Error: Both linkPatterns and searchClusterPatterns are FALSE.
          So no patterns will be searched and no results can be made.")
  }
}


#' shuffleParallel
#' @description Function to preform the shuffling and annotating in parallel
#' @inheritParams shuffleMutations
#' @param resultTable A dataframe were the results of each bootstrap is put into
shuffleParallel <- function(dataTable,               chromHeader,
                            positionHeader,          refHeader,
                            altHeader,               contextHeader,
                            sampleIdHeader,          nBootstrap,
                            maxDistance,             reverseComplement,
                            searchPatterns,          searchRefHeader,
                            searchAltHeader,         searchContextHeader,
                            searchIdHeader,          searchDistanceHeader,
                            searchReverseComplement, searchClusterPatterns,
                            no.cores,                resultTable,
                            linkPatterns){

  # Prepare for parallel loop -----------------------------------------------
  clusters <- parallel::makeCluster(no.cores)
  doParallel::registerDoParallel(clusters)

  # Preform bootstrap --------------------------------------------------------
  resultTables <- foreach::foreach(i = iterators::icount(nBootstrap),
                                   .verbose      = TRUE,
                                   .packages     = c("magrittr", "cMut"),
                                   .export       = c("createShuffleTable")) %dorng% {

                                     # Create a table with shuffled mutations and contexts ------------------------
                                     shuffleTable <- createShuffleTable(dataTable      = dataTable,      chromHeader   = chromHeader,
                                                                        positionHeader = positionHeader, refHeader     = refHeader,
                                                                        altHeader      = altHeader,      contextHeader = contextHeader,
                                                                        sampleIdHeader = sampleIdHeader)

                                     # Identify, annotate and group clustered mutations --------------------------
                                     clusterTablePerMut <- identifyClusters(dataTable               = shuffleTable,    maxDistance          = maxDistance,
                                                                            positionHeader          = "pos",           linkPatterns         = linkPatterns,
                                                                            searchPatterns          = searchPatterns,  searchRefHeader      = searchRefHeader,
                                                                            searchAltHeader         = searchAltHeader, searchContextHeader  = searchContextHeader,
                                                                            searchIdHeader          = searchIdHeader,  searchDistanceHeader = searchDistanceHeader,
                                                                            searchReverseComplement = FALSE)

                                     clusterTable <- groupClusters(dataTable               = clusterTablePerMut, patternIntersect      = TRUE,
                                                                   showWarning             = FALSE,              searchClusterPatterns = searchClusterPatterns,
                                                                   searchRefHeader         = searchRefHeader,    searchAltHeader       = searchAltHeader,
                                                                   searchIdHeader          = searchIdHeader,     searchDistanceHeader  = searchDistanceHeader,
                                                                   searchReverseComplement = FALSE,              searchPatterns        = searchPatterns[nchar(
                                                                                                                                                          dplyr::pull(searchPatterns,
                                                                                                                                                                      searchRefHeader)) != 1, ]
                                                                   )


                                     # Add the frequencies of patterns to the resultTable -----------------------
                                     subResultTable <- resultTable
                                     subResultTable <- createSummaryPatterns(clusterTable,
                                                                             searchPatterns = subResultTable,
                                                                             searchIdHeader)

                                     total <- list("total",
                                                   nrow(clusterTablePerMut[clusterTablePerMut$is.clustered == TRUE,]))
                                     subResultTable <- rbind(subResultTable,
                                                             total)
                                   }
  parallel::stopCluster(clusters)
  return(resultTables)
}


#' summarizeBootstrap
#' @description A function to summarize the bootstrap results
#' @inheritParams shuffleMutations
#' @inheritParams shuffleParallel
#' @param resultTables List with all the results of the bootstrap
#' @importFrom rlang :=
summarizeBootstrap <- function(resultTable,  resultTables,
                               searchIdHeader) {

  # Prepare for calculating the frequencies ---------------------------------
  total       <- list("total",0)
  resultTable <- rbind(resultTable,
                       total)

  # Combine bootstrap results:
  for(table in resultTables){
    resultTable <- data.table::as.data.table(tibble::tibble(!!rlang::sym(searchIdHeader) := dplyr::pull(table,
                                                                                                        !!rlang::sym(searchIdHeader)),
                                                            frequency = resultTable$frequency + table$frequency))
  }

  total <- resultTable[dplyr::pull(resultTable,
                                   searchIdHeader) == "total", 2][[1]]
  resultTable <- resultTable[dplyr::pull(resultTable,
                                         searchIdHeader) != "total", ]

  # Calculate the percentage of the frequecies -----------------------------------
  if(total != 0){
    return(dplyr::mutate(resultTable,
                         percentage = purrr::map_dbl(.data$frequency,
                                                     function(x) {
                                                       x / total * 100
                                                     })))
  } else {
    return(dplyr::mutate(resultTable,
                         percentage = rep.int(0,
                                              nrow(resultTable))))
  }
}


#' addShuffleMutationComment
#' @description A function to add comment about the result table.
#' @param results The results from the shuffleMutation function
#' @param searchIdHeader Name of the the column with mutation processes
addShuffleMutationComment <- function(results, searchIdHeader) {
  comment(results) <-
    paste0("
    Information about the summary table columns:
    ",searchIdHeader,"     : Column with the pattern IDs from
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
  ")
  return(results)
}

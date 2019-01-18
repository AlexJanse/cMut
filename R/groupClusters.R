#' groupClusters
#' @description A function that will group the clusters and if wanted find the
#'   intersection of patterns between the mutations within a cluster. And is
#'   also able to search for cluster patterns.
#' @param dataTable A table with columns containing cluster IDs, reference and
#'   alternative nucleotide. See the output of the
#'   \code{\link{identifyClusters}} function for more information
#'   about the table.
#' @param clusterIdHeader Contains the name of the column with the cluster IDs.
#' @param refHeader Contains the name of the column with the reference
#'   nucleotides.
#' @param altHeader Contains the name of the column with the alternative
#'   nucleotides.
#' @param patternIntersect A Boolean if the table contains patterns and these
#'   needed to be processed as well.
#' @param patternHeader A string with the column name of the column with the
#'   found patterns from the \code{\link{identifyClusters}}. Only in
#'   use when patternIntersect is TRUE.
#' @param showWarning A Boolean if there need to be a warning if nrow is 0.
#' @param searchClusterPatterns A Boolean if it's needed to search to cluster
#'   patterns (e.g. GA > TT).
#' @inheritParams identifyClusters
#' @inheritParams linkPatterns
#' @param reverseComplement A Boolean to tell if the \code{ref}, \code{alt}
#'   needs to be searched with the reverse complement. Irrelevant if
#'   \code{searchClusterPatterns = FALSE} or \code{searchReverseComplement =
#'   TRUE}.
#' @export
#' @import magrittr
#' @import foreach
#' @seealso See \code{\link{mutationPatterns}} help page for a full explanation
#'   of the differences between mutation patterns and cluster patterns.
#' @examples
#' # Example of a table containing the right columns and data for the
#' # identifiAndAnnotateClusters function:
#' test <- testDataSet
#'
#' # Example of using this function with data that contain patterns:
#' mutations <- identifyClusters(dataTable      = test,
#'                               maxDistance    = 20000,
#'                               chromHeader    = "chrom",
#'                               sampleIdHeader = "sampleIDs",
#'                               positionHeader = "start",
#'                               linkPatterns   = TRUE)
#' clusters <- groupClusters(dataTable             = mutations,
#'                           patternIntersect      = TRUE,
#'                           searchClusterPatterns = FALSE)
#' # searchClusterPatterns = FALSE to emphasise the effect for the patterns from the identify function
#'
#'
#' # Example of using this function when it is needed to search for
#' # cluster patterns. Use ?mutationPatterns to learn about the
#' # difference between mutation patterns and cluster patterns.
#' clusters <- groupClusters(dataTable             = mutations,
#'                           patternIntersect      = TRUE,
#'                           searchClusterPatterns = TRUE)
#'
#' # For more information about the table:
#' cat(comment(clusters))
groupClusters <- function(dataTable,                               clusterIdHeader      = "clusterId",
                          refHeader               = "ref",         altHeader            = "alt",
                          contextHeader           = "surrounding", mutationSymbol       = ".",
                          asTibble                = TRUE,          patternIntersect     = TRUE,
                          searchClusterPatterns   = TRUE,          patternHeader        = "linkedPatterns",
                          showWarning             = TRUE,          searchPatterns       = NULL,
                          searchRefHeader         = "ref",         searchAltHeader      = "alt",
                          searchIdHeader          = "process",     searchDistanceHeader = "maxDistance",
                          searchReverseComplement = TRUE,          renameReverse        = FALSE,
                          reverseComplement       = FALSE) {


  # Check data --------------------------------------------------------------
  stopifnot(nrow(dataTable) > 0)
  if (length(setdiff(c(clusterIdHeader, refHeader, altHeader),
                     names(dataTable))) > 0){
    stop ("Error: Please check if the header parameters match with
          the column names of the sent table.")
  }
  if (!is.null(searchPatterns)) {
    if (length(setdiff(c(searchAltHeader, searchRefHeader, searchIdHeader),
                       names(searchPatterns))) > 0){
      stop("Error: Please check if the search header parameters match with
          the column names of the sent searchPatterns table.")
    }
  }
  stopifnot(is.logical(c(asTibble,         searchClusterPatterns,
                         patternIntersect, showWarning,
                         renameReverse,    searchReverseComplement)))

  # Build table -------------------------------------------------------------
  table <- createGroupTable(table         = dataTable,     clusterIdHeader = clusterIdHeader,
                            refHeader     = refHeader,     altHeader       = altHeader,
                            contextHeader = contextHeader, mutationSymbol  = mutationSymbol,
                            showWarning   = showWarning)

  # Find the pattern intersect if asked -------------------------------------
  if (patternIntersect) {
    table <- dplyr::mutate(table, foundPatterns = purrr::map(.data$cMuts,
                                                             getPatternIntersect,
                                                             patternHeader))

    # Add column to confirm if patterns are found during the intersection:
    table <- dplyr::mutate(table, has.intersect = purrr::map_lgl(.data$foundPatterns,
                                                                 function(x) {
                                                                   length(x) > 0 & x[[1]] != ""
                                                                   }))
  }

  # Search for cluster patterns if asked ------------------------------------
  if (searchClusterPatterns) {

    table <- callSearchclusterPatterns(patternIntersect = patternIntersect, searchReverseComplement = searchReverseComplement,
                                       searchPatterns   = searchPatterns,   searchRefHeader         = searchRefHeader,
                                       searchIdHeader   = searchIdHeader,   searchAltHeader         = searchAltHeader,
                                       renameReverse    = renameReverse,    searchDistanceHeader    = searchDistanceHeader,
                                       table = table)
  }

  table <- addClusterTableComment(table                 = table,
                                  patternIntersect      = patternIntersect,
                                  clusterIdHeader       = clusterIdHeader,
                                  searchClusterPatterns = searchClusterPatterns)


  # Return the results in the class of choice -------------------------------
  if (asTibble) {
    return(tibble::as.tibble(table))
  } else {
    return(as.data.frame(table))
  }

}


#' createGroupTable
#' @description A function to group the table by the cluster IDs and add
#'   annotation to it
#' @inheritParams groupClusters
#' @param table A table with columns containing cluster IDs, reference and
#'   alternative nucleotide. See the output of the
#'   \code{\link{identifyClusters}} function for more information
#'   about the table.
#' @importFrom rlang .data
createGroupTable <- function(table,         clusterIdHeader,
                             refHeader,     altHeader,
                             contextHeader, mutationSymbol,
                             showWarning) {

  # Build table -------------------------------------------------------------
  # Make sure that there are no factor columns:
  table <- convertFactor(table)

  # Group the data by the cluster ID:
  table <- dplyr::group_by_(table, clusterIdHeader)

  # Create subtables with the mutations per cluster and store them in the cMuts column:
  table <- tidyr::nest(table, .key = "cMuts")

  # Remove the row with non-clustered mutations:
  table <- dplyr::filter(table, !!rlang::sym(clusterIdHeader) != "")

  # Extract the reference, alternative and surrounding nucleotides from the subtables:
  table <- dplyr::mutate(table,
                         refs         = purrr::map(.data$cMuts,
                                                   ~as.character(dplyr::pull(.,
                                                                             refHeader))),
                         alts         = purrr::map(.data$cMuts,
                                                   ~as.character(dplyr::pull(.,
                                                                             altHeader))),
                         surroundings = purrr::map(.data$cMuts,
                                                   ~as.character(dplyr::pull(.,
                                                                             contextHeader))))

  # Collapse or fuse the data to go from a vector to a single string
  table <- dplyr::mutate(table,
                         refs = purrr::map_chr(.data$refs,
                                               function(x) {
                                                 paste0(x, collapse = "")
                                               }),
                         alts = purrr::map_chr(.data$alts,
                                               function(x) {
                                                 paste0(x,collapse = "")
                                               }))
  table <- dplyr::mutate(table,
                         surroundings = purrr::map_chr(.data$surroundings,
                                                       fuseSurroundings,
                                                       mutationSymbol,
                                                       showWarning))

  # Create a column with the mutation type of all mutations within the cluster:
  table <- dplyr::mutate(table,
                         plusStrand = purrr::map2_chr(.data$refs,
                                                      .data$alts,
                                                      formatClusterMutations),
                         minusStrand = purrr::map2_chr(.data$refs,
                                                       .data$alts,
                                                       formatClusterMutations,
                                                       convert = TRUE))
  # Create a column with the normalised mutation type:
  table <- dplyr::mutate(table,
                         clusterType = purrr::map2_chr(.data$plusStrand,
                                                       .data$minusStrand,
                                                       getClusterType))

  # Create a column with all distances within the cluster
  #  (necessary for the searchClusterPatterns function):
  table <- dplyr::mutate(table, distance = purrr::map(.data$cMuts,
                                                      function(x) {
                                                        list(x$distance)
                                                      }))
  table <- dplyr::mutate(table, distance = purrr::map(.data$distance,
                                                      function(x) {
                                                        x[[1]]
                                                      }))
  # Let know if no rows are found -------------------------------------------
  if (showWarning) {
    if (nrow(table) == 0) {
      warning ("No rows found. Please make sure the cluster IDs are present and try again.")
    }
  }

  return(table)
}


#' fuseSurroundings
#' @description A function to combine surroundings by using nucleotide symbols.
#'   See \code{\link{dnaAlphabet}} for an explanation what these symbols are.
#' @inheritParams groupClusters
#' @param surroundings The list with the surrounding nucleotides found in a
#'   cluster
fuseSurroundings <- function(surroundings, mutationSymbol,showWarning) {

  # Check surroundings length:
  if (showWarning & length(unique(sapply(surroundings, length))) > 1){
    warning ("Warning: The column with the context/surrounding nucleotides varies in size.
             The surrounding column of the result table is therefore not reliable.
             Since this is only an annotation it will not have consequences for the cMut functions.
             The affected clusters will contain one or more X's in the surroundings column.")
  }

  # Create empty vector to be filled:
  fused <- c()

  # For each position get the corresponding
  #  nucleotide of all the surroundings and
  #  find the nucleotide symbol that match them all:
  for (index in 1:nchar(surroundings[[1]])) {
    nucPos <- c()
    for (sur in surroundings) {
      nucPos <- c(nucPos,substr(sur,index,index))
    }
    fused <- c(fused, nucToSymbol(nucPos,mutationSymbol))
  }

  return(paste0(fused,collapse = ""))
}


#' nucToSymbol
#' @description Function to convert a vector nucleotides to a single symbol
#' @param nuc vector with the nucleotides
#' @param mutationSymbol symbol that represent the mutation site
nucToSymbol <- function(nuc, mutationSymbol){

  # Check if the nucleotides are not only the mutation symbol:
  if (length(setdiff(nuc,mutationSymbol)) == 0) {
    return(mutationSymbol)
  }

  # Find the symbol that match with all the nucleotides:
  for (rowIndex in 1:nrow(cMut::dnaAlphabet)) {
    if (length(setdiff(nuc,
                       strsplit(cMut::dnaAlphabet[rowIndex,]$represent,",")[[1]])) == 0) {
      return(cMut::dnaAlphabet[rowIndex,]$symbol)
    }
  }

  # If no match is found it will return an X:
  return("X")
}

#' formatClusterMutations
#' @description Function to make a mutation description symbol (e.g. G>C).
#' @param refs A list containing the reference nucleotides.
#' @param alts A list containing the alternative nucleotides.
#' @param convert A boolean that tells if the nucleotides need to be converted
#'   to the reverse complement.
formatClusterMutations <- function(refs, alts, convert=FALSE) {

  # Check parameters:
  stopifnot(length(refs) == length(alts))
  stopifnot(is.logical(convert))

  # Create the description symbols per mutation:
  if (convert) {
    seq <- paste0(alts,".",refs)
    rev <- strsplit(as.character(
                      Biostrings::reverseComplement(
                        Biostrings::DNAString(seq))),
                    "\\.")[[1]]

    result <- paste(rev[1],rev[2],sep = ">")
  } else {
    result <- paste0(refs, ">", alts)
  }


  # Combine the mutation pairs:
  finalResult <- dplyr::if_else(convert,
                         do.call(paste, as.list(rev(result))),
                         do.call(paste, as.list(result)))
  return(finalResult)
}

#' getClusterType
#' @description Function to sort and create a cluster type (e.g. C>G C>G / G>C
#'   G>C).
#' @param plusStrand A string with the plusStrand mutation description.
#' @param minusStrand A string with the minusStrand mutation description.
getClusterType <- function(plusStrand, minusStrand) {
  s <- sort(c(plusStrand, minusStrand))
  stringr::str_glue("{s[1]} / {s[2]}")
}

#' getPatternIntersect
#' @description A function to find the intersect of patterns between mutations.
#' @param clusterList A tibble with the mutation information.
#' @param patternHeader A string with the column header of the patterns.
getPatternIntersect <- function(clusterList,patternHeader){

  clusterList <- as.data.frame(clusterList)
  patterns <- c()

  # Find the intersection of the found patterns between the mutations:
  for (index in 1:nrow(clusterList)) {
   if (index == 1) {
     patterns <- clusterList[index,patternHeader][[1]]
   }
   patterns <- intersect(patterns, clusterList[index,patternHeader][[1]])
  }

  if (length(patterns) == 0) {
    return(c(""))
  }

  return(patterns)
}

#' callSearchclusterPatterns
#' @description A function to add and change columns of the groupClusters
#'   function with the results of the \code{\link{searchClusterPatterns}}
#'   function.
#' @inheritParams groupClusters
#' @inheritParams createGroupTable
callSearchclusterPatterns <- function(table,                patternIntersect,
                                      searchPatterns,       searchRefHeader,
                                      searchAltHeader,      searchIdHeader,
                                      searchDistanceHeader, searchReverseComplement,
                                      renameReverse) {

  # Create a foundPattern column if needed:
  if (!patternIntersect) {
    table <- dplyr::mutate(table,
                           foundPatterns = c(""))
  }

  # If needed get the default pattern table and/or add the reverse complement:
  if (is.null(searchPatterns)) {
    searchPatterns <- getSearchPatterns(searchReverseComplement,
                                        renameReverse = renameReverse)
  } else if (searchReverseComplement) {
    searchPatterns <- dplyr::bind_rows(searchPatterns,getRevComTable(table         = searchPatterns,
                                                                     refHeader     = searchRefHeader,
                                                                     altHeader     = searchAltHeader,
                                                                     idHeader      = searchIdHeader,
                                                                     renameReverse = renameReverse))
  }

  # Filter the search table so only cluster mutations are kept:
  searchPatterns <- searchPatterns[nchar(dplyr::pull(searchPatterns,
                                                     searchRefHeader)) > 1 |
                                     nchar(dplyr::pull(searchPatterns,
                                                       searchRefHeader)) == 0,]

  # Call the searchClusterPatterns function:
  table <- searchClusterPatterns(table,
                                 searchPatterns,
                                 searchRefHeader,
                                 searchAltHeader,
                                 searchDistanceHeader,
                                 searchIdHeader )

  return(table)
}

#' addClusterTableComment
#' @description Adds explanation to the table about the columns
#' @inheritParams groupClusters
#' @inheritParams createGroupTable
addClusterTableComment <- function(table,           patternIntersect,
                                   clusterIdHeader, searchClusterPatterns) {
  comment(table) <-
    paste0("Information about the columns:
   ",clusterIdHeader,"            : Column with the cluster ID.
   cMuts                : Column with the tables containing the
                          mutations annotation of that cluster.
   refs                 : Column with the reference nucleotides
                          of the mutations within the cluster.
   alts                 : Column with the alternative nucleotides
                          of the mutations within the cluster.
   surroundings         : Column with the combines surroundings
                          of the mutations within the cluster.
   plusStrand           : Column with the reference to variant
                          nucleotides of the mutations within the
                          cluster.
   minusStrand          : Column with the variant to reference
                          nucleotides of the mutations within the
                          cluster.
   clusterType          : The normalised version of the plusStrand,
                          minusStrand columns.
   distance             : Column with the distances found between
                          the mutations within the cluster.",
   ifelse(patternIntersect | searchClusterPatterns,
          paste0("
   foundPatterns        : Column with the patterns that are
                          connected with the cluster.",
   ifelse(patternIntersect,"
   has.intersect        : Column with a Boolean if there were
                          patterns of the mutations within the
                          cluster that overlap eachother.",""),
  ifelse(searchClusterPatterns,
   "
   has.clusterPatterns  : Column with a Boolean if there were
                          cluster patterns found.
    ","")),
    ""))
  return(table)
}

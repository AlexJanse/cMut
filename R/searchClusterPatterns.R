#' searchClusterPatterns
#' @param groupedClusters Table gained by using the groupClusters function.
#' @inheritParams groupClusters
#' @import magrittr
#' @export
searchClusterPatterns <- function(groupedClusters,
                                  searchPatterns = mutationPatterns,
                                  searchRefHeader = "ref",
                                  searchAltHeader = "alt",
                                  searchDistanceHeader = "maxDistance",
                                  searchIdHeader = "process"){

  # Check parameters -----------------------------------------------
  stopifnot(all(nchar(searchPatterns[,searchRefHeader]) > 1))
  stopifnot(all(nchar(searchPatterns[,searchAltHeader]) > 1))
  stopifnot(all(nchar(searchPatterns[,searchAltHeader]) == nchar(searchPatterns[,searchRefHeader])))
  stopifnot(!any(is.na(dplyr::select(searchPatterns,searchRefHeader,
                                     searchAltHeader,searchDistanceHeader))))
  stopifnot(!any(is.null(dplyr::select(searchPatterns,searchRefHeader,
                                       searchAltHeader,searchDistanceHeader))))
  stopifnot("cMuts" %in% names(groupedClusters))

  searchPatterns <- as.data.frame(searchPatterns)
  # Call patterns ------------------------------------------------
  clusterPatterns <- list()
  for(index in 1:nrow(groupedClusters)){
    row <- groupedClusters[index,]
    subclusterPatterns <- c()
    for(pattIndex in 1:nrow(searchPatterns)){
      id <- searchPatterns[pattIndex, searchIdHeader]
      refPat <- searchPatterns[pattIndex,searchRefHeader]
      altPat <- searchPatterns[pattIndex,searchAltHeader]
      maxDistance <- searchPatterns[pattIndex, searchDistanceHeader]
      if(grepl(refPat,row$refs) & grepl(altPat,row$alts) & row$distance <= maxDistance){
        subclusterPatterns[length(subclusterPatterns)+1] <- id
      }
    }
    if(is.null(subclusterPatterns)){
      clusterPatterns[[length(clusterPatterns)+1]] <- c("")
    } else {
      clusterPatterns[[length(clusterPatterns)+1]] <- unique(subclusterPatterns)
    }
  }
  checkList <- clusterPatterns != ""

  # Add results to table -----------------------------------------------
  groupedClusters <- dplyr::mutate(groupedClusters, clusterPatterns = clusterPatterns)
  groupedClusters <- dplyr::mutate(groupedClusters,foundPatterns = purrr::map2(foundPatterns,clusterPatterns,function(x,y){
    ifelse(x == "",list(c(y)),ifelse(y == "",list(c(x)),list(c(x,y))))}))
  groupedClusters <- dplyr::mutate(groupedClusters, foundPatterns = purrr::map(foundPatterns,function(x){x[[1]]}))
  groupedClusters <- dplyr::mutate(groupedClusters, has.clusterPatterns = checkList)
  groupedClusters$clusterPatterns <- NULL
  return(groupedClusters)
}

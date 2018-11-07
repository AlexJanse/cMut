#' getSearchPatterns
#' @description A function to get the default search pattern table
#' @param reverse A Boolean if the reverse complementend also needed to be
#'   added.
#' @param asTibble A Boolean if the returned table needs to be a tibble. Else it
#'   will sent a data.frame.
#' @param renameReverse A Boolean if the id of the process needs to be renamed.
#'   This has the effect on the cMut functions that it will no longer treat the
#'   reverse complement and non reverse complement as the same.
#' @seealso See \code{\link{mutationPatterns}} help page for full description
#'   about the columns and values.
#' @export
getSearchPatterns <- function(reverse = TRUE, renameReverse = FALSE ,asTibble = TRUE){
  if (asTibble) {
    searchPatterns <- tibble::as.tibble(mutationPatterns)
  } else {
    searchPatterns <- as.data.frame(mutationPatterns)
}

  if(reverse){
    searchPatterns <- rbind(searchPatterns,getRevComTable(searchPatterns,
                                                          refHeader = "ref",
                                                          altHeader = "alt",
                                                          contextHeader = "surrounding",
                                                          idHeader = "process",
                                                          renameReverse = renameReverse))
  }

  return(searchPatterns)
}

#' getRevComTable
#' @description A function to get the reverse complement of the known mutations
getRevComTable <- function(table, refHeader, altHeader, contextHeader = NULL, idHeader, renameReverse){
  table <- dplyr::mutate(table, !!rlang::sym(refHeader) := purrr::map_chr(!!rlang::sym(refHeader),function(x){ifelse(nchar(x) == 1,revNuc[x],as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))}))
  table <- dplyr::mutate(table, !!rlang::sym(altHeader) := purrr::map_chr(!!rlang::sym(altHeader),function(x){ifelse(nchar(x) == 1,revNuc[x],as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))}))
  if(!is.null(contextHeader)){
    table <- dplyr::mutate(table, !!rlang::sym(contextHeader) := purrr::map_chr(!!rlang::sym(contextHeader),function(x){
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
      }))
  }
  if(renameReverse){
    table <- dplyr::mutate(table, !!rlang::sym(idHeader) := paste0(!!rlang::sym(idHeader)," [Rev.Com.]"))
  } else {
    table <- dplyr::mutate(table, !!rlang::sym(idHeader) := !!rlang::sym(idHeader))
  }

  return(table)
}

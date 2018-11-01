#' getSearchPatterns
#' @description A function to get the default search pattern table
#' @param reverse A Boolean if the reverse complementend also needed to be added
#' @export
getSearchPatterns <- function(reverse = TRUE, asTibble = T){
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
                                                          idHeader = "process"))
  }

  return(searchPatterns)
}

#' getRevComTable
#' @description A function to get the reverse complement of the known mutations
getRevComTable <- function(table, refHeader, altHeader, contextHeader = NULL, idHeader){
  table <- dplyr::mutate(table, !!rlang::sym(refHeader) := purrr::map_chr(!!rlang::sym(refHeader),function(x){ifelse(nchar(x) == 1,revNuc[x],as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))}))
  table <- dplyr::mutate(table, !!rlang::sym(altHeader) := purrr::map_chr(!!rlang::sym(altHeader),function(x){ifelse(nchar(x) == 1,revNuc[x],as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))}))
  if(!is.null(contextHeader)){
    table <- dplyr::mutate(table, !!rlang::sym(contextHeader) := purrr::map_chr(!!rlang::sym(contextHeader),function(x){
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
      }))
  }
  table <- dplyr::mutate(table, !!rlang::sym(idHeader) := !!rlang::sym(idHeader))

  return(table)
}

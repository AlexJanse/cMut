#' linkPatterns
#' @description A function that is able to find matches between a submitted unkown mutation and a known mutation table
#' @param ref A string containing the reference nucleotide
#' @param alt A string containing the alternative nucleotide
#' @param context A string containing the surrounding nucleotides (e.g. G.C)
#' @param mutationSymbol A string with the symbol that stands for the mutated nucleotide location (e.g. "." in G.C)
#' @param searchPatterns A tibble with the known mutation patterns. See data/mutationPatterns.rds for an example
#' @param searchRefHeader A string with the column name of the reference nucleotide of the searchPatterns table
#' @param searchAltHeaderA string with the column name of the alternative nucleotide of the searchPatterns table
#' @param searchContextHeader A string with the column name of the context nucleotide of the searchPatterns table
#' @export
#' @import magrittr
linkPatterns <- function(ref, alt, context, mutationSymbol = ".",
                         searchPatterns = NULL, searchRefHeader = "ref",
                         searchAltHeader = "alt", searchContextHeader = "surrounding"){

  # check and adjust parameters ---------------------------------------------------------------
  stopifnot(grepl(mutationSymbol,context))
  stopifnot(is.character(ref) & is.character(alt) & is.character(context))

  if(is.null(searchPatterns)){
    searchPatterns <- tibble::as.tibble(readRDS("data/mutationPatterns.rds"))
  } else {
    # check if the assigned headers are present in the given table
    stopifnot(any(grepl(searchAltHeader,names(searchPatterns))))
    stopifnot(any(grepl(searchRefHeader,names(searchPatterns))))
    stopifnot(any(grepl(searchContextHeader,names(searchPatterns))))
  }

  dnaAlphabet <- tibble::as.tibble(readRDS("data/dnaAlphabet.rds"))
  dnaAlphabet <- dplyr::enquo(dnaAlphabet)
  ref <- dplyr::enquo(ref)
  alt <- dplyr::enquo(alt)

  # Create results -----------------------------------------------------------------
  x <- searchPatterns %>%
    dplyr::mutate(match = purrr::map_lgl(!!rlang::sym(searchRefHeader),compare,
                           getAlphaMatches(!!ref,!!dnaAlphabet))) %>%
    dplyr::filter(match == T) %>%
    dplyr::mutate(match = purrr::map_lgl(!!rlang::sym(searchAltHeader),compare,
                           getAlphaMatches(!!alt,!!dnaAlphabet))) %>%
    dplyr::filter(match == T) %>%
    dplyr::mutate(match = purrr::map_lgl(!!rlang::sym(searchContextHeader), compareContext))


  return(x)
}

compare <- function(nucleotide,symbols){
  result <- symbols[symbols$symbol == nucleotide,]

  if(nrow(result) == 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# TODO: verder werken
compareContext <- function(context1, context2, alphabet, mutationSymbol){
  context1List <- strsplit(context1,mutationSymbol)
  context1Before <- context1List[[1]][1]
  context2List <- strsplit(context2,mutationSymbol)
  context2Before <- context2List[[1]][1]
  if(nchar(context2Before) <= nchar(context1Before)){
    nNuc <- nchar(context2Before)
  } else {
    nNuc <- -1
  }

  if(nNuc == -1){
    return(FALSE)
  }

  for(index in seq.int(0,minNcharBefore)){
    if(index > 0){
      symbols <- getAlphaMatches(context1Before[nchar(context1Before)+1-index], alphabet)
      match <- compare(context2Before[nchar(context2Before)+1-index],symbols)
      if(!match){
        return(FALSE)
      }
    }
  }

}

getAlphaMatches <- function(nuc, alphabet){
  return(alphabet[grepl(nuc,alphabet$represent),1])
}


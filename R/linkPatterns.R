#' linkPatterns
#' @description A function that is able to find matches between a submitted
#'   unkown mutation and a known mutation table.
#' @param ref A string containing the reference nucleotide.
#' @param alt A string containing the alternative nucleotide.
#' @param context A string containing the surrounding nucleotides (e.g. G.C).
#' @param reverseComplement A boolean to tell if the ref, alt and context needed
#'   to be searched with the reverse complement.
#' @param mutationSymbol A string with the symbol that stands for the mutated
#'   nucleotide location (e.g. "." in G.C).
#' @param searchPatterns A tibble with the known mutation patterns. The
#'   data/mutationPatterns.rds is the default search table.
#' @param searchRefHeader A string with the column name of the reference
#'   nucleotide of the searchPatterns table.
#' @param searchAltHeader A string with the column name of the alternative
#'   nucleotide of the searchPatterns table.
#' @param searchContextHeader A string with the column name of the context
#'   nucleotide of the searchPatterns table.
#' @param searchSource A string with the column name of the ID of the known
#'   mutations.
#' @param searchReverseComplement A boolean to also search in the reverse
#'   complement of the searchPatterns tibble. When matched the returned ID will
#'   have [Rev.Com] attached to it.
#' @param patternsAsList A boolean to tell if the return value needs to be in a
#'   list or not.
#' @return list or string with the matched patterns.
#' @export
#' @import magrittr
#' @examples
#' results <- linkPatterns("C","T","TA.T",searchPatterns = mutationPatterns)
linkPatterns <- function(ref, alt, context, mutationSymbol = ".", reverseComplement = FALSE,
                         searchPatterns = mutationPatterns, searchRefHeader = "ref",
                         searchAltHeader = "alt", searchContextHeader = "surrounding",
                         searchIdHeader = "process", searchReverseComplement = TRUE,
                         patternsAsList = TRUE){

  # check and adjust parameters ---------------------------------------------------------------
  stopifnot(grepl(mutationSymbol,context))
  stopifnot(is.character(ref) & is.character(alt) & is.character(context))
  if(grepl("[^ACGTacgt]",ref)){return(list(""))}
  stopifnot(grepl("[ACGTacgt]",alt))

  # check if the assigned headers are present in the given table
  stopifnot(any(grepl(searchAltHeader,names(searchPatterns))))
  stopifnot(any(grepl(searchRefHeader,names(searchPatterns))))
  stopifnot(any(grepl(searchContextHeader,names(searchPatterns))))
  searchPatterns <- convertFactor(searchPatterns)


  # Add the reverse complement of the known table to the search table -----------------------------------------------
  if(searchReverseComplement){
    searchPatterns <- rbind(searchPatterns,getRevComTable(searchPatterns,searchRefHeader,searchAltHeader,searchContextHeader,searchIdHeader))
  }

  # Use the reverse complement of the unknown mutation ---------------------------------------------------------------
  if(reverseComplement == T){
    ref <- getReverseComplement(ref)
    alt <- getReverseComplement(alt)
    context <- getReverseComplement(context)
  }

  dnaSymbols <- tibble::as.tibble(dnaAlphabet)
  dnaSymbols <- dplyr::enquo(dnaSymbols)


  ref <- casefold(ref,upper = T)
  ref <- dplyr::enquo(ref)

  alt <- casefold(alt, upper = T)
  alt <- dplyr::enquo(alt)

  context <- casefold(context, upper = T)
  context <- dplyr::enquo(context)

  mutationSymbol <- dplyr::enquo(mutationSymbol)

  # Create results -----------------------------------------------------------------
  results <- dplyr::mutate(searchPatterns, match = purrr::map_lgl(!!rlang::sym(searchRefHeader),compare,
                           getAlphaMatches(!!ref,!!dnaSymbols)))
  results <- dplyr::filter(results, match == T)
  results <- dplyr::mutate(results, match = purrr::map_lgl(!!rlang::sym(searchAltHeader),compare,
                           getAlphaMatches(!!alt,!!dnaSymbols)))
  results <- dplyr::filter(results, match == T)
  results <- dplyr::mutate(results, match = purrr::map_lgl(!!rlang::sym(searchContextHeader), compareContext, context, mutationSymbol, dnaSymbols))
  results <- dplyr::filter(results, match == T)

  if(nrow(results) > 0){
    matchedPatterns <- list(unique(dplyr::pull(results,searchIdHeader)))
  } else {
    return(list(""))
  }
  if(!patternsAsList){
    if(length(matchedPatterns) > 0){
      matchedPatterns <- paste(matchedPatterns,collapse = "; ")
    } else {
      matchedPatterns <- list("")
    }
  } else if (length(matchedPatterns) == 0){
    matchedPatterns <- list("")
  }

  return(matchedPatterns)
}

#' compare
#' @description A function to see if the nucleotide is present in the symbol
#'   data frame
#' @param nucleotide A string with a single nucleotide
#' @param symbols A data frame with the symbols to match with the nucleotide
compare <- function(nucleotide,symbols){
  result <- symbols[symbols$symbol == nucleotide,]

  if(nrow(result) == 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' compareContext
#' @description Function to compare nucleotides context
#' @param context2 A string with a context (e.g. C.G) This one will act as the
#'   known mutation context for linkPatterns function.
#' @param context1 A string with a context (e.g.  C.C)
#' @param mutationSymbol A character that stands for the mutation nucleotide
#'   (e.g. "." for C.G)
#' @param alphabet A tibble with the nucleotides and their symbols
#' @return Boolean if the contexts match or not
compareContext <- function(context2, context1, mutationSymbol, alphabet){
  context1 <- rlang::get_expr(context1)
  context2 <- rlang::get_expr(context2)
  mutationSymbol <- rlang::get_expr(mutationSymbol)
  alphabet <- rlang::get_expr(alphabet)
  context1List <- strsplit(context1,paste0("\\",mutationSymbol))
  context1Before <- context1List[[1]][1]
  context1After <- context1List[[1]][2]

  context2List <- strsplit(context2,paste0("\\",mutationSymbol))
  context2Before <- context2List[[1]][1]
  context2After <- context2List[[1]][2]

  counter = 0
  for(context in list(c(context1Before,context2Before),c(context1After,context2After))){
    counter <- counter+1
    nNuc <- getNnuc(context[1],context[2])
    if(nNuc == -1){
      return(FALSE)
    } else if(nNuc != 0){
      for(index in seq.int(1,nNuc)){

        if(counter%%2 == 0){
          indexContext1 <- index
          indexContext2 <- index
        } else {
          indexContext1 <- nchar(context1Before)+1-index
          indexContext2 <- nchar(context2Before)+1-index
        }

        symbols <- getAlphaMatches(substr(context[1],indexContext1,indexContext1), alphabet)
        match <- compare(substr(context[2],indexContext2,indexContext2),symbols)

        if(!match){
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)

}

#' getNnuc
#' @description A function to determine how many nucleotides has to be compared
#' @return Returns a number with amount of nucleotides that has to be compared
getNnuc <- function(context1,context2){
  if(is.na(context2)){
    return(0)                   # Because the known mutation doesn't have a nucleotide so nothing can be compared
  } else if(is.na(context1)){
    return(-1)                  # Because the known mutation does have nucleotides
                                # but the unkown muation doesn't it can't be compared
                                # and will need to be treated as a mismatch because
                                # the unkown muation only have non nucleotides if the
                                # mutation has reach the end or the beginning of a
                                # chromosoom so match can't be possible
  } else if (nchar(context2) <= nchar(context1)){
    return(nchar(context2))
  } else {
    return(-1)                  # If the kown muation is bigger than the unkown
                                # muation it is not possible to succesfully compair
                                # the two as you would need to cut down the known
                                # mutation which could lead to false positives
  }
}

#' getAlphaMatches
#' @description A function the symbols that contains the nucleotide
getAlphaMatches <- function(nuc, alphabet){
  return(alphabet[grepl(nuc,alphabet$represent),1])
}


#' getRevComTable
#' @description A function to get the reverse complement of the known mutations
getRevComTable <- function(table, refHeader, altHeader, contextHeader, idHeader){
  table <- dplyr::mutate(table, !!rlang::sym(refHeader) := purrr::map_chr(!!rlang::sym(refHeader),function(x){revNuc[x]}))
  table <- dplyr::mutate(table, !!rlang::sym(altHeader) := purrr::map_chr(!!rlang::sym(altHeader),function(x){revNuc[x]}))
  table <- dplyr::mutate(table, !!rlang::sym(contextHeader) := purrr::map_chr(!!rlang::sym(contextHeader),function(x){
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))}))
  table <- dplyr::mutate(table, !!rlang::sym(idHeader) := !!rlang::sym(idHeader))

  return(table)
}

#' getSearchPatterns
#' @description A function to get the default search pattern table
#' @param reverse A Boolean if the reverse complementend also needed to be added
#' @export
getSearchPatterns <- function(reverse = TRUE){
  searchPatterns <- tibble::as.tibble(mutationPatterns)
  if(reverse){
    searchPatterns <- rbind(searchPatterns,getRevComTable(searchPatterns,
                                                          refHeader = "ref",
                                                          altHeader = "alt",
                                                          contextHeader = "surrounding",
                                                          idHeader = "process"))
  }
  return(searchPatterns)
}



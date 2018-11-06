#' linkPatterns
#' @description A function that is able to find matches between a submitted
#'   unkown mutation and a known mutation pattern table.
#' @param ref A string containing the reference nucleotide. (e.g. "C")
#' @param alt A string containing the alternative nucleotide. (e.g. "T")
#' @param context A string containing the surrounding nucleotides. (e.g. "G.C")
#' @param distance A number that tells the distance to the nearest mutation in a
#'   cluster. This parameter should stay NULL unless the \code{searchPatterns}
#'   table has a distance column and it's needed to be take into acount while
#'   comparing.
#' @param reverseComplement A Boolean to tell if the \code{ref}, \code{alt} and
#'   \code{context} needed to be searched with the reverse complement.
#' @param mutationSymbol A string with the symbol that stands for the mutated
#'   nucleotide location in the \code{context}. (e.g. "." in "G.C")
#' @param searchPatterns A tibble with the known mutation patterns. The
#'   \code{\link{mutationPatterns}} is the default search table.
#' @param searchRefHeader A string with the column name of the one with the
#'   reference nucleotide in the searchPatterns table.
#' @param searchAltHeader A string with the column name of the one with the
#'   alternative nucleotide in the searchPatterns table.
#' @param searchContextHeader A string with the column name of the one with the
#'   context nucleotide in the searchPatterns table.
#' @param searchMutationSymbol A string with symbol that stands for the mutated
#'   nucleotide location in the column of the \code{searchContextHeader}. (e.g.
#'   "." in "G.C")
#' @param searchSource A string with the column name of the ID of the known
#'   pattern.
#' @param searchDistanceHeader A string with the column name of the one with the
#'   maximum distance between clustered mutations. Not needed if the distance
#'   parameter is NULL. NA's within this column are allowed.
#' @param searchReverseComplement A boolean to also search the patterns in the
#'   reverse complement of the searchPatterns tibble.
#' @param showWarning A Boolean if warnings about the parameters are alowed.
#' @return list with the matched patterns. If nothing's found, return
#'   an empty list.
#' @export
#' @import magrittr
#' @examples
#' results <- linkPatterns("C","T","TA.T")
#'
#' # To see the default searchPattern table and it's information:
#' ?mutationPatterns
#' mutationPatterns
linkPatterns <- function(ref, alt, context, distance = NULL ,mutationSymbol = ".", reverseComplement = FALSE,
                         searchPatterns = mutationPatterns, searchRefHeader = "ref",
                         searchAltHeader = "alt", searchContextHeader = "surrounding",
                         searchIdHeader = "process", searchDistanceHeader = "maxDistance",
                         searchMutationSymbol = ".",searchReverseComplement = TRUE,
                         showWarning = TRUE){
  # check and adjust parameters ---------------------------------------------------------------
  searchPatterns <- convertFactor(searchPatterns)
  checkLinkPatternsParameters(ref, alt, context, distance, mutationSymbol, reverseComplement,
                              searchPatterns, searchRefHeader,
                              searchAltHeader, searchContextHeader,
                              searchIdHeader, searchDistanceHeader,
                              searchMutationSymbol,searchReverseComplement,
                              showWarning)

  if(grepl("[^ACGTacgt]",ref)){return(list(""))}

  # check if the assigned headers are present in the given table
  stopifnot(any(grepl(searchAltHeader,names(searchPatterns))))
  stopifnot(any(grepl(searchRefHeader,names(searchPatterns))))
  stopifnot(any(grepl(searchContextHeader,names(searchPatterns))))
  searchPatterns <- convertFactor(searchPatterns)
  searchPatterns <- searchPatterns[nchar(dplyr::pull(searchPatterns,searchRefHeader)) == 1,]

  # Add the reverse complement of the known table to the search table -----------------------------------------------
  if(searchReverseComplement){
    searchPatterns <- dplyr::bind_rows(searchPatterns,getRevComTable(table = searchPatterns,
                                                          refHeader = searchRefHeader,
                                                          altHeader = searchAltHeader,
                                                          contextHeader = searchContextHeader,
                                                          idHeader = searchIdHeader))
  }

  # Use the reverse complement of the unknown mutation ---------------------------------------------------------------
  if(reverseComplement){
    ref <- getReverseComplement(ref)
    alt <- getReverseComplement(alt)
    context <- getReverseComplement(context)
  }
  dnaSymbols <- tibble::as.tibble(dnaAlphabet)
  dnaSymbols <- dplyr::enquo(dnaSymbols)


  ref <- casefold(ref,upper = T)
  ref <- dplyr::enquo(ref)
  refSymbols <- getAlphaMatches(ref,dnaSymbols)

  alt <- casefold(alt, upper = T)
  alt <- dplyr::enquo(alt)
  altSymbols <- getAlphaMatches(alt,dnaSymbols)

  context <- casefold(context, upper = T)
  context <- dplyr::enquo(context)

  mutationSymbol <- dplyr::enquo(mutationSymbol)
  searchMutationSymbol <- as.character(searchMutationSymbol)
  searchMutationSymbol <- dplyr::enquo(searchMutationSymbol)
  # Create results -----------------------------------------------------------------
  results <- dplyr::mutate(searchPatterns, match = purrr::map_lgl(!!rlang::sym(searchRefHeader),function(x){
                                                                  compare(x,refSymbols)
                                                                  }))
  results <- dplyr::mutate(results[results$match == T,], match = purrr::map_lgl(!!rlang::sym(searchAltHeader),function(x){
                                                                                compare(x,altSymbols)
                                                                                }))

  results <- dplyr::mutate(results[results$match == T,], match = purrr::map_lgl(!!rlang::sym(searchContextHeader),
                                                                                function(x){
                                                                                compareContext(x,context, mutationSymbol, dnaSymbols, searchMutationSymbol)}
                                                                                ))
  results <- results[results$match == T,]
  if(!is.null(distance)){
    results <- dplyr::mutate(results, match = purrr::map_lgl(!!rlang::sym(searchDistanceHeader),function(x){ifelse(is.na(x),TRUE,x >= distance)}))
    results <- results[results$match == T,]
  }

  if(nrow(results) > 0){
    matchedPatterns <- list(unique(dplyr::pull(results,searchIdHeader)))
  } else {
    return(list(""))
  }
  return(matchedPatterns)
}

#' compare
#' @description A function to see if the nucleotide is present in the symbol
#'   data frame
#' @param nucleotide A string with a single nucleotide
#' @param symbols A data frame with the symbols to match with the nucleotide
compare <- function(nucleotide,symbols){
  if(nrow(symbols[symbols$symbol == nucleotide,]) == 0){
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
compareContext <- function(context2, context1, mutationSymbol, alphabet,searchMutationSymbol){
  context1 <- rlang::get_expr(context1)
  context2 <- rlang::get_expr(context2)
  mutationSymbol <- rlang::get_expr(mutationSymbol)

  alphabet <- rlang::get_expr(alphabet)
  context1List <- strsplit(context1,paste0("\\",mutationSymbol))
  context1Before <- context1List[[1]][1]
  context1After <- context1List[[1]][2]

  searchMutationSymbol <- rlang::get_expr(searchMutationSymbol)
  context2List <- strsplit(context2,paste0("\\",searchMutationSymbol))
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
  nuc <- rlang::get_expr(nuc)
  alphabet <- rlang::get_expr(alphabet)
  return(alphabet[grepl(nuc,alphabet$represent),1])
}

#' getReverseComplement
#' @description A function to get the reverse complement of a sequence
getReverseComplement <- function(x){
  return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
}

#' checkLinkPatternsParameters
#' @description A function to check the parameters of linkPatterns
checkLinkPatternsParameters <- function(ref, alt, context, distance, mutationSymbol, reverseComplement,
                                        searchPatterns, searchRefHeader,
                                        searchAltHeader, searchContextHeader,
                                        searchIdHeader, searchDistanceHeader,
                                        searchMutationSymbol,searchReverseComplement, showWarning){
  if(!all(grepl(paste0("\\",mutationSymbol),context))){
    stop("Please check if the mutationSymbol match with
         the symbol used in the context.")
  }
  if(!all(grepl(paste0("\\",searchMutationSymbol),
                dplyr::pull(
                  searchPatterns[nchar(dplyr::pull(searchPatterns,searchContextHeader)) > 0,],
                  searchContextHeader)
  )
  )
  ){
    stop("Please check if the searchMutationSymbol match with
         the symbol used in the context column of the searchPatterns table.")
  }
  stopifnot(nchar(ref) == 1 & nchar(alt) == 1)
  stopifnot(is.null(distance) | is.numeric(distance))
  stopifnot(is.character(ref) & is.character(alt) & is.character(context))
  if(showWarning){
    if(grepl(paste0("[^ACGTacgt\\",mutationSymbol,"]"), context)){
      warning(paste0("The context contain the symbol \"",stringr::str_extract(context, paste0("[^ACGTacgt\\",mutationSymbol,"]")),
              "\", this is not an A,C,G or T or match with the mutationSymbol: \"",mutationSymbol,"\".
              The results might therefore not be as expected."))
    }
  }
  stopifnot(grepl("[ACGTacgt]",alt) & grepl("[ACGTacgt]",ref))
}

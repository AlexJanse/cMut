#' getSearchPatterns
#' @description A function to get the default search pattern table
#' @param reverse A Boolean if the reverse complementend also needed to be
#'   added.
#' @param asTibble A Boolean if the returned table needs to be a tibble. Else it
#'   will sent a data.frame.
#' @param renameReverse A Boolean if the id of the process needs to be renamed.
#'   This has the effect on the cMut functions that it will no longer treat the
#'   reverse complement and non reverse complement as the same.
#' @note Please note that if there are patterns where the reverse complement is
#'   the same as the original sequence then the reverse complement won't be
#'   made. This is done to avoid double counting when \code{renameReverse =
#'   TRUE} is used.
#' @seealso See \code{\link{mutationPatterns}} help page for full description
#'   about the columns and values.
#' @export
getSearchPatterns <- function(reverse       = TRUE,
                              renameReverse = FALSE,
                              asTibble      = TRUE) {


  # Get default table in the desired class ----------------------------------
  if (asTibble) {
    searchPatterns <- tibble::as.tibble(mutationPatterns)
  } else {
    searchPatterns <- as.data.frame(mutationPatterns)
  }


  # Add the reverse complement of the default table if asked ----------------
  if (reverse) {
    searchPatterns <- rbind(searchPatterns,getRevComTable(searchPatterns,
                                                          refHeader     = "ref",
                                                          altHeader     = "alt",
                                                          contextHeader = "surrounding",
                                                          idHeader      = "process",
                                                          renameReverse = renameReverse))
  }


  # Return the table --------------------------------------------------------
  return(searchPatterns)
}

#' getRevComTable
#' @description A function to get the reverse complement of the sent pattern
#'   table
#' @inheritParams identifyClusters
#' @param table Table were the reverse complement is needed
#' @param idHeader name of the column with the mutation IDs
#' @importFrom rlang :=
getRevComTable <- function(table,
                           refHeader,
                           altHeader,
                           contextHeader = NULL,
                           idHeader,
                           renameReverse){


  # Get the reverse complement of the data ----------------------------------

  # Use the revNuc table (table with nucleotides as names and their reverse
  #   complement as values) for single nucleotides since it is faster than the
  #   reverseComplement function.
  revCom <- dplyr::mutate(table,
                          !!rlang::sym(refHeader) := purrr::map_chr(!!rlang::sym(refHeader),
                                                                    function(x){
                                                                      ifelse(nchar(x) == 1,
                                                                             revNuc[x],
                                                                             as.character(
                                                                               Biostrings::reverseComplement(
                                                                                 Biostrings::DNAString(x))))
                                                                      }))
  revCom <- dplyr::mutate(revCom, !!rlang::sym(altHeader) := purrr::map_chr(!!rlang::sym(altHeader),
                                                                            function(x){
                                                                              ifelse(nchar(x) == 1,
                                                                                     revNuc[x],
                                                                                     as.character(
                                                                                       Biostrings::reverseComplement(
                                                                                         Biostrings::DNAString(x))))
                                                                              }))

  # Get reverse complemnt of the context if available -----------------------
  if (!is.null(contextHeader)) {
    revCom <- dplyr::mutate(revCom,
                            !!rlang::sym(contextHeader) := purrr::map_chr(!!rlang::sym(contextHeader),
                                                                          function(x){
                                                                            as.character(
                                                                              Biostrings::reverseComplement(
                                                                                Biostrings::DNAString(x)))
                                                                            }))
  }


  # Add the IDs and remove duplicates with the original table ---------------
  revCom <- dplyr::mutate(revCom, !!rlang::sym(idHeader) := !!rlang::sym(idHeader))
  revCom <- setdiff.data.frame(revCom,table)


  # Rename the reverse complement IDs if necessary --------------------------
  if (renameReverse) {
    revCom <- dplyr::mutate(revCom,
                            !!rlang::sym(idHeader) := paste0(!!rlang::sym(idHeader),
                                                             " [Rev.Com.]"))
  }

  return(revCom)
}

#' setdiff.data.frame
#' @description Quick function to apply \code{\link{setdiff}} on dataframes
#' @param x data frame that is used to return the differences from y
#' @param y data frame that is used to compare with x
setdiff.data.frame <- function(x, y) {
  x[!duplicated(rbind(y, x))[-seq_len(nrow(y))], ]
}

#' Mutation patterns of mutators
#'
#' A small data set of mutationpatterns of known mutators in human DNA. If
#' you're using this as a reference to build a searchPattern table, then be sure
#' to include the * marked columns. The \*\* marked columns means that at least
#' one of the two is needed. Hereby is the data itself important and not
#' necessary the header name.
#' As seen in the table, there are two kinds of patterns:
#' \itemize{
#'   \item patterns with 1 reference, 1 alternative, with surrounding
#'   nucleotides and a maximum distance with a number or NA. These are used in
#'   the \code{\link{identifyAndAnnotateClusters}} and \code{\link{linkPatterns}}
#'   functions.
#'   \item Patterns with more than 1 or no reference, more than 1 alternative, no
#'   surrounding nucleotides and always a maximum distance. These are used in
#'   the \code{\link{groupClusters}} function. These are seperated evalued from the
#'   other patterns because these patterns depent on the order and distance of the
#'   mutations within a cluster.
#' }
#'
#' @format A data frame with 10 rows and 4 columns \describe{
#'   \item{process}{*A process that is connected with the mutation pattern.}
#'   \item{ref}{*Nucleotide that is found in the reference genome.}
#'   \item{alt}{*Nucleotide that the reference is transformed/mutated in.}
#'   \item{surrounding}{**The surrounding nucleotides where the mutation takes
#'   place. The "." represent the location of the mutation. The symbol itself
#'   can be anything as long as it's adjusted correctly when using functions
#'   that ask for the searchMutationSymbol.}
#'   \item{maxDistance}{**The maximum distance to the next mutation within a
#'   cluster. See the table to see where NA's are allowed and where not.}
#'   \item{reference}{A column with the reference of the process and it's mutation
#'   pattern.}
#'   }
#'
"mutationPatterns"

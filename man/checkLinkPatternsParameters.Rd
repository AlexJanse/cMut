% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linkPatterns.R
\name{checkLinkPatternsParameters}
\alias{checkLinkPatternsParameters}
\title{checkLinkPatternsParameters}
\usage{
checkLinkPatternsParameters(ref, alt, context, distance, mutationSymbol,
  reverseComplement, searchPatterns, searchRefHeader, searchAltHeader,
  searchContextHeader, searchIdHeader, searchDistanceHeader,
  searchMutationSymbol, searchReverseComplement, showWarning,
  renameReverse)
}
\arguments{
\item{ref}{A string containing the reference nucleotide. (e.g. "C")}

\item{alt}{A string containing the alternative nucleotide. (e.g. "T")}

\item{context}{A string containing the surrounding nucleotides. (e.g. "G.C")}

\item{distance}{A number that tells the distance to the nearest mutation in a
cluster. This parameter should stay NULL unless the \code{searchPatterns}
table has a distance column and it's needed to be take into account while
comparing.}

\item{mutationSymbol}{A string with the symbol that stands for the mutated
nucleotide location in the \code{context}. (e.g. "." in "G.C")}

\item{reverseComplement}{A Boolean to tell if the \code{ref}, \code{alt} and
\code{context} needed to be searched with the reverse complement. Irrelevant
if \code{searchReverseComplement = TRUE}.}

\item{searchPatterns}{A tibble with the known mutation patterns. The
\code{\link{mutationPatterns}} is the default search table.}

\item{searchRefHeader}{A string with the column name of the one with the
reference nucleotide in the searchPatterns table.}

\item{searchAltHeader}{A string with the column name of the one with the
alternative nucleotide in the searchPatterns table.}

\item{searchContextHeader}{A string with the column name of the one with the
context nucleotide in the searchPatterns table.}

\item{searchIdHeader}{A string with the column name of the one with the
pattern IDs.}

\item{searchDistanceHeader}{A string with the column name of the one with the
maximum distance between clustered mutations. Not needed if the distance
parameter is NULL. NA's within this column are allowed.}

\item{searchMutationSymbol}{A string with symbol that stands for the mutated
nucleotide location in the column of the \code{searchContextHeader}. (e.g.
"." in "G.C")}

\item{searchReverseComplement}{A boolean to also search the patterns in the
reverse complement of the searchPatterns tibble.}

\item{showWarning}{A Boolean if warnings about the parameters are allowed.}

\item{renameReverse}{A Boolean if the id of the process needs to be renamed.
This has the effect on the cMut functions that it will no longer treat the
reverse complement and non reverse complement as the same. This parameter
will irrelevant if \code{searchReverseComplement} is FALSE.}
}
\description{
A function to check the parameters of linkPatterns
}

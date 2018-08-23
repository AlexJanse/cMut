#' titel functie
#' omschrijving
#' @param x
#' @return x
#' @example
#' voorbeeld
#' @export
identifyAndAnnotateClusters <- function(x,
                                        chromHeader = "Chr",
                                        sampleIdHeader = "sampleID",
                                        positionHeader = "Pos") {


  # Sort data ---------------------------------------------------------------
  arrange(x, x$chromHeader, x$sampleIdHeader, x$positionHeader)



}


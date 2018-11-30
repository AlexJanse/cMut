#' identifyClusters
#' @description A function that is able to search for clusters in DNA mutation
#'   table (as explained by argument dataTable) in a single chromosome for a single
#'   sample.
#' @inheritParams identifyAndAnnotateClusters
#' @return A atomic character vector; contains per index the chromosome nr,
#'   sample name and any cluster ID's. For example: Chr1 testSample 1
identifyClusters <- function(dataTable,
                             maxDistance,
                             chromHeader    = "Chr",
                             sampleIdHeader = "sampleID",
                             positionHeader = "Pos"){

  # Create proximal distance matrix -----------------------------------------
  proxDistMatrix <- 1 / as.matrix(dist(dataTable[ ,positionHeader]))

  proxDistMatrix[proxDistMatrix <= 1 / (maxDistance + 2)] <- 0

  # Create graph from the matrix --------------------------------------------
  graph <- igraph::graph_from_adjacency_matrix(proxDistMatrix,
                                               weighted = TRUE,
                                               mode     = "undirected",
                                               diag     = FALSE)

  # Create cluster IDs ------------------------------------------------------
  clusterNo <- igraph::clusters(graph)$membership
  clusterId <- paste(dataTable[1, chromHeader],
                     dataTable[1, sampleIdHeader],
                     clusterNo)
  clusterId[!(clusterId %in% clusterId[duplicated(clusterId)])] <- ""

  return(clusterId)
}

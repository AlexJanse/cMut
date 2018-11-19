#' identifyClusters
#' @description A function that is able to search for clusters in DNA mutation
#'   table (as explained by argument x) in a single chromosoom for a single
#'   sample.
#' @inheritParams identifyAndAnnotateClusters
#' @return A atomic character vector; contains per index the chromosoom nr,
#'   sample naam and any cluster ID's. For example: Chr1 testSample 1
identifyClusters <- function(x,
                             maxDistance,
                             chromHeader="Chr",
                             sampleIdHeader="sampleID",
                             positionHeader="Pos"){

  # Create proximal distance matrix -----------------------------------------
  proxDistMatrix <- 1/as.matrix(dist(x[,positionHeader]))

  proxDistMatrix[proxDistMatrix <= 1/(maxDistance+2)] <- 0

  # Create graph from the matrix --------------------------------------------
  graph <- igraph::graph_from_adjacency_matrix(proxDistMatrix,
                                       weighted = TRUE,
                                       mode="undirected",
                                       diag=FALSE)

  # Create cluster IDs ------------------------------------------------------
  clusterNo <- igraph::clusters(graph)$membership
  clusterId <- paste(x[1,chromHeader],
                     x[1,sampleIdHeader],
                     clusterNo)
  clusterId[!(clusterId %in% clusterId[duplicated(clusterId)])] <- ""

  return(clusterId)
}

#' Calculate Interlayer Closeness for Plant Trait Multilayer Network Nodes
#'
#' @description
#' Computes the interlayer closeness (IC) for each node in a plant trait multilayer network (PTMN).
#' Interlayer closeness measures how efficiently a node can reach nodes in other layers along the
#' shortest possible path, serving as an indicator of its capacity for cross-layer information or
#' resource transfers. Nodes with high IC can efficiently influence the entire multilayer
#' network and play a central role in rapid, coordinated system-wide responses in dynamic
#' environments, thus facilitating swift phenotypic adjustments.
#'
#' @param data A data frame containing network edge information with required columns:
#'   \code{node.from}, \code{layer.from}, \code{node.to}, \code{layer.to}
#'
#' @return A data frame with two columns:
#'   \item{node}{Character vector of node names}
#'   \item{interlayer_closeness}{Numeric vector of interlayer closeness values, sorted in descending order}
#'
#' @details
#' The interlayer closeness is calculated using the formula:
#' \deqn{IC_i = \frac{n_{inter}^i}{\sum_{\beta \neq \alpha} \sum_{j \in \beta} d_{ij}}}
#' where \eqn{d_{ij}} is the shortest path length from node \eqn{v_i} to node \eqn{v_j} in another layer
#' and \eqn{n_{inter}^i} is the total number of reachable nodes in other layers.
#'
#' This parameter is part of a comprehensive set of quantitative metrics for plant trait multilayer
#' networks (PTMNs) that facilitate the effective identification of hub traits and key cross-layer
#' functional modules. These metrics are essential for understanding the coordinated
#' adaptation of plant traits across functional levels.
#'
#' @examples
#' \dontrun{
#' data(forest_invader_tree)
#' data(forest_invader_traits)
#' traits <- forest_invader_traits[, 6:73]
#' layers <- list(
#'   shoot_dynamics = c("LeafDuration", "LeafFall50", "LeafRate_max",
#'                      "Chl_shade50", "LAgain", "FallDuration",
#'                      "LeafOut", "Chl_sun50", "EmergeDuration",
#'                      "LeafTurnover"),
#'   leaf_structure = c("PA_leaf", "Mass_leaf", "Lifespan_leaf",
#'                      "Thick_leaf", "SLA", "Lobe", "LDMC",
#'                      "Stomate_size", "Stomate_index"),
#'   leaf_metabolism = c("J_max", "Vc_max", "Asat_area", "CC_mass",
#'                       "LSP", "AQY", "CC_area", "Rd_area",
#'                       "Asat_mass", "WUE", "Rd_mass", "PNUE"),
#'   leaf_chemistry = c("N_area", "Chl_area", "DNA", "Phenolics",
#'                      "Cellulose", "N_mass", "N_litter", "Chl_ab",
#'                      "Chl_mass", "N_res", "C_litter", "C_area",
#'                      "C_mass", "Ash", "Lignin", "Solubles",
#'                      "Decomp_leaf", "Hemi"),
#'   root = c("NPP_root", "SS_root", "SRL", "RTD", "RDMC",
#'            "NSC_root", "Decomp_root", "Starch_root",
#'            "C_root", "N_root", "Lignin_root"),
#'   stem = c("Latewood_diam", "Metaxylem_diam", "Earlywood_diam",
#'            "NSC_stem", "Vessel_freq", "SS_stem", "Cond_stem",
#'            "Starch_stem")
#' )
#' graph <- PTMN(traits, layers_list = layers, method = "pearson")
#' interlayer_closeness(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @importFrom igraph graph_from_data_frame distances
#' @export
interlayer_closeness <- function(data) {
  required_cols <- c("node.from", "layer.from", "node.to", "layer.to")
  if (!all(required_cols %in% colnames(data))) {
    stop("The data should contain the following: node.from, layer.from, node.to, layer.to")
  }
  # Create network diagrams
  edges <- data.frame(from = data$node.from, to = data$node.to)
  # Get all unique nodes and their layers
  nodes_from <- data.frame(node = data$node.from, layer = data$layer.from)
  nodes_to <- data.frame(node = data$node.to, layer = data$layer.to)
  all_nodes <- unique(rbind(nodes_from, nodes_to))
  # Create igraph objects
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = all_nodes)
  # Calculate the shortest path between all pairs of nodes using the distances function
  dist_matrix <- igraph::distances(g, weights = NA)
  # Initialise the results data frame
  result <- data.frame(node = character(), interlayer_closeness = numeric(), stringsAsFactors = FALSE)
  # Calculate the Interlayer Closeness of each node
  for (i in 1:nrow(all_nodes)) {
    node_name <- all_nodes$node[i]
    node_layer <- all_nodes$layer[i]
    # Get nodes from other layers
    other_layers_indices <- which(all_nodes$layer != node_layer)
    # If there are no nodes in other layers, set the Interlayer Closeness to NA
    if (length(other_layers_indices) == 0) {
      result <- rbind(result, data.frame(node = node_name, interlayer_closeness = NA))
      next
    }
    shortest_paths <- dist_matrix[node_name, all_nodes$node[other_layers_indices]]
    # (n_i^inter)
    reachable_nodes <- sum(!is.infinite(shortest_paths), na.rm = TRUE)
    if (reachable_nodes == 0) {
      result <- rbind(result, data.frame(node = node_name, interlayer_closeness = 0))
      next
    }
    sum_distances <- sum(shortest_paths[!is.infinite(shortest_paths)], na.rm = TRUE)
    # IC_i = n_i^inter / ∑d_ij
    interlayer_closeness <- reachable_nodes / sum_distances
    result <- rbind(result, data.frame(node = node_name, interlayer_closeness = interlayer_closeness))
  }
  result <- result[order(-result$interlayer_closeness), ]
  return(result)
}

#' Calculate Average Interlayer Path Length (AIPL) for Plant Trait Multilayer Networks
#'
#' This function calculates the Average Interlayer Path Length (AIPL) for a Plant Trait
#' Multilayer Network (PTMN). AIPL describes the average path length between pairs of
#' nodes from different network layers and reflects how efficiently information or
#' resources can be transmitted across the network. Lower AIPL values indicate
#' tighter, more efficient cross-layer connections, supporting the rapid transmission
#' of environmental cues and resources among network layers and thereby enhancing the
#' overall adaptability of the system.
#'
#' @param data A data frame containing the PTMN edge list with required columns:
#'   \code{node.from}, \code{layer.from}, \code{node.to}, and \code{layer.to}.
#'   This data frame represents the network structure with both intralayer and
#'   interlayer connections.
#'
#' @return A numeric value representing the average interlayer path length across
#'   all reachable node pairs between different layers. Returns \code{NaN} if no
#'   interlayer paths exist or all interlayer node pairs are unreachable.
#'
#' @details
#' The AIPL is calculated using the formula from multilayer network theory:
#'
#' \deqn{AIPL = \frac{1}{\sum_{α=1}^{M-1} \sum_{β=α+1}^{M} \sum_{i∈α}\sum_{j∈β}1_{d_{ij} < ∞}} \sum_{α=1}^{M-1} \sum_{β=α+1}^{M} \sum_{i∈α}\sum_{j∈β} d_{ij}}
#'
#' where \eqn{d_{ij}} is the shortest path length from node i (in layer α) to node j
#' (in layer β) and \eqn{1_{d_{ij} < ∞}} is an indicator function of node reachability.
#'
#' The function:
#' \enumerate{
#'   \item Creates an undirected graph from the edge list
#'   \item Calculates shortest path distances between all node pairs
#'   \item Considers only node pairs from different layers (interlayer connections)
#'   \item Includes only finite path lengths (reachable pairs)
#'   \item Returns the mean of all valid interlayer path lengths
#' }
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
#' average_interlayer_path_length(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @importFrom igraph graph_from_data_frame distances
#' @export
average_interlayer_path_length <- function(data) {
  required_cols <- c("node.from", "layer.from", "node.to", "layer.to")
  if (!all(required_cols %in% colnames(data))) {
    stop("The data should contain the following: node.from, layer.from, node.to, layer.to")
  }
  # Create network edges
  edges <- data.frame(from = data$node.from, to = data$node.to)
  # Get all unique nodes and their layers
  nodes_from <- data.frame(node = data$node.from, layer = data$layer.from)
  nodes_to <- data.frame(node = data$node.to, layer = data$layer.to)
  all_nodes <- unique(rbind(nodes_from, nodes_to))
  # Create igraph objects
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = all_nodes)
  # Calculate the shortest path distance between all pairs of nodes
  dist_matrix <- igraph::distances(g, weights = NA)
  # Initialise the list of interlayer path lengths
  interlayer_paths <- numeric()
  # Calculate the shortest path between all pairs of nodes between layers
  for (i in 1:(nrow(all_nodes)-1)) {
    for (j in (i+1):nrow(all_nodes)) {
      node_i <- all_nodes$node[i]
      layer_i <- all_nodes$layer[i]
      node_j <- all_nodes$node[j]
      layer_j <- all_nodes$layer[j]
      # Consider only node pairs between different layers
      if (layer_i != layer_j) {
        path_length <- dist_matrix[node_i, node_j]
        # Add only a finite path length (reachable node pairs)
        if (is.finite(path_length)) {
          interlayer_paths <- c(interlayer_paths, path_length)
        }
      }
    }
  }
  # Calculate average interlayer path length
  avg_path_length <- mean(interlayer_paths, na.rm = TRUE)
  result <- avg_path_length
  return(result)
}

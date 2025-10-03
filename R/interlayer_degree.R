#' Calculate Interlayer Degree for Nodes in Plant Trait Multilayer Networks
#'
#' This function calculates the interlayer degree for each node in a Plant Trait Multilayer Network (PTMN).
#' The interlayer degree represents the number of connections a node has to nodes in different layers,
#' which is crucial for understanding cross-layer coordination and functional integration across plant organs.
#'
#' @param data A data frame containing multilayer network edge information with columns:
#'   \code{node.from}, \code{layer.from}, \code{node.to}, \code{layer.to}. This represents
#'   the extended edge list from a PTMN construction.
#'
#' @return A data frame with three columns:
#'   \item{node}{Character. The name of the trait node}
#'   \item{layer}{Character. The layer (organ or functional system) containing the node}
#'   \item{interlayer_degree}{Numeric. The number of interlayer connections for each node}
#'   The results are ordered by interlayer degree in descending order.
#'
#' @details
#' In PTMNs, interlayer edges represent interactions between traits in different layers,
#' illustrating the integration of multiple functional systems.
#' Nodes with higher interlayer degrees indicate traits that are more integrated across
#' plant organs or functional systems, suggesting coordinated adaptation mechanisms.
#'
#' The function counts both incoming and outgoing interlayer connections for each node,
#' as PTMNs are constructed as undirected networks reflecting the reciprocal nature
#' of trait relationships in plants.
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
#' interlayer_degree(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @export
interlayer_degree <- function(data) {
  required_cols <- c("node.from", "layer.from", "node.to", "layer.to")
  if(!all(required_cols %in% colnames(data))) {
    stop("The data should contain the following: node.from, layer.from, node.to, layer.to")
  }
  # Get all unique nodes and their layers
  all_nodes_from <- data.frame(node = data$node.from, layer = data$layer.from)
  all_nodes_to <- data.frame(node = data$node.to, layer = data$layer.to)
  all_nodes <- unique(rbind(all_nodes_from, all_nodes_to))
  result <- data.frame(node = character(), layer = character(), interlayer_degree = numeric(), stringsAsFactors = FALSE)
  # Calculate interlayer degree for each node
  for(i in 1:nrow(all_nodes)) {
    current_node <- all_nodes$node[i]
    current_layer <- all_nodes$layer[i]
    edges_from <- data[data$node.from == current_node & data$layer.to != current_layer, ]
    edges_to <- data[data$node.to == current_node & data$layer.from != current_layer, ]
    interlayer_degree <- nrow(edges_from) + nrow(edges_to)
    result <- rbind(result, data.frame(node = current_node,
                                       layer = current_layer,
                                       interlayer_degree = interlayer_degree))
  }
  result <- result[order(-result$interlayer_degree), ]
  return(result)
}

#' Calculate Interlayer Clustering Coefficient for Plant Trait Multilayer Networks (PTMNs)
#'
#' This function calculates the interlayer clustering coefficient (ICC) for each node in a
#' plant trait multilayer network. The ICC measures interconnection among a node's interlayer
#' neighbors, indicating the node's capacity to form cooperative, cross-layer functional
#' modules within the local network.
#'
#' @param data A data frame containing multilayer network edge information with required
#'   columns: \code{node.from}, \code{layer.from}, \code{node.to}, and \code{layer.to}.
#'   This is typically the output from the \code{\link{PTMN}} function.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{node}{Character. The node identifier (trait name)}
#'   \item{layer}{Character. The layer identifier (functional category)}
#'   \item{interlayer_degree}{Numeric. Number of interlayer connections for the node}
#'   \item{actual_connections}{Numeric. Actual number of connections between interlayer neighbors}
#'   \item{potential_connections}{Numeric. Maximum possible connections between interlayer neighbors}
#'   \item{interlayer_clustering_coefficient}{Numeric. The ICC value calculated as actual_connections / potential_connections}
#' }
#' The results are ordered by decreasing interlayer clustering coefficient values.
#'
#' @details
#' The interlayer clustering coefficient is defined as:
#' \deqn{ICC_i = \frac{2E_i}{IK_i(IK_i - 1)}}
#' where \eqn{E_i} is the number of actual edges among the node's interlayer neighbors
#' and \eqn{IK_i} is the interlayer degree.
#'
#' For nodes with an interlayer degree of zero or one, ICC is defined as zero.
#' Nodes with high ICC act as "organizers" of highly integrated, cross-functional modules,
#' enhancing system-level responsiveness and phenotypic coordination.
#'
#' The function constructs PTMNs as undirected networks, reflecting the reciprocal nature
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
#' interlayer_clustering_coefficient(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @export
interlayer_clustering_coefficient <- function(data) {
  required_cols <- c("node.from", "layer.from", "node.to", "layer.to")
  if(!all(required_cols %in% colnames(data))) {
    stop("The data should contain the following: node.from, layer.from, node.to, layer.to")
  }
  all_nodes_from <- data.frame(node = data$node.from, layer = data$layer.from)
  all_nodes_to <- data.frame(node = data$node.to, layer = data$layer.to)
  all_nodes <- unique(rbind(all_nodes_from, all_nodes_to))
  result <- data.frame(node = character(),
                       layer = character(),
                       interlayer_degree = numeric(),
                       actual_connections = numeric(),
                       potential_connections = numeric(),
                       interlayer_clustering_coefficient = numeric(),
                       stringsAsFactors = FALSE
  )
  # Calculate inter-layer clustering coefficients for each node
  for(i in 1:nrow(all_nodes)) {
    current_node <- all_nodes$node[i]
    current_layer <- all_nodes$layer[i]
    # Find out the interlayer neighbours of the current node (connections from other layers)
    edges_from <- data[data$node.from == current_node & data$layer.from == current_layer & data$layer.to != current_layer, ]
    edges_to <- data[data$node.to == current_node & data$layer.to == current_layer & data$layer.from != current_layer, ]
    # Get a list of neighbouring nodes between tiers
    neighbors_from_other_layers <- c()
    if(nrow(edges_from) > 0) {
      neighbors_from_other_layers <- c(neighbors_from_other_layers,
                                       paste(edges_from$node.to, edges_from$layer.to, sep="-"))
    }
    if(nrow(edges_to) > 0) {
      neighbors_from_other_layers <- c(neighbors_from_other_layers,
                                       paste(edges_to$node.from, edges_to$layer.from, sep="-"))
    }
    interlayer_degree <- length(neighbors_from_other_layers)
    # Clustering coefficients are defined as 0 when the degree of interstratification is 0 or 1
    if(interlayer_degree <= 1) {
      result <- rbind(result, data.frame(node = current_node,
                                         layer = current_layer,
                                         interlayer_degree = interlayer_degree,
                                         actual_connections = 0,
                                         potential_connections = 0,
                                         interlayer_clustering_coefficient = 0))
      next
    }
    # Calculate the actual number of connections between neighbours between tiers
    actual_connections <- 0
    if(interlayer_degree >= 2) {
      neighbor_info <- do.call(rbind, strsplit(neighbors_from_other_layers, "-"))
      neighbor_nodes <- neighbor_info[,1]
      neighbor_layers <- neighbor_info[,2]
      # Check that a connection exists between each pair of neighbours
      for(j in 1:(length(neighbor_nodes)-1)) {
        for(k in (j+1):length(neighbor_nodes)) {
          conn1 <- data[(data$node.from == neighbor_nodes[j] & data$layer.from == neighbor_layers[j] &
                           data$node.to == neighbor_nodes[k] & data$layer.to == neighbor_layers[k]), ]
          conn2 <- data[(data$node.from == neighbor_nodes[k] & data$layer.from == neighbor_layers[k] &
                           data$node.to == neighbor_nodes[j] & data$layer.to == neighbor_layers[j]), ]
          if(nrow(conn1) > 0 || nrow(conn2) > 0) {
            actual_connections <- actual_connections + 1
          }
        }
      }
    }
    potential_connections <- interlayer_degree * (interlayer_degree - 1) / 2
    interlayer_clustering_coefficient <- ifelse(potential_connections > 0, actual_connections / potential_connections, 0)
    result <- rbind(result, data.frame(node = current_node,
                                       layer = current_layer,
                                       interlayer_degree = interlayer_degree,
                                       actual_connections = actual_connections,
                                       potential_connections = potential_connections,
                                       interlayer_clustering_coefficient = interlayer_clustering_coefficient))
  }
  result <- result[order(-result$interlayer_clustering_coefficient), ]
  return(result)
}

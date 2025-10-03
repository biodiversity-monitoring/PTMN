#' Identify Cross-Layer Modules in Plant Trait Multilayer Networks
#'
#' This function identifies network modules that span multiple functional layers
#' in Plant Trait Multilayer Networks (PTMNs). Cross-layer modules represent
#' coordinated trait groups that integrate multiple plant organs or functional
#' systems, which are crucial for understanding multidimensional plant adaptations
#' and phenotypic integration.
#'
#' @param data A data frame containing the PTMN edge information with columns:
#'   \code{node.from}, \code{node.to}, \code{layer.from}, and \code{layer.to}.
#'   This is typically the output from the \code{PTMN} function.
#'
#' @return A named list where each element represents a cross-layer module.
#'   Each list element contains the node names (trait names) that belong to
#'   that cross-layer module. Only modules spanning multiple layers are returned.
#'   Returns an empty list if no cross-layer modules are detected.
#'
#' @details
#' The function uses the Girvan-Newman community detection algorithm to identify
#' network modules, then filters to retain only those modules that contain nodes
#' from multiple functional layers. Cross-layer modules are particularly
#' important in plant ecology as they represent coordinated trait groups that
#' facilitate integrated responses across different organs and functional systems.
#'
#' In the context of plant invasion ecology, species with more cross-layer modules
#' may exhibit enhanced phenotypic integration and coordination, potentially
#' contributing to invasion success.
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
#' cross_layer_groups(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @importFrom igraph graph_from_data_frame communities
#' @export
cross_layer_groups <- function(data) {
  gg <- data[,c("node.from","node.to")]
  g <- igraph::graph_from_data_frame(gg, directed = FALSE)
  ceb <- suppressWarnings(igraph::cluster_edge_betweenness(g, weights = NA))
  ceb_groups <- igraph::communities(ceb)
  node_layer_dict <- unique(rbind(data.frame(node = data$node.from, layer = data$layer.from),
                                  data.frame(node = data$node.to, layer = data$layer.to)))
  cross_layer_groups <- list()
  for (group_name in names(ceb_groups)) {
    group_nodes <- ceb_groups[[group_name]]
    group_layers <- unique(node_layer_dict$layer[node_layer_dict$node %in% group_nodes])
    # If the grouping corresponds to more than 1 layer, it is a cross-layer grouping
    if (length(group_layers) > 1) {
      cross_layer_groups[[group_name]] <- group_nodes
    }
  }
  return(cross_layer_groups)
}

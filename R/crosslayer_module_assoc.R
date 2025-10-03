#' Calculate Module Interlayer Association (MIA) for Plant Trait Multilayer Networks
#'
#' This function computes the Module Interlayer Association (MIA), a key topological
#' parameter for Plant Trait Multilayer Networks (PTMNs). MIA quantifies the proportion
#' of interlayer edges within each network module, reflecting the degree of cross-layer
#' integration within modules. High MIA values suggest the formation of tightly
#' connected cross-layer modules that facilitate coordinated adaptation across different
#' plant organs and functional systems.
#'
#' @param data A data frame containing the PTMN edge information with columns:
#'   \code{node.from}, \code{node.to}, \code{layer.from}, and \code{layer.to}.
#'   This is typically the output from the \code{PTMN} function.
#'
#' @return A numeric value representing the average Module Interlayer Association (MIA).
#'   Returns \code{NA} if no interlayer modules are found. Values range from 0 to 1,
#'   where higher values indicate stronger cross-layer integration within modules.
#'
#' @details
#' The function uses the Girvan-Newman algorithm to detect network modules, then calculates
#' the proportion of interlayer edges within each module that spans multiple layers.
#' The MIA is computed as the average of interlayer association values across all
#' cross-layer modules. This metric is particularly useful for understanding how plant
#' traits coordinate across different organs (roots, stems, leaves) and functional
#' systems in response to environmental pressures.
#'
#' In the context of plant invasion ecology, higher MIA values have been associated
#' with greater invasion success, as they may indicate more efficient coordination
#' among functional systems.
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
#' crosslayer_module_assoc(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @export
crosslayer_module_assoc <- function(data) {
  gg <- data[,c("node.from","node.to")]
  g <- graph_from_data_frame(gg, directed = FALSE)
  ceb <- suppressWarnings(igraph::cluster_edge_betweenness(g, weights = NA))
  ceb_groups <- igraph::communities(ceb)
  interlayer_assoc_list <- c()
  for (group_name in names(ceb_groups)) {
    group_nodes <- ceb_groups[[group_name]]
    # Determine how many layers the module's nodes come from
    node_layers_in_module <- unique(c(data$layer.from[data$node.from %in% group_nodes], data$layer.to[data$node.to %in% group_nodes]))
    if (length(node_layers_in_module) < 2) next
    # Find all edges in the module
    ind_in_group <- (data$node.from %in% group_nodes) & (data$node.to %in% group_nodes)
    module_edges <- data[ind_in_group, ]
    if (nrow(module_edges) == 0) next
    total_edges <- nrow(module_edges)
    interlayer_edges <- sum(module_edges$layer.from != module_edges$layer.to)
    assoc <- interlayer_edges / total_edges
    interlayer_assoc_list <- c(interlayer_assoc_list, assoc)
  }
  # Calculate the mean of the interlayer association for the interlayer module
  if (length(interlayer_assoc_list) == 0) {
    return(NA)
  } else {
    return(mean(interlayer_assoc_list))
  }
}

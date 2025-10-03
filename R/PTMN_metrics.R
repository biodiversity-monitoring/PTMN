#' Calculate Plant Trait Multilayer Network (PTMN) Metrics
#'
#' This function computes comprehensive network metrics for plant trait multilayer networks,
#' including both node-level and global network properties. It quantifies the structural
#' characteristics and connectivity patterns that reveal coordinated adaptation mechanisms
#' across multiple plant organs and functional systems.
#'
#' @param graph A multilayer network object created by the \code{\link{PTMN}} function,
#'   representing trait relationships across different plant organs or functional systems.
#'
#' @return A list containing two data frames:
#'   \describe{
#'     \item{node}{Node-level metrics including:
#'       \itemize{
#'         \item \code{node}: Name of the trait node
#'         \item \code{layer}: Layer (organ/functional system) containing the node
#'         \item \code{interlayer_degree}: Number of connections a node has across layers
#'         \item \code{interlayer_closeness}: Measure of how close a node is to all other nodes across layers
#'         \item \code{interlayer_clustering_coefficient}: Local clustering coefficient considering interlayer connections
#'       }
#'     }
#'     \item{global}{Global network metrics including:
#'       \itemize{
#'         \item \code{Interlayer_edge_density}: Proportion of possible interlayer connections that are realized
#'         \item \code{Average_interlayer_path_length}: Average shortest path length between nodes across layers
#'         \item \code{Average_interlayer_clustering_coefficient}: Network-wide clustering coefficient
#'         \item \code{Module_interlayer_association}: Degree of modular organization across layers
#'       }
#'     }
#'   }
#'
#' @details
#' The PTMN framework addresses limitations of traditional single-layer plant trait networks
#' by capturing trait relationships across multiple plant organs and functional systems.
#' These metrics quantify network topology using specially designed parameters that facilitate
#' identification of hub traits and key cross-layer functional modules, essential for
#' understanding coordinated adaptation of plant traits.
#'
#' Node-level metrics help identify traits that serve as integration points across different
#' plant systems, while global metrics characterize the overall network structure and
#' connectivity patterns. Higher interlayer connectivity may indicate greater functional
#' integration and potentially enhanced adaptive capacity.
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
#' PTMN_metrics(graph)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @export
PTMN_metrics <- function(graph) {
  # Node-level metrics
  Interlayer_degree <- interlayer_degree(graph)
  Interlayer_closeness <- interlayer_closeness(graph)
  Interlayer_clustering <- interlayer_clustering_coefficient(graph)
  merge1 <- merge(Interlayer_degree,Interlayer_closeness, by = "node", all = TRUE)
  merge2 <- merge(merge1, Interlayer_clustering, by = c("node", "layer", "interlayer_degree"), all = TRUE)
  node_metrics <- merge2[,c("node","layer","interlayer_degree","interlayer_closeness","interlayer_clustering_coefficient")]
  node_metrics <- node_metrics[order(-node_metrics$interlayer_degree), ]
  # Global metrics
  Interlayer_edge_density <- interlayer_edge_density(graph)
  Average_interlayer_path_length <- average_interlayer_path_length(graph)
  Average_interlayer_clustering_coefficient <- average_interlayer_clustering_coefficient(graph)
  Module_interlayer_association <- crosslayer_module_assoc(graph)
  global_metrics <- data.frame(Interlayer_edge_density = Interlayer_edge_density,
                               Average_interlayer_path_length = Average_interlayer_path_length,
                               Average_interlayer_clustering_coefficient = Average_interlayer_clustering_coefficient,
                               Module_interlayer_association = Module_interlayer_association)
  colnames(global_metrics) <- c("Interlayer_edge_density", "Average_interlayer_path_length", "Average_interlayer_clustering_coefficient", "Module_interlayer_association")
  return(list(node = node_metrics, global = global_metrics))
}

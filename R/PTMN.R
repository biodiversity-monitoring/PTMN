#' Plant Trait Multilayer Network Analysis
#'
#' @description
#' Constructs a phylogenetic trait multilayer network by calculating correlations
#' between traits while accounting for phylogenetic relationships. The function
#' creates a network where nodes represent traits organized into different layers,
#' and edges represent significant correlations between traits.
#'
#' @param traits_matrix A data frame or matrix containing trait values where rows
#'   represent species and columns represent traits. Row names should match the
#'   tip labels in the phylogenetic tree when phylogenetic correction is used.
#' @param layers_list A named list where each element contains the names of traits
#'   belonging to that layer. Names of the list elements represent layer names.
#' @param rThres Numeric. Correlation threshold for filtering weak correlations.
#'   Correlations with absolute values below this threshold are set to zero.
#'   Default is 0.2.
#' @param pThres Numeric. P-value threshold for statistical significance after
#'   FDR correction. Correlations with p-values above this threshold are excluded.
#'   Default is 0.05.
#' @param method Character. Correlation method to use. Either "pearson" or
#'   "spearman". Default is "pearson".
#' @param phylo_correction Logical. Whether to apply phylogenetic correction
#'   using phylogenetic independent contrasts. Default is FALSE.
#' @param phylo_tree A phylo object from the ape package. Required when
#'   phylo_correction = TRUE. Should contain all species present in traits_matrix.
#'
#' @return
#' A data frame containing the network edges with the following columns:
#' \itemize{
#'   \item \code{node.from}: Name of the source trait
#'   \item \code{layer.from}: Layer name of the source trait
#'   \item \code{node.to}: Name of the target trait
#'   \item \code{layer.to}: Layer name of the target trait
#'   \item \code{correlation}: Correlation coefficient between the traits
#' }
#' Returns an empty data frame with the same structure if no significant
#' correlations are found.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data consistency
#'   \item Calculates correlation matrix using standard methods or phylogenetic
#'         independent contrasts
#'   \item Filters correlations based on correlation and p-value thresholds
#'   \item Applies FDR correction to p-values
#'   \item Constructs an igraph network object
#'   \item Maps traits to their respective layers
#'   \item Returns edge list with layer information
#' }
#'
#' When phylogenetic correction is enabled, the function uses phylogenetic
#' independent contrasts to account for the non-independence of species data
#' due to shared evolutionary history.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(forest_invader_tree)
#' data(forest_invader_traits)
#' traits <- forest_invader_traits[, 6:73]
#'
#' # Define trait layers
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
#' # Run PTMN analysis without phylogenetic correction
#' PTMN(traits, layers_list = layers, method = "pearson")
#'
#' # Run PTMN analysis with phylogenetic correction
#' PTMN(traits, layers_list = layers, method = "pearson",
#'      phylo_correction = TRUE, phylo_tree = forest_invader_tree)
#'}
#'
#' @importFrom igraph graph_from_adjacency_matrix simplify delete_vertices
#'   degree E V ecount as_data_frame graph_from_data_frame
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom ape keep.tip
#'
#' @export
PTMN <- function(traits_matrix, layers_list, rThres = 0.2, pThres = 0.05,
                 method = "pearson", phylo_correction = FALSE, phylo_tree = NULL) {
  method <- match.arg(method, choices = c("pearson", "spearman"))
  # Input validation
  if (!is.data.frame(traits_matrix) && !is.matrix(traits_matrix)) {
    stop("traits_matrix must be a dataframe or matrix")
  }
  if (!is.list(layers_list)) {
    stop("layers_list must be a list where each element contains the names of the traits belonging to that layer")
  }
  layer_names <- names(layers_list)
  if (is.null(layer_names) || any(layer_names == "")) {
    stop("layers_list must have named elements indicating the name of each layer")
  }
  # Check all traits are present in traits_matrix
  all_traits <- unlist(layers_list)
  if (!all(all_traits %in% colnames(traits_matrix))) {
    missing_traits <- all_traits[!all_traits %in% colnames(traits_matrix)]
    stop(paste("The following traits are not present in traits_matrix:",
               paste(missing_traits, collapse = ", ")))
  }
  # CHECK: Ensure all traits in matrix are defined in layers
  all_traits_in_matrix <- colnames(traits_matrix)
  all_traits_in_layers <- unlist(layers_list)
  unlayered_traits <- all_traits_in_matrix[!all_traits_in_matrix %in% all_traits_in_layers]
  if (length(unlayered_traits) > 0) {
    stop(paste("The following traits are in traits_matrix but not defined in layers_list:",
               paste(unlayered_traits, collapse = ", ")))
  }
  # Phylogenetic correction parameter validation
  if (phylo_correction && is.null(phylo_tree)) {
    stop("phylo_tree must be provided when phylo_correction = TRUE")
  }
  if (phylo_correction) {
    if (!inherits(phylo_tree, "phylo")) {
      stop("phylo_tree must be a phylo object from the ape package")
    }
    # Check if species names in phylogenetic tree match traits_matrix row names
    if (!all(rownames(traits_matrix) %in% phylo_tree$tip.label)) {
      missing_species <- rownames(traits_matrix)[!rownames(traits_matrix) %in% phylo_tree$tip.label]
      stop(paste("The following species are not present in phylo_tree:",
                 paste(missing_species, collapse = ", ")))
    }
    # Ensure data and phylogenetic tree species order are consistent
    common_species <- intersect(rownames(traits_matrix), phylo_tree$tip.label)
    traits_matrix <- traits_matrix[common_species, ]
    phylo_tree <- ape::keep.tip(phylo_tree, common_species)
  }
  # Calculate correlation matrix
  if (phylo_correction) {
    # Calculate phylogenetic independent contrasts correlation
    correlation_matrix <- phylo_correlation(traits_matrix, phylo_tree, method = method)
  } else {
    # Standard correlation analysis
    correlation_matrix <- Hmisc::rcorr(as.matrix(traits_matrix), type = method)
  }
  r <- correlation_matrix$r
  r[abs(r) < rThres] <- 0
  p <- correlation_matrix$P
  p <- stats::p.adjust(p, method = 'fdr')
  p[p >= pThres] <- -1
  p[p < pThres & p >= 0] <- 1
  p[p == -1] <- 0
  z <- r * p
  diag(z) <- 0
  g <- igraph::graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')
  g <- igraph::simplify(g)
  g <- igraph::delete_vertices(g, names(igraph::degree(g)[igraph::degree(g) == 0]))
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  gdf <- igraph::as_data_frame(g, what = "edges")
  gdf <- gdf[,c("from","to","correlation")]
  nodes_name <- unique(c(gdf$from, gdf$to))
  nodes <- data.frame(id = nodes_name, num = 1)
  g <- igraph::graph_from_data_frame(gdf, vertices = nodes, directed = FALSE)
  if (igraph::ecount(g) == 0) {
    return(data.frame(
      node.from = character(0),
      layer.from = character(0),
      node.to = character(0),
      layer.to = character(0),
      correlation = numeric(0)
    ))
  }
  edges_df <- igraph::as_data_frame(g, what = "edges")
  # Create node to layer mapping
  node_to_layer <- list()
  for (layer_name in layer_names) {
    for (trait in layers_list[[layer_name]]) {
      node_to_layer[[trait]] <- layer_name
    }
  }
  # Create extended edges dataframe with layer information
  extended_edges <- data.frame(
    node.from = edges_df$from,
    layer.from = sapply(edges_df$from, function(node) node_to_layer[[node]]),
    node.to = edges_df$to,
    layer.to = sapply(edges_df$to, function(node) node_to_layer[[node]]),
    correlation = edges_df$correlation,
    stringsAsFactors = FALSE
  )
  return(extended_edges)
}

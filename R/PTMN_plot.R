#' Plot Plant Trait Multilayer Networks (PTMN)
#'
#' @description
#' This function creates visualizations of plant trait multilayer networks (PTMNs) constructed using the PTMN framework.
#' PTMNs systematically integrate multilayer network theory with plant functional trait analysis, enabling quantitative
#' assessments of trait relationships across plant organs and functional systems.
#'
#' @param data A data frame containing the PTMN edge list with columns: node.from, node.to, layer.from, layer.to.
#'   This should be the output from the \code{PTMN} function.
#' @param style Character string specifying the layout style. Options are:
#'   \itemize{
#'     \item "default" - Standard network layout with cross-layer module highlighting
#'     \item "circle" - Circular layout with nodes arranged by layer
#'   }
#' @param vertex_size Numeric value controlling the size of network nodes (vertices). Default is 15.
#' @param show_labels Logical indicating whether to display node labels. Default is TRUE.
#' @param vertex_label_cex Numeric value controlling the size of vertex labels. Default is 0.7.
#' @param vertex_label_dist Numeric value controlling the distance of labels from vertices. Default is 0.2.
#' @param edge_width Numeric value controlling the width of network edges. Default is 2.
#' @param vertex_label_font Integer specifying the font style for vertex labels (1=plain, 2=bold, 3=italic, 4=bold italic). Default is 2.
#' @param node_alpha Numeric value between 0 and 1 controlling the transparency of nodes. Default is 1 (opaque).
#' @param module_alpha Numeric value between 0 and 1 controlling the transparency of module highlighting. Default is 0.3.
#' @param show_legend Logical indicating whether to display the layer legend. Default is TRUE.
#' @param legend_pos Character string specifying legend position. Options include "bottomright", "bottomleft", "topright", "topleft". Default is "bottomright".
#' @param legend_cex Numeric value controlling the size of legend text. Default is 1.
#' @param legend_pt_size Numeric value controlling the size of legend symbols. Default is 2.5.
#' @param x_intersp Numeric value controlling horizontal spacing in legend. Default is 2.
#' @param y_intersp Numeric value controlling vertical spacing in legend. Default is 2.
#' @param title_font Integer specifying the font style for legend title. Default is 2 (bold).
#' @param title_cex Numeric value controlling the size of legend title. Default is 1.2.
#'
#' @details
#' The PTMN visualization distinguishes between intralayer edges (black lines) showing relationships between traits
#' within the same organ or functional system, and interlayer edges (red lines) representing interactions between
#' traits in different layers. Each layer corresponds to specific functional or structural units
#' such as individual plant organs (leaves, stems, roots) or functional systems.
#'
#' In the default style, cross-layer modules are highlighted with shaded areas, representing tightly connected
#' functional groups that span multiple layers. The circular layout arranges nodes by layer in a circular pattern,
#' making layer organization more visually apparent.
#'
#' Node colors are automatically assigned based on layer membership using the "ggsci::nrc_npg" color palette,
#' ensuring visual distinction between different functional layers.
#'
#' @return
#' No return value. This function is called for its side effect of creating a network plot.
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
#' PTMN_plot(graph, style = "default", vertex_size = 8,
#'           vertex_label_cex = 0.5, edge_width = 2,
#'           show_legend = FALSE)
#'}
#'
#' @seealso
#' \code{\link{PTMN}} for constructing plant trait multilayer networks
#'
#' @importFrom igraph graph_from_data_frame V layout_in_circle
#' @importFrom paletteer paletteer_d
#' @importFrom graphics plot legend
#' @importFrom scales alpha
#' @export
PTMN_plot <- function(data,
                      style = c("default", "circle"),
                      vertex_size = 15,
                      show_labels = TRUE,
                      vertex_label_cex = 0.7,
                      vertex_label_dist = 0.2,
                      edge_width = 2,
                      vertex_label_font = 2,
                      node_alpha = 1,
                      module_alpha = 0.3,
                      show_legend = TRUE,
                      legend_pos = "bottomright",
                      legend_cex = 1,
                      legend_pt_size = 2.5,
                      x_intersp = 2,
                      y_intersp = 2,
                      title_font = 2,
                      title_cex = 1.2
) {
  style <- match.arg(style)
  gg <- data[,c("node.from","node.to")]
  g <- igraph::graph_from_data_frame(gg, directed = FALSE)
  node_info <- unique(rbind(
    data.frame(node = data$node.from, layer = data$layer.from),
    data.frame(node = data$node.to, layer = data$layer.to)))
  nlayer <- length(unique(node_info$layer))
  layer_color_map <- alpha(paletteer_d("ggsci::nrc_npg", n = nlayer), node_alpha)
  node_info$layer_code <- as.numeric(as.factor(node_info$layer))
  igraph::V(g)$layer_code <- node_info$layer_code[match(igraph::V(g)$name, node_info$node)]
  igraph::V(g)$color <- layer_color_map[igraph::V(g)$layer_code]
  cross_groups <- cross_layer_groups(data)
  module_colors <-  alpha("#E2EEF2", 0.3)

  edge_layer_from <- data$layer.from
  edge_layer_to <- data$layer.to
  edge_colors <- ifelse(edge_layer_from != edge_layer_to, "red", "black")
  if (show_labels) {
    igraph::V(g)$label <- igraph::V(g)$name
  } else {
    igraph::V(g)$label <- NA
  }
  unique_layers <- unique(node_info[, c("layer", "layer_code")])
  unique_layers <- unique_layers[order(unique_layers$layer_code), ]
  layer_names <- unique_layers$layer
  if (style == "default") {
    plot(g,
         vertex.color = igraph::V(g)$color,
         vertex.size = vertex_size,
         vertex.label = igraph::V(g)$label,
         vertex.label.cex = vertex_label_cex,
         vertex.label.dist = vertex_label_dist,
         vertex.label.font = vertex_label_font,
         vertex.frame.color = "black",
         vertex.label.color = "black",
         edge.width = edge_width,
         edge.color = edge_colors,
         mark.groups = cross_groups,
         mark.col = module_colors,
         mark.border = "black"
    )
  }
  if (style == "circle") {
    ordered_nodes <- node_info[order(node_info$layer_code), "node"]
    layout <- layout_in_circle(g, order = match(ordered_nodes, V(g)$name))
    plot(g, layout = layout,
         vertex.color = V(g)$color,
         vertex.size = vertex_size,
         vertex.label = V(g)$label,
         vertex.label.cex = vertex_label_cex,
         vertex.label.dist = vertex_label_dist,
         vertex.label.color = "black",
         vertex.label.font = vertex_label_font,
         edge.width = edge_width,
         edge.color = edge_colors
    )
  }
  if (show_legend) {
    layer_colors_for_legend <- paletteer_d("ggsci::nrc_npg", n = nlayer)
    legend(legend_pos,
           legend = layer_names,
           col = "black",
           pt.bg = layer_colors_for_legend,
           pch = 21,
           pt.cex = legend_pt_size,
           title = "Layers",
           cex = legend_cex,
           bty = "n",
           x.intersp = x_intersp,
           y.intersp = y_intersp,
           title.font = title_font,
           title.cex = title_cex)
  }
}

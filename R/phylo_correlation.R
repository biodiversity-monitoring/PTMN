#' Phylogenetically Corrected Correlation Analysis
#'
#' @description
#' This function calculates phylogenetically independent correlations between traits
#' using phylogenetic independent contrasts (PICs). It accounts for phylogenetic
#' relationships when computing correlations, which is important when analyzing
#' trait data from related species.
#'
#' @param traits_matrix A numeric matrix or data frame where rows represent species
#'   and columns represent traits. Row names should contain species names that
#'   match the tip labels in the phylogenetic tree.
#' @param phylo_tree A phylogenetic tree object of class "phylo" (from the ape package).
#'   The tree should contain the same species as in the traits matrix.
#' @param method Character string specifying the correlation method to use.
#'   Options are "pearson" (default) or "spearman".
#'
#' @return A list containing two matrices:
#'   \item{r}{A symmetric correlation matrix with phylogenetically corrected
#'     correlation coefficients}
#'   \item{P}{A symmetric matrix of p-values corresponding to the correlation tests}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the phylogenetic tree is binary and converts it if necessary
#'   \item Matches species between the trait matrix and phylogenetic tree
#'   \item For each pair of traits, calculates phylogenetic independent contrasts (PICs)
#'   \item Computes correlations between PICs instead of raw trait values
#'   \item Handles missing data by using only complete cases for each trait pair
#' }
#'
#' The phylogenetic independent contrasts method removes the effects of phylogenetic
#' relatedness, allowing for proper statistical inference about trait correlations.
#' This is crucial when analyzing data from related species, as standard correlation
#' methods may be biased due to phylogenetic non-independence.
#'
#' @note
#' \itemize{
#'   \item The function requires at least 3 species with complete data for each trait pair
#'   \item Non-binary trees are automatically converted to binary using \code{multi2di()}
#'   \item The function handles missing values by performing pairwise complete case analysis
#' }
#'
#' @examples
#' data(forest_invader_tree)
#' forest_invader_tree
#' data(forest_invader_traits)
#' traits <- forest_invader_traits[,6:73]
#' head(traits)
#' phylo_correlation(traits, forest_invader_tree, method = "pearson")
#'
#' @references
#' Felsenstein, J. (1985). Phylogenies and the comparative method. The American Naturalist, 125(1), 1-15.
#'
#' Harvey, P. H., & Pagel, M. D. (1991). The comparative method in evolutionary biology. Oxford University Press.
#'
#' @export
#' @importFrom ape pic multi2di is.binary is.rooted keep.tip
#' @importFrom stats cor.test complete.cases
phylo_correlation <- function(traits_matrix, phylo_tree, method = "pearson") {
  method <- match.arg(method, choices = c("pearson", "spearman"))
  # Check if the tree is binary, if not, convert to binary using multi2di()
  if (!ape::is.binary(phylo_tree)) {
    cat("Warning: Phylogenetic tree is not binary. Converting to binary tree using multi2di()...\n")
    phylo_tree <- ape::multi2di(phylo_tree)
    cat("Tree converted to binary successfully.\n")
  }
  # Validate tree structure after potential conversion
  if (!ape::is.rooted(phylo_tree)) {
    cat("Warning: Tree is not rooted. Some analyses may require a rooted tree.\n")
  }
  n_traits <- ncol(traits_matrix)
  trait_names <- colnames(traits_matrix)
  # Initialize correlation and p-value matrices
  r_matrix <- matrix(NA, n_traits, n_traits)
  p_matrix <- matrix(NA, n_traits, n_traits)
  rownames(r_matrix) <- colnames(r_matrix) <- trait_names
  rownames(p_matrix) <- colnames(p_matrix) <- trait_names
  # Set diagonal to 1 (self-correlation)
  diag(r_matrix) <- 1
  diag(p_matrix) <- 0
  # Ensure data and tree species order consistency
  species_order <- phylo_tree$tip.label
  # Check if all species in the tree are present in the trait matrix
  missing_species1 <- setdiff(species_order, rownames(traits_matrix))
  if (length(missing_species1) > 0) {
    cat("Warning: The following species are in the tree but not in the trait matrix:",
        paste(missing_species1, collapse = ", "), "\n")
  }
  missing_species2 <- setdiff(rownames(traits_matrix), species_order)
  if (length(missing_species2) > 0) {
    cat("Warning: The following species are in the trait matrix but not in the tree:",
        paste(missing_species2, collapse = ", "), "\n")
  }
  # Only keep species that are present in both tree and trait matrix
  common_species <- intersect(species_order, rownames(traits_matrix))
  if (length(common_species) < 3) {
    stop("Insufficient overlap between tree species and trait matrix species (need at least 3).")
  }
  traits_matrix <- traits_matrix[common_species, ]
  species_order <- common_species
  # Calculate phylogenetically independent correlations for each pair of traits
  for (i in 1:(n_traits-1)) {
    for (j in (i+1):n_traits) {
      trait1 <- traits_matrix[, i]
      trait2 <- traits_matrix[, j]
      # Remove missing values
      complete_cases <- complete.cases(trait1, trait2)
      if (sum(complete_cases) < 3) {
        r_matrix[i, j] <- r_matrix[j, i] <- NA
        p_matrix[i, j] <- p_matrix[j, i] <- NA
        next
      }
      # Get species with complete data
      species_with_data <- species_order[complete_cases]
      trait1_clean <- trait1[complete_cases]
      trait2_clean <- trait2[complete_cases]
      # Add species names to trait vectors
      names(trait1_clean) <- species_with_data
      names(trait2_clean) <- species_with_data
      # Prune phylogenetic tree
      tryCatch({
        tree_clean <- ape::keep.tip(phylo_tree, species_with_data)
        # Check if pruned tree is valid
        if (is.null(tree_clean) || length(tree_clean$tip.label) < 3) {
          r_matrix[i, j] <- r_matrix[j, i] <- NA
          p_matrix[i, j] <- p_matrix[j, i] <- NA
          next
        }
        # Ensure the pruned tree is still binary (in case pruning created polytomies)
        if (!ape::is.binary(tree_clean)) {
          tree_clean <- ape::multi2di(tree_clean)
        }
        # Ensure data order matches the pruned tree
        trait1_final <- trait1_clean[tree_clean$tip.label]
        trait2_final <- trait2_clean[tree_clean$tip.label]
        # Calculate phylogenetic independent contrasts
        pic1 <- ape::pic(trait1_final, tree_clean)
        pic2 <- ape::pic(trait2_final, tree_clean)
        # Calculate correlations between PICs
        cor_test <- cor.test(pic1, pic2, method = method)
        r_matrix[i, j] <- r_matrix[j, i] <- cor_test$estimate
        p_matrix[i, j] <- p_matrix[j, i] <- cor_test$p.value
      }, error = function(e) {
        cat("Warning: Analysis failed for traits", trait_names[i], "and", trait_names[j], ":", e$message, "\n")
        r_matrix[i, j] <<- r_matrix[j, i] <<- NA
        p_matrix[i, j] <<- p_matrix[j, i] <<- NA
      })
    }
  }
  return(list(r = r_matrix, P = p_matrix))
}

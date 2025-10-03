#' Phylogenetic tree for forest invader trait analysis
#'
#' A phylogenetic tree object containing 76 woody plant species (44 native and 32 non-native)
#' used in the analysis of functional trait syndromes for forest invaders in North American
#' deciduous forests. The tree is based on species-level phylogeny with aged branch lengths
#' and includes native, naturalized, and invasive species classifications.
#'
#' @format A phylo object with:
#' \describe{
#'   \item{edge}{A 142 x 2 matrix defining the tree topology}
#'   \item{edge.length}{A numeric vector of 142 branch lengths}
#'   \item{Nnode}{Number of internal nodes (67)}
#'   \item{node.label}{Character vector with node labels including taxonomic groups}
#'   \item{tip.label}{Character vector of 76 species names}
#' }
#'
#' @details
#' The phylogenetic tree was constructed using correlation matrix of phylogenetic distances
#' calculated with an aged tree based on the species-level phylogeny. The tree includes
#' representatives from major angiosperm families and is used for phylogenetically
#' independent contrasts in trait analysis. Branch lengths reflect evolutionary time,
#' and the tree structure allows for proper phylogenetic corrections in comparative
#' analyses of functional traits between native and invasive species.
#'
#' Species include shrubs and vines from deciduous forests, with invasive species
#' representing those managed as high-impact invaders in northeastern United States,
#' and naturalized species representing non-native species that have established
#' but are not considered highly invasive.
#'
#' @source
#' Tree construction based on Jo et al. (2016) phylogeny with additional species
#' placement using updated familial/generic phylogenies.
#'
#' @references
#' Fridley, J. D., Bauerle, T. L., Craddock, A., Ebert, A. R., Frank, D. A., Heberling, J. M., Hinman, E. D., Jo, I., Martinez, K. A., Smith, M. S., Woolhiser, L. J., & Yin, J. (2022). Fast but steady: An integrated leaf‐stem‐root trait syndrome for woody forest invaders. Ecology Letters, 25(4), 900-912.
#' \url{ https://doi.org/10.1111/ele.13967 }
#'
#' @examples
#' data(forest_invader_tree)
#' print(forest_invader_tree)
"forest_invader_tree"

# Load Data
traits <- read.csv("data-raw/shrub_traits_final.csv",row.names=1)
trcont <- traits[,sapply(traits,class)!="character"]  # continuous traits only

############################################################
## Data transformation
# Inspection of distributions suggests the majority of traits show substantial right skew,
# which are log-transformed prior to analysis to reduce leverage of high trait values in linear models (eg PCA):
for(i in 1:dim(trcont)[2]) {  # figures not shown here
  par(mfrow=c(1,2),mar=c(5,5,1,1))
  plot(density(trcont[,i],na.rm=T),col="black",main=names(trcont[i]))
  plot(density(log(trcont[,i]),na.rm=T),col="red2",main=paste("log ",names(trcont[i]),sep=""))
}
# Columns of high right skew to log-transform
log.cols <- c(8,9,11,12,14,17,23,24,26,28,29,30,33,36,37,40,41:45,47:49,51,53:57,59,63,64,66,67,70,72,76)
traits2 <- traits
for(i in log.cols) { traits2[,i] = log(traits2[,i]) }  # traits2 is trait matrix included logged columns

############################################################
## Principal Components Analysis
# PPCA is performed in R in the pcaMethods library, installed via the Bioconductor platform (<https://www.bioconductor.org/install/>).
# install.packages("BiocManager")
# BiocManager::install("pcaMethods")
library(pcaMethods)    # https://academic.oup.com/bioinformatics/article/23/9/1164/272597

# Input matrix of the continuous 68 traits for 76 species, including those log-transformed, and standardize columns to zero mean and unit variance.
pcamat <- scale(traits2[,sapply(traits,class)!="character"])

# To determine how many principal components are to be used to impute missing values,
# pcaMethods includes the 'kEstimate' function that calculates error via cross validation,
# with choices of two different types of error (NRMSEP, normalized root mean square error,
# and Q2 distance, more suitable for smaller datasets; Stacklies et al. 2007).
# In this case both error methods give the same answer of two PCs. (This process takes about 15 min.)
# not run:
# estimate <- kEstimate(pcamat,method="ppca",em="q2",evalPcs=1:10)
# estimate$bestNPcs  #optimal number of PCs = 2

# run PPCA
ppca1 <- pca(pcamat,method="ppca",seed=123,nPcs=2)
ppca1    # PPCA eigenvalues: PC1 = 18%, PC2 = 14% (32% of total variance)
impmat <- completeObs(ppca1) # complete dataset (with imputed values for NAs)

library(ape)
library(geiger)

tree <- read.tree(file="data-raw/finaltree.tre") # phylogeny
traits3 <- traits2 # includes log transformations
targets <- sub(" ","_",traits3$taxon)
# simplify infraspecific taxon names
targets[27] <- "Euonymus_hamiltonianus"
targets[58] <- "Sambucus_nigra"
targets[39] <- "Lonicera_japonica"
targets[49] <- "Lonicera_tatarica"
targets[51] <- "Lonicera_xbella"
targets[71] <- "Viburnum_opulus"
targets[70] <- "Viburnum_nudum"
targets[62] <- "Symphoricarpos_albus"
rownames(traits3) = targets
name.check(tree,traits3)

trait <- traits[,c("taxon", "family", "growthform", "nativity", "status")]
trait$taxon[27] <- "Euonymus hamiltonianus"
trait$taxon[58] <- "Sambucus nigra"
trait$taxon[39] <- "Lonicera japonica"
trait$taxon[49] <- "Lonicera tatarica"
trait$taxon[51] <- "Lonicera xbella"
trait$taxon[71] <- "Viburnum opulus"
trait$taxon[70] <- "Viburnum nudum"
trait$taxon[62] <- "Symphoricarpos albus"
trait

forest_invader_traits <- cbind(trait, impmat)
colnames(forest_invader_traits) <- c("spesies", "family", "growthform", "nativity", "status",
                                     "Mass_leaf", "SLA", "LDMC", "Lobe", "PA_leaf", "Thick_leaf",
                                     "Stomate_size", "Stomate_index", "Rd_area", "Rd_mass", "Asat_area", "Asat_mass",
                                     "Vc_max", "J_max", "AQY", "LSP", "PNUE", "WUE", "CC_area", "CC_mass", "N_mass",
                                     "N_area", "N_litter", "N_res", "C_mass", "C_area", "C_litter", "Solubles", "Hemi",
                                     "Lignin", "Cellulose", "Ash", "Phenolics", "Chl_area", "Chl_mass", "Chl_ab",
                                     "Decomp_leaf", "DNA", "SRL", "RTD", "RDMC", "N_root", "C_root", "Lignin_root",
                                     "Decomp_root", "LAgain", "LeafRate_max", "NPP_root", "SS_root", "SS_stem", "Starch_root",
                                     "Starch_stem", "NSC_root", "NSC_stem", "Vessel_freq", "Metaxylem_diam",
                                     "Earlywood_diam", "Latewood_diam", "Cond_stem", "LeafFall50", "Lifespan_leaf",
                                     "EmergeDuration", "LeafDuration", "LeafTurnover", "FallDuration", "Chl_sun50",
                                     "Chl_shade50", "LeafOut")
row.names(forest_invader_traits) <- forest_invader_traits$spesies
forest_invader_traits$nativity <- ifelse(forest_invader_traits$nativity == "native", 0, 1)
forest_invader_traits

tree$tip.label <- gsub("_", " ", tree$tip.label)
forest_invader_tree <- tree
forest_invader_tree

# Save data to the package
usethis::use_data(forest_invader_traits, internal = FALSE, compress = "xz", overwrite = TRUE)
usethis::use_data(forest_invader_tree, internal = FALSE, compress = "xz", overwrite = TRUE)


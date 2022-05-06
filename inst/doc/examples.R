## -----------------------------------------------------------------------------
# Load the package into R
library(rdiversity)

# Initialise data
partition <- data.frame(a=c(1,1),b=c(2,0),c=c(3,1))
row.names(partition) <- c("cows", "sheep")

## -----------------------------------------------------------------------------
# Generate metacommunity object
meta <- metacommunity(partition = partition)

## -----------------------------------------------------------------------------
# Initialise similarity matrix
s <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
row.names(s) <- c("cows", "sheep")
colnames(s) <- c("cows", "sheep")

# Generate similarity object 
s <- similarity(similarity = s, dat_id = "my_taxonomic")

# Generate metacommunity object
meta <- metacommunity(partition = partition, similarity = s)

## -----------------------------------------------------------------------------
# Initialise distance matrix
d <- matrix(c(0, 0.7, 0.7, 0), nrow = 2)
row.names(d) <- c("cows", "sheep")
colnames(d) <- c("cows", "sheep")

# Generate distance object
d <- distance(distance = d, dat_id = "my_taxonomic")

# Convert the distance object to similarity object (by means of a linear or exponential transform)
s <- dist2sim(dist = d, transform = "linear")

# Generate metacommunity object
meta <- metacommunity(partition = partition, similarity = s)

## -----------------------------------------------------------------------------
# Initialise data
partition <- data.frame(a=c(1,1),b=c(2,0),c=c(3,1))
row.names(partition) <- c("cows", "sheep")

# Generate a metacommunity object
meta <- metacommunity(partition)

# Calculate diversity
norm_sub_alpha(meta, 0:2)

## -----------------------------------------------------------------------------
# Initialise data
partition <- data.frame(a=c(1,1),b=c(2,0),c=c(3,1))
row.names(partition) <- c("cows", "sheep")

# Generate a metacommunity object
meta <- metacommunity(partition)

# Calculate the species-level component for normalised alpha
component <- norm_alpha(meta)

# Calculate normalised alpha at the subcommunity-level 
subdiv(component, 0:2)

# Likewise, calculate normalised alpha at the metacommunity-level 
metadiv(component, 0:2)

## -----------------------------------------------------------------------------
# Calculate all subcommunity diversity measures
subdiv(meta, 0:2)

# Calculate all metacommunity diversity measures
metadiv(meta, 0:2)

## -----------------------------------------------------------------------------
# Taxonomic lookup table
Species <- c("tenuifolium", "asterolepis", "simplex var.grandiflora", "simplex var.ochnacea")
Genus <- c("Protium", "Quararibea", "Swartzia", "Swartzia")
Family <- c("Burseraceae", "Bombacaceae", "Fabaceae", "Fabaceae")
Subclass <- c("Sapindales", "Malvales", "Fabales", "Fabales")
lookup <- cbind.data.frame(Species, Genus, Family, Subclass)

# Partition matrix
partition <- matrix(rep(1, 8), nrow = 4)
colnames(partition) <- LETTERS[1:2]
rownames(partition) <- lookup$Species


## -----------------------------------------------------------------------------
values <- c(Species = 0, Genus = 1, Family = 2, Subclass = 3, Other = 4)

## -----------------------------------------------------------------------------
d <- tax2dist(lookup, values)

## -----------------------------------------------------------------------------
s <- dist2sim(d, "linear")

## -----------------------------------------------------------------------------
meta <- metacommunity(partition, s)

## -----------------------------------------------------------------------------
meta_gamma(meta, 0:2)

## -----------------------------------------------------------------------------
# Example data
tree <- ape::rtree(4)
partition <- matrix(1:12, ncol=3)
partition <- partition/sum(partition)

## -----------------------------------------------------------------------------
d <- phy2dist(tree)

## -----------------------------------------------------------------------------
s <- dist2sim(d, "linear")

## -----------------------------------------------------------------------------
meta <- metacommunity(partition, s)

## -----------------------------------------------------------------------------
meta_gamma(meta, 0:2)

## -----------------------------------------------------------------------------
tree <- ape::rtree(4)
partition <- matrix(1:12, ncol=3)
partition <- partition/sum(partition)
colnames(partition) <- letters[1:3]
row.names(partition) <- paste0("sp",1:4)
tree$tip.label <- row.names(partition)

## -----------------------------------------------------------------------------
s <- phy2branch(tree, partition)

## -----------------------------------------------------------------------------
meta <- metacommunity(partition, s)

## -----------------------------------------------------------------------------
meta_gamma(meta, 0:2)

## ---- eval = FALSE------------------------------------------------------------
#  library(rdiversity)
#  vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
#  #read in twice: first for the column names then for the data
#  tmp_vcf <- readLines(vcf_file)
#  vcf_data <- read.table(vcf_file, stringsAsFactors = FALSE)
#  # filter for the columns names
#  vcf_names <- unlist(strsplit(tmp_vcf[grep("#CHROM",tmp_vcf)],"\t"))
#  names(vcf_data) <- vcf_names
#  partition <- cbind.data.frame(A = c(rep(1, 9), rep(0, 9)), B = c(rep(0, 9), rep(1, 9)))
#  partition <- partition/sum(partition)

## ---- eval = FALSE------------------------------------------------------------
#  d <- gen2dist(vcf)

## ---- eval = FALSE------------------------------------------------------------
#  s <- dist2sim(d, transform = 'l')

## ---- eval = FALSE------------------------------------------------------------
#  rownames(partition) <- rownames(s@similarity)
#  meta <- metacommunity(partition, s)

## ---- eval = FALSE------------------------------------------------------------
#  norm_meta_beta(meta, 0:2)

## -----------------------------------------------------------------------------
partition <- matrix(sample(6), nrow = 3)
rownames(partition) <- paste0("sp", 1:3)
partition <- partition / sum(partition)

d <- matrix(c(0,.75,1,.75,0,.3,1,.3,0), nrow = 3)
rownames(d) <- paste0("sp", 1:3)
colnames(d) <- paste0("sp", 1:3)
d <- distance(d, "my_taxonomy")
s <- dist2sim(d, "linear")

meta <- metacommunity(partition, s)

## -----------------------------------------------------------------------------
partition <- matrix(sample(6), nrow = 3)
rownames(partition) <- paste0("sp", 1:3)
partition <- partition / sum(partition)

s <- matrix(c(1,.8,0,.8,1,.1,0,.1,1), nrow = 3)
rownames(s) <- paste0("sp", 1:3)
colnames(s) <- paste0("sp", 1:3)
s <- similarity(s, "my_functional")

meta <- metacommunity(partition, s)

## -----------------------------------------------------------------------------
tree <- ape::rtree(5)
tree$tip.label <- paste0("sp", 1:5)

partition <- matrix(rep(1,10), nrow = 5)
row.names(partition) <- paste0("sp", 1:5)
partition <- partition / sum(partition)
s <- phy2branch(tree, partition)
meta <- metacommunity(partition, s)

new_partition <- matrix(sample(10), nrow = 5)
row.names(new_partition) <- paste0("sp", 1:5)
new_partition <- new_partition / sum(new_partition)

new_meta <- repartition(meta, new_partition)


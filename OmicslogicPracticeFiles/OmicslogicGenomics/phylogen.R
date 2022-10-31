getwd()
install.packages("ape")
# Load libraries for analysis
library(stats)
library(ape)
install.packages("ade4")
library(ade4)

# loading data
dna <- read.dna(file = "sample.fasta", format = "fasta")

#generating phylogenetic distances
Dist <- dist.dna(dna, model = "TN93")

# generating heatmap
temp <- as.data.frame(as.matrix(Dist)) 
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)

#coloring the heatmap
# convert into matrix
temp1 <-as.matrix(Dist)
# Define position of plot
par(mar = c(0.05, 4, 3.2, 0.05))


# Plot Distance matrix in the form of colored heatmap
image(x = 1:40, y = 1:40,temp1, col = rev(heat.colors(100)), xaxt = "n", yaxt = "n", xlab = "Samples", ylab = "Samples")
axis(side = 2, at = 1:40, lab = rownames(dna), las = 2, cex.axis = 0.6)
axis(side = 3, at = 1:40, lab = rownames(dna), las = 3, cex.axis = 0.6)

# Build phylogenetic tree using NJ method
tree1 <- nj(Dist)
# plot tree
plot(tree1, cex = 0.6)

# Build phylogenetic tree (dendrogram) using hclust method
tree2 <- hclust(Dist)
# plot
plot(tree2, labels = NULL, hang = 0.1, check = TRUE,cex=0.6,  axes = TRUE, frame.plot = FALSE, ann = TRUE, main = "", sub = NULL, xlab = NULL, ylab = "Height")
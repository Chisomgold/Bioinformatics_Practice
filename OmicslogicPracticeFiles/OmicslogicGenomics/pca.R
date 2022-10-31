#running principal component analysis

library(ggplot2)
library(ggfortify)
library(scatterplot3d)

#input file
input_file = 'https://raw.githubusercontent.com/PineBiotech/omicslogic/master/CellLines_15Genes_1.txt'

#define the data
full_table <- read.table (input_file, header = TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
expressions <- as.matrix(full_table[2:nrow(full_table), 2:ncol(full_table)]) 

#transpose the matrix
expressionst <- t(expressions)

#prepare the 'pca' object
pca <- (expressionst)
pca <- prcomp(expressionst)

#print the pca object:

#plot the pca using native scatterplot function plot()
plot(pca$x, pca$y)

#Create plot
plot2 <- autoplot(pca, label = TRUE, label.size = 3 , colour = "red" )

#Show plot
plot(plot2)

#print SD of PCA object
print(pca$sd)

#assign sample names
Samples <- colnames(expressions)

#run PCA and test length (number of rows)
pca <- prcomp((t(expressions)), scale. = TRUE, center= TRUE)
pca_res <- data.frame(Samples, pca$x)

#print the length of pca_res object to make sure you have 15 rows as expected
print(length(pca_res))

#Save PCA table
write.table(pca_res, file="PCA_table1.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=TRUE)

#MODIFYING PLOT APPEARANCE
#Load libraries
#Load data
ExpressionTable <- read.table('https://raw.githubusercontent.com/PineBiotech/bioinformatics/master/15gene_transposed.txt', sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=TRUE, row.names=1)

#Take data without groups
df <- ExpressionTable[2:16]

#Generate PCA components
pca_res <- prcomp(df, scale.=TRUE)

#Plot PCA plot
plot <- (pca_res)
plot(plot)

# add color based on group
#Create Frame
plot1 <- autoplot(pca_res, data=ExpressionTable, colour='Group')
plot(plot1)


#Add labels
plot2 <- autoplot(pca_res, data=ExpressionTable, colour='Group', label=TRUE)
plot(plot2)

#Create Frame
plot3 <- autoplot(pca_res, data=ExpressionTable, colour='Group', frame=TRUE)
plot(plot3)


#Print top scaled values of pca
head(pca_res$scale)

#3D viz
PCA_table <- as.matrix(pca$x) 
scatterplot3d(PCA_table[, 1 : 3 ], angle= 55, label=TRUE)

#adding frametype
autoplot(pca, data=expressionst, frame= TRUE , label = TRUE , label.size = 3 , frame.type= 'norm' )

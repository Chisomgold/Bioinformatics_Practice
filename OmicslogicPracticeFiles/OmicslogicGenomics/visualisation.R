setwd("C:/Users/hp/OneDrive/Desktop/OmicslogicGenomics")

#loading data
df <- read.table('https://raw.githubusercontent.com/pine-bio-support/DataScience/main/Final_cell_lines_RNA-expression_FPKM_values_1000genes_with_NA.txt', sep ='\t', header = TRUE, row.names=1)

#write the first 10 rows of the dataset:
head(df)
# Check dimension of data (how many rows and columns)
dim(df)
#check datatype
str(df)

# Extract the sample names (column names)
Samples <- colnames(df)
print(Samples)

#Extract row ids from data
Genes  <- rownames(df)
print(Genes)

# remove NA values to clean the data
df1 <-  na.exclude(df)
# To check the dimension of cleaned data
dim(df1)

# Descriptive statistics, i.e., compute the summary statistics of data & Visualisation
summ <- summary(df1)
# print summary statistics
print(summ)
# Write summary statistics of data  into a file and export it
write.table(summ, file="stat_sum-FPKM-data.txt", col.names=TRUE, sep="\t")

# Draw boxplot for all samples (for FPKM values)
boxplot(df1)
# Add the main title, axis titles, and color 
boxplot(df1, main="Boxplot for FPKM data", xlab="Samples", ylab="Gene expression (FPKM)", col="red", las=2, cex.axis = 0.65)
#resulting plot doesn't look good, so

#perform log scale transformation of data
log_df <- log(df1+1)
#Check the dimension of data
dim(log_df)
# check the summary statistics of data after log transformation
summ_log <- summary(log_df)
# Print summary statistics of data 
summ_log
# Write summary statistics to a file and export it
write.table(summ_log, file="stat_sum-log-data.txt", col.names=TRUE, sep="\t")
# Draw boxplot for all samples (for FPKM values)
boxplot(log_df, main="Boxplot for log-transformed Data", xlab="Samples", ylab="Gene expression (log[FPKM])", col="red", las=2, cex.axis = 0.7)

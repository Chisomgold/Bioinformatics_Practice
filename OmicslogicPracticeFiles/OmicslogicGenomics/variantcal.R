#Analysing vcf files
# Load the libraries
library('tidyverse')
# Load sample vcf into a dataframe object
sample_vcf_tp53 <- read.table('https://raw.githubusercontent.com/pine-bio-support/Merge-VCF-files/main/aneuploid_samples_freebayes_tp53.vcf', header = FALSE, comment.char = "#", sep = "\t")
# Display first few lines of the vcf dataframe
head(sample_vcf_tp53)
#Renaming vcf file column headers
# Define column names
names(sample_vcf_tp53) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
head(sample_vcf_tp53)

# Select only the first 7 columns and ignore the rest
sample_vcf_tp53 <- select(sample_vcf_tp53, c('CHROM','POS','ID','REF','ALT','QUAL','FILTER'))
head(sample_vcf_tp53)

###References
#Load Reference VCF to a Dataframe
clinvar_vcf_tp53 <- read.table('https://raw.githubusercontent.com/pine-bio-support/Merge-VCF-files/main/clinVar_all_tp53_edt.vcf',
                               header = FALSE, comment.char = "#", sep = "\t")
head(clinvar_vcf_tp53)

# Define column names
names(clinvar_vcf_tp53) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
head(clinvar_vcf_tp53)

#Extract clinical significance
clinvar_vcf_tp53 <- clinvar_vcf_tp53 %>% separate(INFO, c('INFO', 'Significance'), sep=';CLNSIG=')
head(clinvar_vcf_tp53)

#Find common variants between sample and reference dataframe
sample_tp53_clnvar <- inner_join(sample_vcf_tp53, clinvar_vcf_tp53, by = c("CHROM" = "CHROM", "POS"="POS","REF" = "REF"))
head(sample_tp53_clnvar)

#Select only relevant columns
sample_tp53_clnvar <- select(sample_tp53_clnvar, c("CHROM", "POS", "REF", "ALT.x", "QUAL.x", "ALT.y", "Significance"))
#Rename Columns
sample_tp53_clnvar <- rename(sample_tp53_clnvar, "ALT.sample"=ALT.x, "ALT.clnvar"=ALT.y, "QUAL"=QUAL.x)
head(sample_tp53_clnvar)

#Save the output
write.table(sample_tp53_clnvar,file="sample_tp53_clnvar_annotated.txt", sep='\t',  quote = F, row.names = FALSE)

#identifying pathogenic variants
# Tabulate the frequency of diverse clinical significant variants
table(sample_tp53_clnvar$Significance)

# Extract clinically pathogenic variants.
sample_tp53_clnvar_pathogenic <- filter(sample_tp53_clnvar, Significance == "Pathogenic")
sample_tp53_clnvar_pathogenic

#viewing table
view(sample_tp53_clnvar_pathogenic)
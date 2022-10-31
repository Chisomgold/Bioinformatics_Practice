getwd()
#setting working directory
setwd("C:/Users/hp/OneDrive/Desktop/OmicslogicGenomics/")

#DNA seqs analysis
DNAstring1 <- "ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT"
RNA <- chartr("T", "U", DNAstring1)
print(RNA)

#loading biostring library for translation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

#using biostrings for translation
library(Biostrings) 
dna <- DNAString(DNAstring1) 
translate(dna)

#trying alignment
BiocManager::install("DECIPHER")
library(DECIPHER)

dna2 <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC")
dna3 <- DNAString("CAAAGCGGTCCATGAGATGC")

myalign1 <- pairwiseAlignment(dna, dna2)
myalign2 <- pairwiseAlignment(dna, dna3)

#viewing alignment results
myalign1
myalign2 #this is minus because the default alignment is global which is not 
# effective for highly diverse seqs

#performing local alignment
localalign <- pairwiseAlignment(dna, dna3, type="local")
print(localalign) #results or score is better

#performing overlap alignment
overlapalign <- pairwiseAlignment(dna, dna3, type="overlap")
print(overlapalign)

#trying multiple sequence alignment
#first make a list of seqs
DNASeq <- list()
DNASeq[1] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGAT")
DNASeq[2] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAGATGCAATGCCCTTGTATATATATGTGATTTAC")
DNASeq[3] <- DNAString("ATTTTAAGTAGTTAAGCCAGTGCCCGATGCAAAGCGGTCCATGAATGCAATGCCCTTGTATATATATGTGATAA")


seqs <- DNAStringSet(unlist(DNASeq))

#perform the alignment of multiple sequences
aligned <- AlignSeqs(seqs, verbose = FALSE)

#print out the aligned sequences 
print(aligned)

#to save results as fasta
myFASTA <- "my.fasta"
writeXStringSet(aligned, myFASTA)
#view or save the alignment in an HTML file for a browser view 
TF <- tempfile("plot___", fileext = ".html")
TF <- BrowseSeqs(aligned, highlight=0, htmlFile = TF)
file.copy(TF, './')

#in case of loading seqs from a fasta file
fas <- "myfasta.fasta"
seqs <- readDNAStringSet(fas)
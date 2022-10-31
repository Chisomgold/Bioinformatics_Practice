# making a function that will generate report on variants or mutations
#take in seq positions of interest and run automatically

#defining positions of interest
myPos <- list(1878,14408)
myPos[1] <- 8782
myPos[2] <- 14408
myPos[3] <- 28144

generateReport <- function(myPos) {
  j = 1
  for (j in 1:PosLength) {
    
    #select 50 NT section around the variant of interest
    startnt <- as.numeric(myPos[j])-25
    endnt <- as.numeric(myPos[j])+25
    
    #select only data where any sample has a variant,
    msa_sel <- msa[startnt:endnt,] #select specified rows from the letter data frame
    msa_sel[] <- lapply(msa_sel, as.character) #if number of characters in each sample does not match, we can convert the values to characters
    msa_sel$Variant <- ifelse(msa_sel[1] == msa_sel[2] & msa_sel[1] == msa_sel[3] & msa_sel[1] == msa_sel[4] & msa_sel[1] == msa_sel[5] & msa_sel[1] == msa_sel[6], "", "V")
    msa_selT <- t(msa_sel) #transpose data frame
    msa_sel1 <- msa1[startnt:endnt,] #select specified rows from the letter data frame
    msa_sel1T <- t(msa_sel1) #transpose data frame
    
    #transform the data for ggplot
    library(reshape2)
    melted_mat <- melt(msa_selT)
    melted_mat1 <- melt(msa_sel1T)
    
    #create a tile for the MSA plot
    plottitle <- as.character(c("Complete MSA plot for the selected sequence region","(",startnt," - ",endnt,")"))
    plottitle1 <- paste( unlist(plottitle), collapse='')
    
    #for graphics display
    dev.off()
    #draw the plot
    (plot1 <- ggplot(melted_mat) + #add data and fill cells
        geom_tile(data = melted_mat, aes(x=Var2, y=Var1, fill=value), colour = "black") + #add black border around cells
        geom_text(data = melted_mat, aes(x=Var2, y=Var1,label = value), size=4, family = "Avenir Next") + 
        scale_fill_manual(values = c("-" = "lightgray", "A" = "lightgreen", "C" = "pink", "G" = "lightblue", "T" = "yellow", "V"="red", ""="white")) +
        coord_equal() +
        xlab("Nucleotide Position") +
        ylab("Samples") +
        labs(fill = "NT") +
        labs(title=plottitle1, size=2) +
        scale_x_continuous(breaks = seq(startnt, endnt, by = 2), expand = c(0, 0)) +
        theme(axis.text = element_text(size=10, family = "Avenir Next"),
              plot.title = element_text(size = 10, face = "bold", family = "Avenir Next"),
              legend.title = element_text(color = "gray", size = 8),
              legend.text = element_text(color = "gray"),
              axis.title = element_text(size=6, vjust = 2, face="italic"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4),
              axis.text.y = element_text(size = 6),
              panel.border = element_rect(colour = "black", fill=NA, size=1)
        ))
    
    #make the first plot
    plot(plot1)
    
    ###second plot
    
    #Report on the number of mutations in each sample:
    x <- ncol(msa_sel)-1
    
    #create an empty data frame to store the new data
    daf2 <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
    colnames(daf2) <- colnames(msa_sel)
    
    #create a title for the second plot
    varType <- as.factor(msa_sel[msa_sel$Variant == "V",])
    varType1 <- as.character(varType)
    plot2title <- as.character(c("Variant at ",myPos[j], "( ", varType1, " )"))
    plot2title1 <- paste( unlist(plot2title), collapse='')
    
    i=2
    
    for (i in 2:x) {
      daf2[,i] <- ifelse(msa_sel[,1]== msa_sel[,i],0,1)
    }
    
    
    #drop the variant column 
    
    daf3 <- subset(daf2, select = -Variant )
    
    
    #count the number of mutations in each sample
    summ1 <- colSums(daf3)
    
    #make a barplot
    par(mar=c(11,4,4,4))
    
    
  barplot1 <- barplot(summ1, las=2, cex.axis = 0.5, cex.lab=0.5, cex.names=0.5, col="magenta", ylab="Number of mutations",space=0.2, main=plot2title1)
  }
}

generateReport(myPos)
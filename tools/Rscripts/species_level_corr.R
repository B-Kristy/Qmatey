# Load libraries
library(plyr); library(dplyr)
library(car)
library(MASS)
library(data.table)
library(ggplot2)
library(ggcorrplot)
library(reshape)
library(stringr)

args <- commandArgs(TRUE)

species <- read.table("proj_taxainfo_mean.txt", header = T, sep="\t", check.names=FALSE, fill=TRUE)

for (i in c(5,10)) {
  metag <- subset(species, select=-c(tax_id,genus,family,order,class,phylum,kingdom,superkingdom,taxname))
  metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
  metag$count <- rowSums(metag[,2:ncol(metag)] > "0")
  metag$percent <- ((metag$count)/(ncol(metag)-2))*100
  metag <- subset(metag, percent >= i )
  metag <- subset(metag, select=-c(count,percent))
  metag <- setNames(data.frame(t(metag[,-1])), metag[,1])
  metag_corr <- cor(metag, method = c("spearman")); axis_density <- 2000/nrow(metag_corr)
  metag_corr <- melt(metag_corr)
  metag_pmat <- cor_pmat(metag, method = c("spearman"), exact=FALSE)
  metag_pmat <- melt(metag_pmat)
  corr <- merge(metag_corr, metag_pmat, by=c("X1","X2"))
  colnames(corr) <- c("Var1","Var2","coeff","pvalue")
  corr <- subset(corr, pvalue <= 0.05)
  
  
  plot <- ggplot(corr, aes(x=Var1, y=Var2, fill=coeff, size= ppvalue)) +
    geom_point(aes(size=-pvalue), shape=21) + scale_fill_gradient2(low="red", mid="white", high="cornflowerblue") +
    theme_bw() + coord_equal() + scale_size(guide = 'none') +
    labs(x="",y="",fill="Correlation\nCoeficient",size="p-value") +
    theme(axis.text.x=element_text(size=axis_density, angle=45, vjust=1, hjust=1, margin=margin(0,0,0,0)),
          axis.text.y=element_text(size=axis_density, margin=margin(0,0,0,0)), panel.grid.major=element_line(colour = "grey95"),
          legend.title=element_text(size=15), legend.text=element_text(size=20),legend.key.size = unit(0.5, "in"),
          plot.title = element_text(size=15)) +
    labs(title= "Species-Level Correlogram")
  ggsave(filename=paste("Metagenome",i,"_species_corr_spearman.tiff",sep=""), plot=plot, width=35, height=35, dpi=600, compression = "lzw", limitsize = FALSE)
}

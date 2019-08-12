# Load libraries
library(plyr); library(dplyr)
library(car)
library(MASS)
library(data.table)
library(ggplot2)
library(ggcorrplot)
library(reshape)
library(stringr)

setwd ("/home/rkuster/Desktop/DIV/metagenome/Results")
metag <- read.table("DIV_taxainfo_mean.txt", header=T, sep="\t", fill =TRUE)
metag <- subset(metag, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag$count <- rowSums(metag[,2:ncol(metag)] > "0")
metag$percent <- ((metag$count)/(ncol(metag)-2))*100
metag <- subset(metag, percent >= "20")
metag <- subset(metag, select=-c(count,percent))
metag = setNames(data.frame(t(metag[,-1])), metag[,1])


metag_corr <- cor(metag, method = c("spearman"))
metag_pmat <- cor_pmat(metag, method = c("spearman"))


plot <- ggcorrplot(metag_corr, hc.order = TRUE,
           type = "lower", p.mat = metag_pmat, insig = "blank",
           lab = FALSE, lab_size = 3, 
           method="circle",  sig.level = 0.05, show.diag = FALSE,
           colors = c("red", "white", "cornflowerblue"),
           tl.cex=7,
           title="Correlogram of Leaf microbiome (Beauregard x Tanzania bi-parental population)", 
           ggtheme=theme_bw)

ggsave(filename="Metagenome20_corr_spearman.tiff", plot=plot, width=35, height= 35, dpi=300, compression = "lzw")



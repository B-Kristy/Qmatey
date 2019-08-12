# Load libraries
library(plyr); library(dplyr)
library(car)
library(MASS)
library(data.table)
library(ggplot2)
library(ggcorrplot)
library(reshape)
library(stringr)





setwd("/media/drive_d/NKB/metagenome/results")
metag <- read.table("proj_taxainfo_mean.txt", header=T, sep="\t", fill =TRUE)
metag <- subset(metag, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag$count <- rowSums(metag[,2:ncol(metag)] > "0")
metag$percent <- ((metag$count)/(ncol(metag)-2))*100
metag <- subset(metag, percent >= "20")
metag <- subset(metag, select=-c(count,percent))
metag = setNames(data.frame(t(metag[,-1])), metag[,1])


metag_corr <- cor(metag, method = c("spearman"))
metag_pmat <- cor_pmat(metag, method = c("spearman"), exact=FALSE)


plot <- ggcorrplot(metag_corr, hc.order = TRUE,
           type = "lower", insig = "blank",
           lab = FALSE, lab_size = 3, 
           method="circle",  sig.level = 0.05, show.diag = FALSE,
           colors = c("red", "white", "cornflowerblue"),
           tl.cex=7,
           title="Correlogram of Leaf microbiome (NKB bi-parental population)", 
           ggtheme=theme_bw)
grDevices::dev.set(1)
ggsave(filename="Metagenome20_corr_spearman.tiff", plot=plot, width=35, height= 35, dpi=300, compression = "lzw")


in_list <- commandArgs(trailingOnly = TRUE)
wd <- in_list[1]
input1 <- in_list[2]
taxa_percent <- in_list[3]



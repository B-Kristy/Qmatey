#Load Libraries
library(dplyr)
library(plotly)
library(htmlwidgets)
args <- commandArgs(TRUE)
dfm <- read.table("strain_taxainfo_mean.txt", header=T, sep="\t", check.names =FALSE, fill=TRUE)
dfu <- read.table("strain_taxainfo_unique_sequences.txt", header=T, sep="\t", check.names = FALSE, fill=TRUE)
dfe <- read.table("strain_taxainfo_quantification_accuracy.txt", header=T, sep="\t", check.names = FALSE, fill= TRUE)
percent_thresh <- 5

data.pipe<-function(dfm,dfu,dfe,cut){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  vDrop<-c("family","order", "tax_id","superkingdom","kingdom","class","species","genus")
  mirror<-dfm[,(which(colnames(dfm) %!in% vDrop))]
  mirror$phylum = as.character(mirror$phylum)
  mirror$taxname<-as.character(mirror$taxname)
 
   #set up empty arrays for the data values
  rowNum<-0
  taxid<-c()
  means<-c()
  uniq<-c()
  err<-c()
  phylum<-c()
  lines<-0
  #loop to input data into selected arrays
  for (row in 1:nrow(mirror) ) {
    lines=lines+1
    for (col in 1:(ncol(mirror)-2)) {
      rowNum<- rowNum + 1
      taxid[rowNum]<-mirror[row,ncol(mirror)-1]
      phylum[rowNum]<-mirror[row,ncol(mirror)]
      means[rowNum]<-mirror[row,col]
      uniq[rowNum]<-dfu[row,col+1]
      err[rowNum]<-dfe[row,col+1]
      
    }
    
  }
  
  #create dataframe to hold the data
  reform<-data.frame(matrix(nrow=rowNum, ncol = 5))

  #convert the arrays to vectors 
  taxid<-as.vector(unlist(taxid))
  means<-as.vector(unlist(means))
  uniq<-as.vector(unlist(uniq))
  err<-as.vector(unlist(err))
  phy<-as.vector(unlist(phylum))
  #set variables equal to corresponding vectors of data
  reform$X1<-as.factor(taxid)
  reform$X2<-means
  reform$X3<-uniq
  reform$X4<-err
  reform$X5<-phy
  #rename the variables
  colnames(reform)<-c("taxid","covMean","Unique Reads","errors","Phylum")
  #create total coverage data by multiplying the coverage by the unique reads
  reform$Tcov<-reform$covMean*reform$`Unique Reads`
  
  #cut data that makes up small percent of data
  met<-dfm
  met<-subset.data.frame(met, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
  met<-subset.data.frame(met, select=c(ncol(met),1:(ncol(met)-1)))
  met$count <-rowSums(met[,2:ncol(met)] > "0")
  met$percent <-((met$count)/(ncol(met)-2))*100
  bxD<-subset.data.frame(met, met$percent>= (5)/100)
  reform<-subset(reform, reform$taxid %in% bxD$taxname)
  reform[reform==0]<-NA
  
  return(reform)
}

reform <- data.pipe(dfm, dfu, dfe, cut)
df <- reform


uniq_box<-function(df, dfu){
  library(plotly)
  library(dplyr)
  #chop off taxa names
  uniq <- subset(dfu, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
  fence<-data.frame(matrix(nrow=nrow(uniq), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(uniq)) {
    taxa[pl]<-uniq[pl,ncol(uniq)]
    vec<-as.numeric(uniq[pl,2:(ncol(uniq)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  r<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box",marker=list(color = "transparent" ))%>%layout(title = "Boxplot of Unique Reads without Outliers",xaxis = list(range = c(0,xlim), title = "Unique Reads"), yaxis = list(size = 1, title = "Taxaname"))
  
  u<-plot_ly(reform, x = ~`Unique Reads`, y = ~taxid, color = ~Phylum, type = "box")%>%layout(title = "Boxplot of Unique Reads with Outliers", xaxis = list(range=c(0,xlim), title = "Unique Reads"), yaxis = list(size = 1, title = "Taxaname"))
  list(u,r)
  htmlwidgets::saveWidget(as_widget(r),"strain_level_unique_reads.html", selfcontained=FALSE)
  htmlwidgets::saveWidget(as_widget(u), "strain_level_outliers_unique_reads.html", selfcontained=FALSE)
}
uniq_box(df, dfu)

mean_box<-function(df, dfm){
  library(plotly)
  library(dplyr)
  #chop off taxa names
  mean <- subset(dfm, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
  fence<-data.frame(matrix(nrow=nrow(mean), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(mean)) {
    taxa[pl]<-mean[pl,ncol(mean)]
    vec<-as.numeric(mean[pl,2:(ncol(mean)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  r<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box",marker=list(color = "transparent" ))%>%layout(title = "Boxplot of Mean without Outliers",xaxis = list(range = c(0,xlim), title = "Mean"), yaxis = list(size = 1, title = "Taxaname"))
  
  u<-plot_ly(reform, x = ~`covMean`, y = ~taxid, color = ~Phylum, type = "box")%>%layout(title = "Boxplot of Mean Reads with Outliers", xaxis = list(range=c(0,xlim), title = "Mean Reads"), yaxis = list(size = 1, title = "Taxaname"))
  list(u,r)
  htmlwidgets::saveWidget(as_widget(r),"strain_level_mean_reads.html", selfcontained = FALSE)
  htmlwidgets::saveWidget(as_widget(u), "strain_level_outliers_mean_reads.html", selfcontained = FALSE)
}
mean_box(df, dfm)

error_box<-function(df, dfe){
  library(plotly)
  library(dplyr)
  #chop off taxa names
  error <- subset(dfe, select=-c(tax_id,species,genus,family,order,class,phylum,kingdom,superkingdom))
  fence<-data.frame(matrix(nrow=nrow(error), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(error)) {
    taxa[pl]<-error[pl,ncol(error)]
    vec<-as.numeric(error[pl,2:(ncol(error)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  r<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box",marker=list(color = "transparent" ))%>%layout(title = "Boxplot of Standard Error without Outliers",xaxis = list(range = c(0,xlim), title = "Standard Error"), yaxis = list(size = 1, title = "Taxaname"))
  
  u<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box")%>%layout(title= "Boxplot of Standard Error with Outliers", xaxis = list(titlerange=c(0,xlim), title = "Standard Error"), yaxis = list(size = 1, title = "Taxaname"))
  list(u,r)
  htmlwidgets::saveWidget(as_widget(r),"strain_level_errors.html", selfcontained = FALSE)
  htmlwidgets::saveWidget(as_widget(u), "strain_level_outliers_errors.html", selfcontained = FALSE)
}
error_box(df, dfe)
  

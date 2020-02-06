library(plyr)
library(data.table)
args <- commandArgs(TRUE)
sighits <- args[1]
minUniqRead <- args[2]

if (args[3] == "strain"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,11))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "species"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,13))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "genus"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,14))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "family"){
fileNames <- c(sighits) 
for (fileName in fileNames) {
stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
stats1 <- subset(stats1, select=c(1,15))
groups <- unique(stats1[,2])
stats <- NULL

for (i in 1:length(groups)) {
  outlier <- subset(stats1, stats1[,c(2)] == groups[i])
  OutVals <- boxplot(outlier, plot=FALSE)$out
  outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
  stats <- rbind(stats,outlier)
}
stats1 <- stats
stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
               N    = sum(!is.na(abundance)),
               mean = mean(abundance, na.rm=TRUE),
               sd   = sd(abundance, na.rm=TRUE),
               se   = sd / sqrt(N)
)
stats1 <- subset(stats1, select=c(1,3,2,5))
colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
stats1$Normerr <- (stats1$stderr/stats1$mean)
stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
stats1 <- subset(stats1, select=-c(Normerr))
stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
}

write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
#cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "order"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,16))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "class"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,17))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

if (args[3] == "phylum"){
  fileNames <- c(sighits) 
  for (fileName in fileNames) {
    stats1 <- read.table(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
    stats1 <- subset(stats1, select=c(1,18))
    groups <- unique(stats1[,2])
    stats <- NULL
    
    for (i in 1:length(groups)) {
      outlier <- subset(stats1, stats1[,c(2)] == groups[i])
      OutVals <- boxplot(outlier, plot=FALSE)$out
      outlier <- outlier[ ! outlier$abundance %in% OutVals, ]
      stats <- rbind(stats,outlier)
    }
    stats1 <- stats
    stats1 <- ddply(stats1, c(colnames(stats1[2])), summarise,
                    N    = sum(!is.na(abundance)),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$Normerr <- (stats1$stderr/stats1$mean)
    stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
    stats1 <- subset(stats1, select=-c(Normerr))
    stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
    stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
    stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
  }
  
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
  #cor(stats1, method="pearson", use="complete.obs")
}

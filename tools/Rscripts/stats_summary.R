library(plyr)
library(data.table)
args <- commandArgs(TRUE)
sighits <- args[1]
minUniqRead <- 1

fileNames <- c(sighits)
for (fileName in fileNames) {
stats1 <- read.table(file=fileName, header=T, sep="\t", quote="", check.names = FALSE)
stats1 <- subset(stats1, select=c(1,11))
groups <- unique(stats1[,2])
stats <- NULL

for (i in 1:length(groups)) {
  outlier <- subset(stats1, staxids == groups[i])
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
colnames(stats1) <- c("tax_id","mean","uniq_reads","stderr")
stats1$Normerr <- (stats1$stderr/stats1$mean)
stats1 <- subset(stats1, Normerr <= 0.2 | is.na(Normerr))
stats1 <- subset(stats1, select=-c(Normerr))
stats1$mean[stats1$uniq_reads < minUniqRead] <- NA
stats1$uniq_reads[stats1$uniq_reads < minUniqRead] <- NA
stats1$stderr[stats1$uniq_reads < minUniqRead] <- NA
}

write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE)
#cor(stats1, method="pearson", use="complete.obs")
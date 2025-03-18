#script to pre-process merged CpG files merged_CpG_evidence.cov (Bismark output)
#identify sites that agree on methylation status
#and export in a format that can subsequently be used in methylkit.

#arguments received when running using Rscript:
#1) rep 1 merged_CpG_evidence.cov to import
#2) rep 2 merged_CpG_evidence.cov to import
#3) rep 1 output as agree.cov file
#4) rep 2 output as agree.cov file
#5) rep 1 output as full_report.txt file (includes calculated 0/1 methylation status)
#6) rep 2 output as full_report.txt file (includes calculated 0/1 methylation status)
#7) stats file output named after both reps 

args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
setwd("/nobackup/beegfs/workspace/bh471/scripts/ma_imi/wgbs/")
repData1 <- read.delim(args[1], header=FALSE)
repData2 <- read.delim(args[2], header=FALSE)
#initialise stats
stats <- data.frame(matrix(ncol=2, nrow=0))
colnames(stats) <- c("description", "count")

#process file1
colnames(repData1) <- c("scaffold","start","end","percent","cReads","tReads")
repData1$cov <- repData1$cReads + repData1$tReads
repData1 <- subset(repData1, repData1$cov > 9 | cReads > 4)
stats <- rbind(stats, data.frame(status = "rep1 raw loci", count = nrow(repData1)))
repData1$loc <- paste0(repData1$scaffold,":",repData1$start)
bt <- function(a, b, p) {binom.test(a, b, p, alternative="greater")$p.value}
repData1$pVal <- mapply(bt, repData1$cReads, repData1$cov, 0.007)
repData1$FDR <- p.adjust(repData1$pVal, method = "BH", n = length(repData1$pVal))
repData1$status <- 0
repData1$status[repData1$FDR < 0.05] <- 1
stats <- rbind(stats, data.frame(status = "rep1 methylated loci", count = nrow(subset(repData1, repData1$status == 1))))

#process file2
colnames(repData2) <- c("scaffold","start","end","percent","cReads","tReads")
repData2$cov <- repData2$cReads + repData2$tReads
repData2 <- subset(repData2, repData2$cov > 9 | cReads > 4)
stats <- rbind(stats, data.frame(status = "rep2 raw loci", count = nrow(repData2)))
repData2$loc <- paste0(repData2$scaffold,":",repData2$start)
bt <- function(a, b, p) {binom.test(a, b, p, alternative="greater")$p.value}
repData2$pVal <- mapply(bt, repData2$cReads, repData2$cov, 0.007)
repData2$FDR <- p.adjust(repData2$pVal, method = "BH", n = length(repData2$pVal))
repData2$status <- 0
repData2$status[repData2$FDR < 0.05] <- 1
stats <- rbind(stats, data.frame(status = "rep2 methylated loci", count = nrow(subset(repData2, repData2$status == 1))))

#merge and identify sites that agree on status
mergedData <- merge(repData1, repData2, by = "loc")

colnames(mergedData) <- c("loc","scaffold1","start1", "end1", "percent1", "cReads1","tReads1", "cov1","pVal1", "FDR1", "status1", "scaffold2","start2", "end2", "percent2", "cReads2","tReads2", "cov2","pVal2", "FDR2", "status2")
mergedData$agree <- mergedData$status1 + mergedData$status2
agreeData <- subset(mergedData, mergedData$agree == 0 | mergedData$agree == 2)
stats <- rbind(stats, data.frame(status = "agreeing loci", count = nrow(agreeData)))

consensusRep1 <- select(merge(repData1, agreeData), scaffold,start,end,percent,cReads,tReads)
consensusRep2 <- select(merge(repData2, agreeData), scaffold,start,end,percent,cReads,tReads)
#output
write.table(consensusRep1, args[3], sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(consensusRep2, args[4], sep="\t", row.names = F, col.names = FALSE, quote = FALSE)

#full data inc status
fullRep1 <- select(merge(repData1, agreeData), scaffold,start,end,percent, cReads,tReads,status1)
fullRep2 <- select(merge(repData2, agreeData), scaffold,start,end,percent,cReads,tReads,status2)
write.table(fullRep1, args[5], sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(fullRep2, args[6], sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(stats, file=args[7], sep="\t", row.names = FALSE, quote = FALSE)
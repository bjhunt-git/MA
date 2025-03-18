#this script processes the output of Bismark to a format suitable for analysis with alphabeta
#argument supplied via Rscript is sample ID ie D01
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
setwd("/nobackup/beegfs/workspace/bh471/scripts/ma_imi/wgbs/")
inputPath="/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs/extract/"
outputPath <- "/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs/bismark_processed/"


data <- read.delim(paste0(inputPath, args[1], "_CpG/",args[1],"_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt"), header = FALSE, sep = "\t")
colnames(data) <- c("chromosome","position","strand","count_methylated","count_unmethylated","Ccontext","trinucleotide_context")
data$total <- data$count_methylated + data$count_unmethylated
data <- subset(data, total > 0)
data$fraction <- data$count_methylated / data$total
methData <- subset(data, count_methylated > 0)
unmethData <- subset(data, count_methylated == 0)

#deal with unmethylated data first
unmethData$status <- "U"
finalData <- unmethData

#binomial testing for methylated data
bt <- function(a, b, p = 0.007) {binom.test(a, b, 0.007, alternative="greater")$p.value}
methData$pVal <- mapply(bt, methData$count_methylated, methData$total)
methData$FDR <- p.adjust(methData$pVal, method = "BH", n = length(methData$pVal))
methData$status <- "I"
methData$status[methData$FDR < 0.05] <- "M"

finalData <- rbind(finalData, select(methData, chromosome, position, strand, count_methylated, count_unmethylated, Ccontext, trinucleotide_context, total, fraction, status))
finalData <- finalData[order(finalData$chromosome, finalData$position),]
finalData$posteriorMax = 1

outputData <- select(finalData, seqnames=chromosome, start=position, strand, context=Ccontext, counts.methylated=count_methylated, counts.total=total, posteriorMax, status, rc.meth.lvl=fraction, context.trinucleotide=trinucleotide_context)
write.table(outputData, file = paste0(outputPath, args[1],"_bismark.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

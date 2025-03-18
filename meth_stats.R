#this script generates methylation stats for the MA experiment per sample.

library("tidyverse")
library(data.table)

setwd("C:/Users/bh471/OneDrive - University of Exeter/projects/MA/data/cpg/diffmeth")
outdir1="C:/Users/bh471/OneDrive - University of Exeter/projects/MA/data/cpg/methylation_analysis"
#import the annotated loci data
lociDetail <- read.table("annotated_data_table_methylation.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")


#set up the start of a matrix of values. This will be added to, sample by sample, using cbind as we cycle through samples.
statsFile <- data.frame(description=character(), stringsAsFactors = FALSE)
statsFile[nrow(statsFile)+1,] <- c("No of loci analysed")
statsFile[nrow(statsFile)+1,] <- c("No of methylated loci")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated CpGs")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated reads")

wmpf <- data.frame(description=character(), stringsAsFactors = FALSE)
wmpf[nrow(wmpf)+1,] <- c("CDS")
wmpf[nrow(wmpf)+1,] <- c("intergenic")
wmpf[nrow(wmpf)+1,] <- c("intron")
wmpf[nrow(wmpf)+1,] <- c("promoter")

fracF <- data.frame(description=character(), stringsAsFactors = FALSE)
fracF[nrow(fracF)+1,] <- c("CDS")
fracF[nrow(fracF)+1,] <- c("intergenic")
fracF[nrow(fracF)+1,] <- c("intron")
fracF[nrow(fracF)+1,] <- c("promoter")


propF <- data.frame(description=character(), stringsAsFactors = FALSE)
propF[nrow(propF)+1,] <- c("CDS")
propF[nrow(propF)+1,] <- c("intergenic")
propF[nrow(propF)+1,] <- c("intron")
propF[nrow(propF)+1,] <- c("promoter")

wmpg <- select(lociDetail %>% count(gff_id1), gene=gff_id1)
wmpp <- wmpg
geneCat <- wmpg
fracG <- wmpg
fracP <- wmpg
pCat <- wmpg
catSum <- data.frame(category=character(), stringsAsFactors = FALSE)
catSum[nrow(catSum)+1,] <- c("High")
catSum[nrow(catSum)+1,] <- c("Low")
catSum[nrow(catSum)+1,] <- c("Medium")
catSum[nrow(catSum)+1,] <- c("Unmethylated")

pCatSum <- data.frame(category=character(), stringsAsFactors = FALSE)
pCatSum[nrow(pCatSum)+1,] <- c("High")
pCatSum[nrow(pCatSum)+1,] <- c("Low")
pCatSum[nrow(pCatSum)+1,] <- c("Medium")
pCatSum[nrow(pCatSum)+1,] <- c("Unmethylated")

#prep the import data
importList <- list.files(path=getwd(), pattern = ".CpG_report.merged_CpG_evidence.cov$")
#importList <- head(importList, 3)
rawDataNames <- gsub("_mergeCpG.CpG_report.merged_CpG_evidence.cov","",list.files(getwd(), pattern=".CpG_report.merged_CpG_evidence.cov$", full.names = FALSE), fixed = TRUE)
#rawDataNames <- head(rawDataNames,3)
for(i in seq_along(importList)){
  importFile <- importList[i]
  rawDataFileName <- rawDataNames[[i]]
  rawData <- read.delim(importFile, header = FALSE, stringsAsFactors = FALSE)
  colnames(rawData) <- c("scaffold","start","end","percent","cReads","tReads")
  rawData$id <- paste0(rawData$scaffold,":",rawData$start)
  rawData$total <- rawData$cReads + rawData$tReads
  #only look at reads that are able to be confidently called
  processData <- subset(rawData, total >= 10 | cReads >= 5)
  bt <- function(a, b, p) {binom.test(a, b, p, alternative="greater")$p.value}
  processData$pVal <- mapply(bt, processData$cReads, processData$total, 0.007)
  processData$FDR <- p.adjust(processData$pVal, method = "BH", n = length(processData$pVal))
  processData$status <- 0
  processData$status[processData$FDR < 0.05] <- 1
  processData$adjCReads <- processData$cReads 
  processData$adjCReads[processData$FDR > 0.05] <- 0
  processData$cCount <- 1
  pdAnnotated <- left_join(processData, lociDetail, by = "id")
  pdAnnotated <- select(pdAnnotated, scaffold, start, end, percent, cReads, tReads, id, total, pVal, FDR, status, adjCReads, cCount, feature, featureStart, featureEnd, strand, gff_id1, gff_id2, gff_id3, transcriptNo)

  fraction <- processData %>% summarise(sumCReads = sum(cReads), sumAdjCReads = sum(adjCReads), sumTReads = sum(tReads), sumTotal = sum(total), mCpgCount = sum(status), cpgCount = sum(cCount))

  tempStats <- data.frame(count=numeric(), stringsAsFactors = FALSE)
  tempStats[nrow(tempStats)+1,] <- format(nrow(processData), scientific = FALSE)
  tempStats[nrow(tempStats)+1,] <- format(nrow(subset(processData, status == 1)), scientific = FALSE)
  tempStats[nrow(tempStats)+1,] <- format(fraction$mCpgCount/fraction$cpgCount, scientific = FALSE)
  tempStats[nrow(tempStats)+1,] <- format(fraction$sumAdjCReads/fraction$sumTotal, scientific = FALSE)

  statsFile <- cbind(statsFile, tempStats)
  colnames(statsFile)[colnames(statsFile) == "count"] <- rawDataFileName

  featureLevel <- pdAnnotated %>% group_by(feature) %>% summarise(sumCReads = sum(cReads), sumAdjCReads = sum(adjCReads), sumTReads = sum(tReads), sumTotal = sum(total), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = sum(adjCReads)/sum(total))
  featureLevel$proportionMethLoci <- featureLevel$methCLoci / sum(featureLevel$methCLoci)
  featureLevel$fractionMethLoci <- featureLevel$methCLoci / featureLevel$totalCLoci

  wmpf <- cbind(wmpf, select(featureLevel, weightedMeth))
  colnames(wmpf)[colnames(wmpf) == "weightedMeth"] <- rawDataFileName

  fracF <- cbind(fracF, select(featureLevel,fractionMethLoci))
  colnames(fracF)[colnames(fracF) == "fractionMethLoci"] <- rawDataFileName

  propF <- cbind(propF, select(featureLevel,proportionMethLoci))
  colnames(propF)[colnames(propF) == "proportionMethLoci"] <- rawDataFileName

  #new addition to remove the promoters
  pdAnnotatedGene <- subset(pdAnnotated, feature != "promoter" & feature != "intergenic")
  pdAnnotatedPromoter <- subset(pdAnnotated, feature == "promoter")
  
  geneLevel <- pdAnnotatedGene %>% group_by(gff_id1) %>% summarise(sumCReads = sum(cReads), sumAdjCReads = sum(adjCReads), sumTReads = sum(tReads), sumTotal = sum(total), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = format(sum(adjCReads)/sum(total), scientific = FALSE))
  colnames(geneLevel)[colnames(geneLevel) == "gff_id1"] <- "gene"
  geneLevel$fractionMethLoci <- geneLevel$methCLoci / geneLevel$totalCLoci
  wmpg <- left_join(wmpg, select(geneLevel, gene, weightedMeth), by = "gene")
  colnames(wmpg)[colnames(wmpg) == "weightedMeth"] <- rawDataFileName
  fracG <- left_join(fracG, select(geneLevel, gene, fractionMethLoci), by = "gene")
  colnames(fracG)[colnames(fracG) == "fractionMethLoci"] <- rawDataFileName
  geneLevel$cat <- "U"
  geneLevel$cat[geneLevel$weightedMeth >= 0.05 & geneLevel$weightedMeth < 0.3] <- "L"
  geneLevel$cat[geneLevel$weightedMeth >= 0.3 & geneLevel$weightedMeth < 0.7] <- "M"
  geneLevel$cat[geneLevel$weightedMeth >= 0.7] <- "H"
  sampleCat <- geneLevel %>% count(cat)
  catSum <- cbind(catSum, select(sampleCat, n))
  colnames(catSum)[colnames(catSum) == "n"] <- rawDataFileName
  geneCat <- left_join(geneCat, select(geneLevel, gene, cat), by = "gene")
  colnames(geneCat)[colnames(geneCat) == "cat"] <- rawDataFileName
  
  promoterLevel <- pdAnnotatedPromoter %>% group_by(gff_id1) %>% summarise(sumCReads = sum(cReads), sumAdjCReads = sum(adjCReads), sumTReads = sum(tReads), sumTotal = sum(total), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = format(sum(adjCReads)/sum(total), scientific = FALSE))
  colnames(promoterLevel)[colnames(promoterLevel) == "gff_id1"] <- "gene"
  promoterLevel$fractionMethLoci <- promoterLevel$methCLoci / promoterLevel$totalCLoci
  wmpp <- left_join(wmpp, select(promoterLevel, gene, weightedMeth), by = "gene")
  colnames(wmpp)[colnames(wmpp) == "weightedMeth"] <- rawDataFileName
  fracP <- left_join(fracP, select(promoterLevel, gene, fractionMethLoci), by = "gene")
  colnames(fracP)[colnames(fracP) == "fractionMethLoci"] <- rawDataFileName
  promoterLevel$cat <- "U"
  promoterLevel$cat[promoterLevel$weightedMeth >= 0.05 & promoterLevel$weightedMeth < 0.3] <- "L"
  promoterLevel$cat[promoterLevel$weightedMeth >= 0.3 & promoterLevel$weightedMeth < 0.7] <- "M"
  promoterLevel$cat[promoterLevel$weightedMeth >= 0.7] <- "H"
  sampleCat <- promoterLevel %>% count(cat)
  pCatSum <- cbind(pCatSum, select(sampleCat, n))
  colnames(pCatSum)[colnames(pCatSum) == "n"] <- rawDataFileName
  pCat <- left_join(pCat, select(promoterLevel, gene, cat), by = "gene")
  colnames(pCat)[colnames(pCat) == "cat"] <- rawDataFileName
  

  write.table(geneLevel, file=paste0(outdir1,"/", rawDataFileName,"geneLevel.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
  write.table(promoterLevel, file=paste0(outdir1,"/", rawDataFileName,"promoterLevel.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
  write.table(featureLevel, file=paste0(outdir1,"/", rawDataFileName,"featureLevel.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
}


geneCat$count <- apply(geneCat, 1, function(x)length(unique(na.omit(x)))-1)
geneCat$countU <- apply(geneCat, 1, function(x) length(which(x=="U")))
geneCat$countL <- apply(geneCat, 1, function(x) length(which(x=="L")))
geneCat$countM <- apply(geneCat, 1, function(x) length(which(x=="M")))
geneCat$countH <- apply(geneCat, 1, function(x) length(which(x=="H")))

pCat$count <- apply(pCat, 1, function(x)length(unique(na.omit(x)))-1)
pCat$countU <- apply(pCat, 1, function(x) length(which(x=="U")))
pCat$countL <- apply(pCat, 1, function(x) length(which(x=="L")))
pCat$countM <- apply(pCat, 1, function(x) length(which(x=="M")))
pCat$countH <- apply(pCat, 1, function(x) length(which(x=="H")))

write.table(geneCat, file=paste0(outdir1,"/", "genes_methylation_status_categories.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(pCat, file=paste0(outdir1,"/", "promoter_methylation_status_categories.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(catSum, file=paste0(outdir1,"/", "genes_categories_summary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(pCatSum, file=paste0(outdir1,"/", "promoters_categories_summary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(fracF, file=paste0(outdir1,"/", "feature_fractional_methylation.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(fracG, file=paste0(outdir1,"/", "genes_fractional_methylation.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(fracP, file=paste0(outdir1,"/", "promoters_fractional_methylation.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(propF, file=paste0(outdir1,"/", "features_proportion_methylation.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(statsFile, file=paste0(outdir1,"/", "basic_stats.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(wmpf, file=paste0(outdir1,"/", "weighted_meth_per_feature.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(wmpg, file=paste0(outdir1,"/", "weighted_meth_per_gene.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")
write.table(wmpp, file=paste0(outdir1,"/", "weighted_meth_per_promoter.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep= "\t")

#script to identify differential methylated sites using Fisher's exact test, with annotated output.
#requires input from preprocess_merged_cov.R, which identifies the methylated sites per sample, including only sites where the replicates agree.
#set scan for args
#arguments received
#1) Rep 1 name
#2) Rep 2 name
#3) group (eg line 2 gen 51 unexposed, in this analysis that would be L2_G51_UNSEL)
args = commandArgs(trailingOnly=TRUE)
print(args[1])
print(args[2])
print(args[3])
#import libraries and sort working directories
input.path="/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs/preprocess"
output.path="/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs/diffmeth/"
print(input.path)
print(output.path)

setwd(input.path)
library(methylKit)
library(GenomicFeatures)
library(Rsamtools)
library(dplyr)

#process the 'full_report.txt' files output by preprocess_merged_cov.R to identify sites with methylation in at least one sample, then output new CpG.txt files for each one. 

fullReport1 <- read.delim(paste0(args[1],"_full_report.txt"), header=FALSE)
colnames(fullReport1) <- c("scaffold","start","end","fraction", "cReads","tReads", "status")
fullReport2 <- read.delim(paste0(args[2],"_full_report.txt"), header=FALSE)
colnames(fullReport2) <- c("scaffold","start","end","fraction", "cReads","tReads", "status")
fullReport3 <- read.delim("D59_full_report.txt", header=FALSE)
colnames(fullReport3) <- c("scaffold","start","end","fraction", "cReads","tReads", "status")
fullReport4 <- read.delim("D60_full_report.txt", header=FALSE)
colnames(fullReport4) <- c("scaffold","start","end","fraction", "cReads","tReads", "status")

#create loc field
fullReport1$loc <- paste0(fullReport1$scaffold,":",fullReport1$start)
fullReport2$loc <- paste0(fullReport2$scaffold,":",fullReport2$start)
fullReport3$loc <- paste0(fullReport3$scaffold,":",fullReport3$start)
fullReport4$loc <- paste0(fullReport4$scaffold,":",fullReport4$start)

m1 <- select(fullReport1, loc, status1 = status)
m2 <- select(fullReport2, loc, status2 = status)
m3 <- select(fullReport3, loc, status3 = status)
m4 <- select(fullReport4, loc, status4 = status)

mergedData <- merge(m1, m2, by = "loc")
mergedData <- merge(mergedData, m3, by = "loc")
mergedData <- merge(mergedData, m4, by = "loc")

mergedData$minOneMeth <- mergedData$status1 + mergedData$status2 + mergedData$status3 + mergedData$status4

methLocs <- subset(mergedData, minOneMeth > 0)
unMethLocs <- subset(mergedData, minOneMeth == 0)

cpgReport1 <- select(merge(fullReport1, methLocs, by = "loc"), scaffold,start,end,fraction,cReads,tReads)
cpgReport2 <- select(merge(fullReport2, methLocs, by = "loc"), scaffold,start,end,fraction,cReads,tReads)
cpgReport3 <- select(merge(fullReport3, methLocs, by = "loc"), scaffold,start,end,fraction,cReads,tReads)
cpgReport4 <- select(merge(fullReport4, methLocs, by = "loc"), scaffold,start,end,fraction,cReads,tReads)

dir.create(paste0(args[3],"_methCovFiles"))
write.table(cpgReport1, paste0(args[3],"_methCovFiles","/",args[1],".CpG.final.cov"), sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(cpgReport2, paste0(args[3],"_methCovFiles","/",args[2],".CpG.final.cov"), sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(cpgReport3, paste0(args[3],"_methCovFiles","/D59.CpG.final.cov"), sep="\t", row.names = F, col.names = FALSE, quote = FALSE)
write.table(cpgReport4, paste0(args[3],"_methCovFiles","/D60.CpG.final.cov"), sep="\t", row.names = F, col.names = FALSE, quote = FALSE)

#----differential methylation calculation----
#import the relevant sample files and the ancestral files to make the comparison for that line

importFileList <- as.list(c(paste0(args[3],"_methCovFiles","/D59.CpG.final.cov"),paste0(args[3],"_methCovFiles","/D60.CpG.final.cov"),paste0(args[3],"_methCovFiles","/",args[1],".CpG.final.cov"),paste0(args[3],"_methCovFiles","/",args[2],".CpG.final.cov")))
importSampleList <- as.list(c("D59","D60",args[1],args[2]))
rawDB <- methRead(importFileList, sample.id=importSampleList,
                  assembly="NSSCAFFOLD",
                  treatment=c(0,0,1,1),
                  context="CpG",
                  mincov=5,
                  pipeline = "bismarkCoverage")

filterDB <- filterByCoverage(rawDB,lo.count = 5,
                             lo.perc = NULL,
                             hi.count = 200,
                             hi.perc = NULL)

normDB = normalizeCoverage(filterDB,method="median")

uniteDB=methylKit::unite(normDB)
poolDB=methylKit::pool(uniteDB, sample.ids=c("ancestor","descendant"))

#### loci level ####
diff=calculateDiffMeth(poolDB,overdispersion="MN")
sigDiff=getMethylDiff(diff,difference=25,qvalue=0.01)
lociAnnotation <- read.delim("annotated_data_table_methylation.txt", header = TRUE, sep = "\t")
annotation <- read.delim("ns_jg_blast_annotation.txt", header = TRUE, sep = "\t")
annotation <- dplyr::select(annotation, gene, transcript=SeqName, Description, evalue=e.Value, GOterms=GO.IDs, GONames=GO.Names, enzymeCodes=Enzyme.Codes, enzymeNames = Enzyme.Names, InterProIDs=InterPro.IDs, InterProGOIDs=InterPro.GO.IDs, InterProNames=InterPro.GO.Names)
sigDiffLoci <- getData(sigDiff)
sigDiffLoci$id <- paste0(sigDiffLoci$chr,":",sigDiffLoci$start)
#genes
sigDiffLociAnnotated <- merge(sigDiffLoci, lociAnnotation, by = "id")
sigDiffLociAnnotated$gene <- sigDiffLociAnnotated$gff_id1
sigDiffLociAnnotated <- merge(sigDiffLociAnnotated, annotation, by = "gene", all.x = TRUE)
sigDiffLociAnnotated <- dplyr::select(sigDiffLociAnnotated, id, gene, transcript, pvalue, qvalue, meth.diff, scaffold=chr, start, end, strand=strand.x, feature, featureStart, featureEnd, Description, evalue, GOterms, GONames, enzymeCodes, InterProIDs, InterProGOIDs, InterProNames)
write.table(sigDiffLociAnnotated, file=paste(output.path,args[3],"_all_sigdiff_loci_annotated.txt", sep=""), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#final output of basic numbers
#process the 'full_report.txt' files to identify sites with methylation in at least one sample, then output new CpG.txt files for each one. 
#final output of basic numbers
stats <- data.frame(matrix(ncol=2, nrow=0))
colnames(stats) <- c("description", "count")
stats <- rbind(stats, data.frame(status = "All covered loci", count = nrow(mergedData)))
stats <- rbind(stats, data.frame(status = "Loci methylated in at least one replicate", count = nrow(methLocs)))
stats <- rbind(stats, data.frame(status = "Loci with no methylation in any replicate", count = nrow(unMethLocs)))
stats <- rbind(stats, data.frame(status = "Number of significant methylation differences at 25% change", count = nrow(sigDiff)))
stats <- rbind(stats, data.frame(status = "Hypermethylated differences", count = nrow(subset(sigDiffLoci, meth.diff >0))))
stats <- rbind(stats, data.frame(status = "Hypomethylated differences", count = nrow(subset(sigDiffLoci, meth.diff <0))))
write.table(stats, file=paste(output.path,args[3],"_diff_meth_stats.txt", sep=""), sep="\t", row.names = FALSE)

#scaffold stats
scaffoldCount <-count(mergedData, chromosome)
colnames(scaffoldCount) <- c("scaffold", "scaffoldCount")

diffCount <-count(sigDiffLoci, chr)
colnames(diffCount) <- c("scaffold", "diffCount")
stats2 <- merge(scaffoldCount, diffCount, by = "scaffold", all.x = TRUE)
stats2 <- select(tidyr::separate(stats2, scaffold, c("HiC","scaf","number"), sep="_", remove=FALSE,convert=TRUE), scaffold, scaffoldCount,diffCount, number)
stats2 <- stats2[order(stats2$number),]
write.table(stats2, file=paste(output.path,args[3],"_scaffold_stats.txt", sep=""), sep="\t", row.names = FALSE)

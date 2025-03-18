library(tidyverse)

inputPath="/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/rnaseq/stringtie_new_annotation/"
#following generated in terminal with:
#awk '$3 == "gene"' ns.gff3 | sed 's/ID=//g' | cut -f9 > ns.jg.genelist.txt
geneList <- read.delim("/nobackup/beegfs/workspace/bh471/results/static_data/ns_annotation/ns_jg_softmask_bams/ns.jg.genelist.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(geneList) <- c("gene")
covMatrix <- geneList
fpkmMatrix <- geneList
tpmMatrix <- geneList

importList <- as.list(list.files(inputPath,full.names=T,recursive=T,pattern=".abund.txt"))
samples <- as.data.frame(gsub(".abund.txt","",list.files(inputPath, recursive=T,pattern=".abund.txt$", full.names = FALSE), fixed = TRUE))
colnames(samples) <- "sample"
samples <- select(separate(samples, sample, c("sample","T2"), sep = "\\/", fill = "right", convert=TRUE, remove = TRUE), sample)

for(i in 1:nrow(samples)){
  importFile <- importList[[i]]
  sampleName <- samples[i,]
  tempImport <- read.delim(importFile, header = TRUE, stringsAsFactors = FALSE)
  colnames(tempImport) <- c("gene","name","reference","strand","start", "end", "coverage","FPKM", "TPM")
  covMatrix <- left_join(covMatrix, select(tempImport, gene, coverage), by = "gene")
  colnames(covMatrix)[colnames(covMatrix) == "coverage"] <- sampleName
  fpkmMatrix <- left_join(fpkmMatrix, select(tempImport, gene, FPKM), by = "gene")
  colnames(fpkmMatrix)[colnames(fpkmMatrix) == "FPKM"] <- sampleName
  tpmMatrix <- left_join(tpmMatrix, select(tempImport, gene, TPM), by = "gene")
  colnames(tpmMatrix)[colnames(tpmMatrix) == "TPM"] <- sampleName
  
}
write.table(covMatrix, file=paste0(inputPath,"/","coverage.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(fpkmMatrix, file=paste0(inputPath,"/","FPKM.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tpmMatrix, file=paste0(inputPath,"/","TPM.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
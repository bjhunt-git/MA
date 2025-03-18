#script for annotating individual loci after intersection of loci with bedtools
#input file generated using annotate_loci.sh 

library("tidyverse")
library(data.table)

setwd("C:/Users/bh471/OneDrive - University of Exeter/projects/MA/data/cpg/diffmeth")

annotated <- read.delim("C:/Users/bh471/OneDrive - University of Exeter/projects/MA/data/cpg/diffmeth/annotated.split.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(annotated) <- c("scaffold","bedStart","start","featureScaffold","source","feature","featureStart","featureEnd","score","strand","dot","gff_id1", "gff_id2", "gff_id3")
annotated$end <- annotated$start+1
annotated <- select(annotated, scaffold, start, end, feature, featureScaffold,featureStart, featureEnd, strand, gff_id1, gff_id2, gff_id3)
annotated <- annotated[order( annotated[,1], annotated[,2]),]
annotated$id <- paste0(annotated$scaffold,":",annotated$start)

#split annotated into genes, transcripts and features.
genes <- subset(annotated, feature == "gene" )
transcripts <- subset(annotated, feature == "mRNA" )
features <- subset(annotated, feature != "mRNA" & feature != "gene" & feature != "promoter")
promoters <- subset(annotated, feature == "promoter" )
#nrow(features) + nrow(genes) + nrow(transcripts) + nrow(promoters)
nrow(annotated) == nrow(features) + nrow(genes) + nrow(transcripts) + nrow(promoters)


#generate a final list of sites with the best annotation for each. 
#create the table that will accumulate annotated sites

#prepare variables and columns used for looping through the maximum number of transcripts a gene might have, in order to sequentially annotate catching those that are annotated only by, eg transcript 2. 
features <- separate(features, gff_id2, c("T1","T2"), sep = "\\.", fill = "right", convert=TRUE, remove = FALSE)
features$transcriptNo <- as.integer(gsub("t","",features$T2))
features <- select(features, scaffold,start,end,feature,featureScaffold,featureStart,featureEnd,strand,gff_id1,gff_id2,gff_id3,id,transcriptNo)
#colnames(features)[9] <- "gene_id"
annList <- split(features, f = features$transcriptNo)
loopMax <- max(features$transcriptNo)


finalAnnotation <- data.frame(scaffold=character(),start=integer(),end=integer(),feature=character(),featureScaffold=character(),featureStart=integer(),featureEnd=integer(),strand=character(),gff_id1=character(),gff_id2=character(),gff_id3=character(), id=character(),transcriptNo=character()) 

for (i in 1:loopMax){
  print(i)
  #bottom of the feature hierarch first - tts, tss, start and stop codons
  level1 <- subset(annList[[i]], feature == "transcription_end_site" | feature == "transcription_start_site" | feature == "start_codon" | feature == "stop_codon")
  #countlevel1ID <- level1 %>% group_by(id) %>% tally 
  #the line below removes duplicates based only on id, keeping the first instance and discarding the rest
  level1Uniq <- level1[!duplicated(level1[,c('id')]),] #number?
  #add to the final annotation
  newInLevel1 <- anti_join(level1Uniq, finalAnnotation, by = c("id"))
  finalAnnotation <- rbind(finalAnnotation, newInLevel1)
  
  #next level - UTRs and CDS.
  level2 <- subset(annList[[i]], feature == "five_prime_utr" | feature == "three_prime_utr" | feature == "CDS")
  #countlevel2ID <- level2 %>% group_by(id) %>% tally 
  #the line below removes duplicates based only on id, keeping the first instance and discarding the rest
  level2Uniq <- level2[!duplicated(level2[,c('id')]),] #number?
  #add to the final annotation
  newInLevel2 <- anti_join(level2Uniq, finalAnnotation, by = c("id"))
  finalAnnotation <- rbind(finalAnnotation, newInLevel2)
  
  #level 3 - exons and introns
  level3 <- subset(annList[[i]], feature == "exon" | feature == "intron")
  #countlevel3ID <- level3 %>% group_by(id) %>% tally 
  #the line below removes duplicates based only on id, keeping the first instance and discarding the rest
  level3Uniq <- level3[!duplicated(level3[,c('id')]),] #number?
  #add to the final annotation
  newInLevel3 <- anti_join(level3Uniq, finalAnnotation, by = c("id"))
  finalAnnotation <- rbind(finalAnnotation, newInLevel3)
  
}

#level 4 - promoters and downstream
level4 <- subset(annotated, feature == "promoter" | feature == "downstream")
#colnames(level4)[9] <- "gene_id"
level4$transcriptNo <- NA
level4$gff_id1 <- level4$gff_id3
#level4 <- select(level4, "scaffold","start","end","feature","featureScaffold","featureStart","featureEnd","strand","gene_id","T2","gff_id1","id", "transcriptNo")



#countlevel4ID <- level4 %>% group_by(id) %>% tally 
#the line below removes duplicates based only on id, keeping the first instance and discarding the rest
level4Uniq <- level4[!duplicated(level4[,c('id')]),] #number?
#add to the final annotation
newInLevel4 <- anti_join(level4Uniq, finalAnnotation, by = c("id"))
finalAnnotation <- rbind(finalAnnotation, newInLevel4)


#finally, the sites that can only be annotated by gene id (shouldn't be any but just in case)
level5 <- subset(annotated, feature == "gene")
#colnames(level5)[9] <- "gene_id"
level5$transcriptNo <- NA
#level5 <- select(level5, "scaffold","start","end","feature","featureScaffold","featureStart","featureEnd","strand","gene_id","T2","gff_id1","id", "transcriptNo")


#countlevel5ID <- level5 %>% group_by(id) %>% tally 
#the line below removes duplicates based only on id, keeping the first instance and discarding the rest
level5Uniq <- level5[!duplicated(level5[,c('id')]),] #number?
#add to the final annotation
newInLevel5 <- anti_join(level5Uniq, finalAnnotation, by = c("id"))
finalAnnotation <- rbind(finalAnnotation, newInLevel5)

coveredSites <- read.delim("C:/Users/bh471/OneDrive - University of Exeter/projects/MA/data/cpg/diffmeth/sorted.loci.bed", header=FALSE, stringsAsFactors = FALSE)
colnames(coveredSites) <- c("siteScaffold","bedStart","siteStart")
coveredSites$siteEnd <- coveredSites$siteStart + 1
coveredSites <- select(coveredSites, siteScaffold, siteStart, siteEnd)
coveredSites$id <- paste0(coveredSites$siteScaffold,":",coveredSites$siteStart)

#left to do - gene start and end for wmpg (can do a separate table?)
finalData <- left_join(coveredSites, finalAnnotation, by = "id")
finalData <- select(finalData, siteScaffold, siteStart, siteEnd, id, feature, featureStart, featureEnd, strand, gff_id1, gff_id2, gff_id3,transcriptNo)
finalData$feature[is.na(finalData$feature)] <- "intergenic"
write.table(finalData, file="annotated_data_table_methylation.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
genesAndFeatures <- select(subset(finalData, feature == "intron" | feature == "exon"), gene_id, feature, gff_id1)

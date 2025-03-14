#script for conducting and visualising DMR cluster analysis in 
#the NS clone of M. persicae. 

library(tidyverse) #general wrangling
library(data.table) #also wrangling and processing the input data
library(RColorBrewer) #plotting

wd <- "datasets"
out_path <- "datasets"
#load the data - all differences from generation 51. 
import_list <- list.files(path=wd, pattern = "G51.+_annotated.txt$", recursive = TRUE, full.names = TRUE)
diffs_data <- lapply(import_list, function(file) {
  # Read the file with headers
  print(file)
  df <- read.delim(file, header = TRUE, fill = TRUE)
  df$sample <- gsub("projects/MA/data/cpg/diffmeth/diff_files/","",gsub("_all_sigdiff_loci_annotated.txt","",file, fixed = TRUE))
  return(df)
})

names(diffs_data) <- gsub("_all_sigdiff_loci_annotated.txt","",list.files(path=wd, pattern="G51.+_annotated.txt$", recursive = TRUE, include.dirs = FALSE, full.names = FALSE), fixed = TRUE)

#exposed and unexposed numbers
unexp <- diffs_data[grepl("_UNSEL", names(diffs_data))]
exp <- diffs_data[grepl("_SEL", names(diffs_data))]

exp_combo <- bind_rows(exp) #84827 differences total in gen 51 across all lines
exp_combo <- select(exp_combo, id, gene, transcript, qvalue, meth.diff, scaffold, start, end)
exp_uniq <- exp_combo[!duplicated(exp_combo[,c('id')]),] #42768 unique sites = 50.4% of total is unique

unexp_combo <- bind_rows(unexp) #129702 differences total in gen 51 across all lines
unexp_combo <- select(unexp_combo, id, gene, transcript, qvalue, meth.diff, scaffold, start, end)
unexp_uniq <- unexp_combo[!duplicated(unexp_combo[,c('id')]),] #54770 unique sites = 42% of total is unique

#get a unique list of differentially methylated sites, randomly select 20k

diffs_combo <- bind_rows(diffs_data) #214529 differences total in gen 51 across all lines
diffs_combo <- select(diffs_combo, id, gene, transcript, qvalue, meth.diff, scaffold, start, end)
diffs_uniq <- diffs_combo[!duplicated(diffs_combo[,c('id')]),] #66333 unique sites
write.table(diffs_uniq,file=paste0(out_path, "/g51_diffs_unique.txt"))
analysis_diffs <- sample(diffs_uniq$id, 10000, replace = FALSE)
write.table(analysis_diffs, file=paste0(out_path, "/analysis_diffs.txt"))
#get the methylation level data for all sites - generated from Bismark output using bedtools_union.sh
union <- read.table(file=paste0(wd,"/union.na.omit.bedGraph"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(union) <- c("scaffold", "start", "end", "D01","D02","D03","D04","D05","D06","D07","D08","D09","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D27","D28","D29","D30","D31","D32","D33","D34","D35","D36","D37","D38","D39","D40","D41","D42","D43","D44","D45","D46","D47","D48","D49","D50","D51","D52","D53","D54","D55","D56","D57","D58","D59","D60")
union$id <- paste0(union$scaffold, ":", union$end)

union_diffs <- union[union$id %in% analysis_diffs,]
union_diffs$id <- paste0(union_diffs$scaffold, ":", union_diffs$end)
union_diffs_heatmap <- select(union_diffs, D01:D56, D59:D60)
colnames(union_diffs_heatmap) <- c("Unexposed line 1 gen 51","Unexposed line 1 gen 51","Unexposed line 1 gen 52","Unexposed line 1 gen 52","Unexposed line 2 gen 51","Unexposed line 2 gen 51","Unexposed line 2 gen 52","Unexposed line 2 gen 52","Unexposed line 3 gen 51","Unexposed line 3 gen 51","Unexposed line 3 gen 52","Unexposed line 3 gen 52","Unexposed line 4 gen 51","Unexposed line 4 gen 51","Unexposed line 4 gen 52","Unexposed line 4 gen 52","Unexposed line 5 gen 51","Unexposed line 5 gen 51","Unexposed line 5 gen 52","Unexposed line 5 gen 52","Unexposed line 5 gen 51","Unexposed line 5 gen 51","Unexposed line 5 gen 52","Unexposed line 5 gen 52","Unexposed line 7 gen 51","Unexposed line 7 gen 51","Unexposed line 7 gen 52","Unexposed line 7 gen 52","Unexposed line 8 gen 51","Unexposed line 8 gen 51","Unexposed line 8 gen 52","Unexposed line 8 gen 52","Exposed line 2 gen 51","Exposed line 2 gen 51","Exposed line 2 gen 52","Exposed line 2 gen 52","Exposed line 4 gen 51","Exposed line 4 gen 51","Exposed line 4 gen 52","Exposed line 4 gen 52","Exposed line 5 gen 51","Exposed line 5 gen 51","Exposed line 5 gen 52","Exposed line 5 gen 52","Exposed line 6 gen 51","Exposed line 6 gen 51","Exposed line 6 gen 52","Exposed line 6 gen 52","Exposed line 7 gen 51","Exposed line 7 gen 51","Exposed line 7 gen 52","Exposed line 7 gen 52","Exposed line 8 gen 51","Exposed line 8 gen 51","Exposed line 8 gen 52","Exposed line 8 gen 52","Ancestral G2","Ancestral G2")
#union_diffs_heatmap_long <- pivot_longer(union_diffs_heatmap, cols = starts_with("D"), names_to = "sample", values_to = "methylation")
#meth diffs plot
pdf(file=paste0(out_path, "/diffs_heatmap.pdf"), width = 10, height = 15)
heatmap(
  as.matrix(union_diffs_heatmap), Rowv=NULL,scale = "row",
  col = heat.colors(12),
  Colv=as.dendrogram(hclust(dist(t(as.matrix(union_diffs_heatmap)))))
)
dev.off()
legend(x="topright", legend=c("min", "ave", "max"), fill=heat.colors(12))

#looking at the persistence of diffs
#import again just in case you've tinkered with things
import_list <- list.files(path=wd, pattern = "G51.+_annotated.txt$", recursive = TRUE, full.names = TRUE)
diffs_data <- lapply(import_list, function(file) {
  # Read the file with headers
  print(file)
  df <- read.delim(file, header = TRUE, fill = TRUE)
  df$sample <- gsub("projects/MA/data/cpg/diffmeth/diff_files/","",gsub("_all_sigdiff_loci_annotated.txt","",file, fixed = TRUE))
  df$group <- substr(df$sample, 1,6)
  df$line <- substr(df$group, 1,2)
  return(df)
})
names(diffs_data) <- gsub("_all_sigdiff_loci_annotated.txt","",list.files(path=wd, pattern="G51.+_annotated.txt$", recursive = TRUE, include.dirs = FALSE, full.names = FALSE), fixed = TRUE)

#and generation 52
import_list <- list.files(path=wd, pattern = "G52.+_annotated.txt$", recursive = TRUE, full.names = TRUE)
diffs_data_52 <- lapply(import_list, function(file) {
  # Read the file with headers
  print(file)
  df <- read.delim(file, header = TRUE, fill = TRUE)
  df$sample <- gsub("projects/MA/data/cpg/diffmeth/diff_files/","",gsub("_all_sigdiff_loci_annotated.txt","",file, fixed = TRUE))
  df$group <- substr(df$sample, 1,6)
  df$line <- substr(df$group, 1,2)
  return(df)
})
names(diffs_data_52) <- gsub("_all_sigdiff_loci_annotated.txt","",list.files(path=wd, pattern="G52.+_annotated.txt$", recursive = TRUE, include.dirs = FALSE, full.names = FALSE), fixed = TRUE)

diffs_combo_52 <- bind_rows(diffs_data_52)
diffs_combo_52 <- select(diffs_combo_52, id, gene, transcript, qvalue, meth.diff, scaffold, start, end)
diffs_uniq_52 <- diffs_combo_52[!duplicated(diffs_combo_52[,c('id')]),] #66333 unique sites
write.table(diffs_uniq_52,file=paste0(out_path, "/g52_diffs_unique.txt"))


combo51 <- bind_rows(diffs_data)
combo51 <- select(combo51, id, gene, transcript, qvalue, meth.diff, scaffold, start, end, feature, sample, group, line)
combo52 <- bind_rows(diffs_data_52)
combo52 <- select(combo52, id, gene, transcript, qvalue, meth.diff, scaffold, start, end, feature, sample, group, line)
both_gens <- rbind(combo51, combo52)

#number of sites that are intergenic
nrow(subset(combo51, feature == "intergenic")) #69071 
nrow(subset(combo52, feature == "intergenic")) #70284

#do this separately to avoid confusion over line number
unexp <- subset(both_gens, grepl("_UNSEL", sample))
lines <- unexp %>% distinct(line)

for (x in lines$line) {
  print(x)
  g51 <- subset(subset(unexp, grepl("_G51", sample)), line == x)
  g52 <- subset(subset(unexp, grepl("_G52", sample)), line == x)
  joined <- left_join(g51, g52, by = "id")
  write.table(joined, file=paste0(wd,"/unexposed_",x,"gen_compare.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote =FALSE)
}

exp <- subset(both_gens, grepl("_SEL", sample))
lines <- exp %>% distinct(line)

for (x in lines$line) {
  print(x)
  g51 <- subset(subset(exp, grepl("_G51", sample)), line == x)
  g52 <- subset(subset(exp, grepl("_G52", sample)), line == x)
  joined <- left_join(g51, g52, by = "id")
  write.table(joined, file=paste0(wd,"/exposed_",x,"gen_compare.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote =FALSE)
}

genesunexp <- (unexp %>% distinct(unexp$gene))
colnames(genesunexp) <- c("gene")
genesexp <- exp %>% distinct(exp$gene)
colnames(genesexp) <- c("gene")
genes <- rbind(genesunexp, genesexp)               
genes <- genes %>% distinct(genes$gene)
write.table(genes, file=paste0(wd,"/gen51_genes.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote =FALSE)

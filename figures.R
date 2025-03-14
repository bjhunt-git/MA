#code for generation of figures in "Sub-lethal insecticide stress alters epimutation rate
#but not genetic mutation rate in the pest insect Myzus persicae"
library(ggplot2)
library(circlize)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(FSA)

wd <- "datasets"

#(epi)mutation rates
dnadata <- read.delim(paste0(wd,"/rates_dna_final.txt"), stringsAsFactors = FALSE)
dnadata$Treatment <- factor(dnadata$Treatment , levels=c("Unexposed", "Exposed"))
dnadata$RateBoost <- dnadata$Rate * 10000000000
ggplot(dnadata, aes(x=Treatment, y = RateBoost, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 2, alpha = 1) +
  scale_fill_manual(values = c("cornflowerblue","firebrick2")) +
  scale_color_manual(values = c("cornflowerblue","firebrick2")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_minimal() +
  labs(x = "Treatment", y = "Rate (e-10)") +
  theme(axis.text.x = element_text(size = 12, family=("serif")),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, family=("serif")),
        legend.position = "")

epidata <- read.delim(paste0(wd,"/rates_epi_final.txt"))
epidata$Treatment <- factor(epidata$Treatment , levels=c("Unexposed", "Exposed"))
epidata$RateBoost <- epidata$Rate * 100000
ggplot(epidata, aes(x=Treatment, y = RateBoost, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 2, alpha = 1) +
  scale_fill_manual(values = c("cornflowerblue","firebrick2")) +
  scale_color_manual(values = c("cornflowerblue","firebrick2")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_minimal() +
  labs(x = "Treatment", y = "Rate (e-5)") +
  theme(axis.text.x = element_text(size = 12, family=("serif")),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, family=("serif")),
        legend.position = "")

#circle plot
circos.clear()
genome.file = read.delim(paste0(wd,"/ns_chr.txt"), sep = "\t", header = FALSE)
circos.genomicInitialize(genome.file, sector.names = c("1", "2", "3","4","5","6"), major.by = 25000000, axis.labels.cex = 0.7, labels.cex = 0.001, tickLabelsStartFromZero = TRUE)
circos.track(ylim = c(0, 1), 
             bg.col = c("pink", "lightgreen", "violet","skyblue","grey","orange"), 
             bg.border = NA, track.height = 0.05)
bed <- read.delim(paste0(wd,"/circos_control_final_snps.txt"), sep = "\t", header = FALSE)
bed$V4 <- bed$V4+0.25
circos.genomicTrack(bed, track.height = 0.05, ylim = c(0.25, 0.75),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col = "cornflowerblue")
                    })
bed2 <- read.delim(paste0(wd,"/circos_exposed_final_snps.txt"), sep = "\t", header = FALSE)
bed2$V4 <- bed2$V4+0.25
circos.genomicTrack(bed2, track.height = 0.05, ylim = c(0.25, 0.75),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col = "firebrick2")
                    })
controlDiffs <- read_delim(paste0(wd,"/control_diffs_circos.txt"), 
                           delim = "\t", escape_double = FALSE, 
                           col_types = cols(start = col_double(), 
                           end = col_double(), meth.diff = col_double()), 
                           trim_ws = TRUE)
controlDiffs <- subset(controlDiffs, scaffold == "HiC_scaffold_1" | scaffold == "HiC_scaffold_2" | scaffold == "HiC_scaffold_3" | scaffold == "HiC_scaffold_4" | scaffold == "HiC_scaffold_5" | scaffold == "HiC_scaffold_6")
exposedDiffs <- read_delim(paste0(wd,"/exposed_diffs_circos.txt"), 
                           delim = "\t", escape_double = FALSE, 
                           col_types = cols(start = col_double(), 
                           end = col_double(), meth.diff = col_double()), 
                           trim_ws = TRUE)
exposedDiffs <- subset(exposedDiffs, scaffold == "HiC_scaffold_1" | scaffold == "HiC_scaffold_2" | scaffold == "HiC_scaffold_3" | scaffold == "HiC_scaffold_4" | scaffold == "HiC_scaffold_5" | scaffold == "HiC_scaffold_6")
circos.genomicDensity(controlDiffs, col = c("cornflowerblue"), track.height = 0.1)
circos.genomicDensity(exposedDiffs, col = c("firebrick2"), track.height = 0.1)


text(0, 0.05, "Unexposed", col = "cornflowerblue", cex = 1.5, family = "serif")
text(0, -0.05, "Exposed", col = "firebrick2", cex = 1.5, family = "serif")

#fitness

fitness <- read.delim(paste0(wd,"/total_offspring.txt"), stringsAsFactors = FALSE) #total offspring
fitness <- read.delim(paste0(wd,"/days_til_first.txt"), stringsAsFactors = FALSE) #days until first nymph
fitness <- read.delim(paste0(wd,"/days_repro.txt"), stringsAsFactors = FALSE) #reproductive days
fitness$Treatment <- factor(fitness$Treatment , levels=c("Ancestral", "Unexposed", "Exposed", "Cage"))

ggplot(fitness, aes(x=Treatment, y = Value, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 2, alpha = 1) +
  scale_fill_manual(values = c("cornflowerblue","firebrick2","green","purple")) +
  scale_color_manual(values = c("cornflowerblue","firebrick2","green","purple")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_minimal() +
  labs(x = "Reproductive days") + #"Total offspring "Days to 1st nymph" "Reproductive days"
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_text(size = 12, family=("serif")),
        axis.title.y = element_blank(),
        legend.position = "none", #right #none
        legend.text = element_text(size = 12, family=("serif")),
        legend.title =element_blank())

#methylation genome level multifigure
#load the data
import_list <- list.files(path=wd, pattern = "CpG_evidence.cov$", recursive = TRUE, full.names = TRUE)
column_headers <- c("scaffold","start","end","percent","cReads","tReads")
anc_data <- lapply(import_list, function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- column_headers
  
  return(df)
})
names(anc_data) <- substr(gsub("_mergeCpG.CpG_report.merged_CpG_evidence.cov","",list.files(path=wd, pattern="CpG_evidence.cov$", recursive = TRUE, include.dirs = FALSE, full.names = FALSE), fixed = TRUE), 10,13)[anc]

combo <- bind_rows(anc_data) 
result_anc <- combo %>%
  group_by(scaffold, start, end) %>%
  summarise(totalC = sum(cReads, na.rm = TRUE), totalT = sum(tReads, na.rm = TRUE), .groups = "drop") #drop the grouping so its a regular tibble
result_anc$meth_fraction <- result_anc$totalC / (result_anc$totalC + result_anc$totalT)
result_anc$strand <- "*"
result_anc$cov <- result_anc$totalC + result_anc$totalT

write.table(result_anc, file=paste0(wd,"/ancestral.CpG.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#### stats ####
#start a stats df
statsFile <- data.frame(description=character(), stringsAsFactors = FALSE)
statsFile[nrow(statsFile)+1,] <- c("No of loci analysed")
statsFile[nrow(statsFile)+1,] <- c("No of methylated loci")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated CpGs")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated reads")

#feature calculations
#I annotated ancestral.CpG.txt using bedtools_features.sh
annotated <- read.delim(file=paste0(wd,"/ancestral.CpG.annotated.txt"), sep = "\t", header = FALSE)

#get the gene ID
annotated <- separate(annotated, V17, c("T1","T2","T3"), sep = "\\;", fill = "right", convert=TRUE, remove = FALSE)
annotated$gene <- annotated$T1

#make a final set of data
annotated <- select(annotated, V1, V2, V3, V4, V5, V6, V8, V11, gene)
colnames(annotated) <- c("scaffold","start","end","totalC","totalT","meth_fraction", "cov","feature", "gene")
annotated$id <- paste0(annotated$scaffold,":",annotated$start)
annoUniq <- annotated[!duplicated(annotated[,c('id')]),] #one result per site

#identify methylated CpGs
bt <- function(a, b, p) {binom.test(a, b, p, alternative="greater")$p.value}
annoUniq$pVal <- mapply(bt, annoUniq$totalC, annoUniq$cov, 0.007)
annoUniq$FDR <- p.adjust(annoUniq$pVal, method = "BH", n = length(annoUniq$pVal))
annoUniq$status <- 0
annoUniq$status[annoUniq$FDR < 0.05] <- 1
annoUniq$adjCReads <- annoUniq$totalC 
annoUniq$adjCReads[annoUniq$FDR > 0.05] <- 0
annoUniq$cCount <- 1
annoUniq$feature[annoUniq$feature == "."] <- "intergenic"

#now to the stats themselves
wmpf <- data.frame(description=character(), stringsAsFactors = FALSE)
wmpf[nrow(wmpf)+1,] <- c("exon")
wmpf[nrow(wmpf)+1,] <- c("intergenic")
wmpf[nrow(wmpf)+1,] <- c("intron")
wmpf[nrow(wmpf)+1,] <- c("promoter")

fracF <- data.frame(description=character(), stringsAsFactors = FALSE)
fracF[nrow(fracF)+1,] <- c("exon")
fracF[nrow(fracF)+1,] <- c("intergenic")
fracF[nrow(fracF)+1,] <- c("intron")
fracF[nrow(fracF)+1,] <- c("promoter")

propF <- data.frame(description=character(), stringsAsFactors = FALSE)
propF[nrow(propF)+1,] <- c("exon")
propF[nrow(propF)+1,] <- c("intergenic")
propF[nrow(propF)+1,] <- c("intron")
propF[nrow(propF)+1,] <- c("promoter")

fraction <- annoUniq %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), mCpgCount = sum(status), cpgCount = sum(cCount))

tempStats <- data.frame(count=numeric(), stringsAsFactors = FALSE)
tempStats[nrow(tempStats)+1,] <- format(nrow(annoUniq), scientific = FALSE)
tempStats[nrow(tempStats)+1,] <- format(nrow(subset(annoUniq, status == 1)), scientific = FALSE)
tempStats[nrow(tempStats)+1,] <- format(fraction$mCpgCount/fraction$cpgCount, scientific = FALSE)
tempStats[nrow(tempStats)+1,] <- format(fraction$sumAdjCReads/fraction$sumTotal, scientific = FALSE)

statsFile <- data.frame(description=character(), stringsAsFactors = FALSE)
statsFile[nrow(statsFile)+1,] <- c("No of loci analysed")
statsFile[nrow(statsFile)+1,] <- c("No of methylated loci")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated CpGs")
statsFile[nrow(statsFile)+1,] <- c("Proportion of methylated reads")
statsFile <- cbind(statsFile, tempStats)
colnames(statsFile)[colnames(statsFile) == "count"] <- "ancestral"

featureLevel <- annoUniq %>% group_by(feature) %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = sum(adjCReads)/sum(cov))
featureLevel$proportionMethLoci <- featureLevel$methCLoci / sum(featureLevel$methCLoci)
featureLevel$fractionMethLoci <- featureLevel$methCLoci / featureLevel$totalCLoci

wmpf <- cbind(wmpf, select(featureLevel, weightedMeth))

fracF <- cbind(fracF, select(featureLevel,fractionMethLoci))
colnames(fracF)[colnames(fracF) == "fractionMethLoci"] <- "ancestral"

propF <- cbind(propF, select(featureLevel,proportionMethLoci))
colnames(propF)[colnames(propF) == "proportionMethLoci"] <- "ancestral"

#### doing the same as above but for TEs ####
te_annotated <- read.delim(file=paste0(wd,"/te.CpG.annotated.txt"), sep = "\t", header = FALSE)
te_annotated <- subset(te_annotated, V11 == "transposable_element")
#get the gene ID
te_annotated <- separate(te_annotated, V17, c("T1","T2"), sep = "\\;", fill = "right", convert=TRUE, remove = FALSE)
te_annotated <- separate(te_annotated, T2, c("T3","T4"), sep = "\\/", fill = "right", convert=TRUE, remove = FALSE)
te_annotated$te <- te_annotated$T1
te_annotated$feature <- te_annotated$T3
te_annotated$specific_feature <- te_annotated$T4

#make a final set of data
te_annotated <- select(te_annotated, V1, V2, V3, V4, V5, V6, V8, feature, te, specific_feature)
colnames(te_annotated) <- c("scaffold","start","end","totalC","totalT","meth_fraction", "cov","feature", "te", "specific_feature")
te_annotated$id <- paste0(te_annotated$scaffold,":",te_annotated$start)
te_annoUniq <- te_annotated[!duplicated(te_annotated[,c('id')]),] #one result per site, we don't care which.
te_annoUniq$feature[te_annoUniq$feature == "Retroposon?"] <- "Retroposon"
te_annoUniq$feature[te_annoUniq$feature == "SINE?"] <- "SINE"
te_annoUniq$specific_feature[te_annoUniq$specific_feature == "TcMar-m44?"] <- "TcMar-m44"
te_annoUniq$specific_feature[te_annoUniq$specific_feature == "MITE?"] <- "MITE"

#identify methylated CpGs
bt <- function(a, b, p) {binom.test(a, b, p, alternative="greater")$p.value}
te_annoUniq$pVal <- mapply(bt, te_annoUniq$totalC, te_annoUniq$cov, 0.007)
te_annoUniq$FDR <- p.adjust(te_annoUniq$pVal, method = "BH", n = length(te_annoUniq$pVal))
te_annoUniq$status <- 0
te_annoUniq$status[te_annoUniq$FDR < 0.05] <- 1
te_annoUniq$adjCReads <- te_annoUniq$totalC 
te_annoUniq$adjCReads[te_annoUniq$FDR > 0.05] <- 0
te_annoUniq$cCount <- 1

te_fraction <- te_annoUniq %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), mCpgCount = sum(status), cpgCount = sum(cCount))

te_tempStats <- data.frame(count=numeric(), stringsAsFactors = FALSE)
te_tempStats[nrow(te_tempStats)+1,] <- format(nrow(te_annoUniq), scientific = FALSE)
te_tempStats[nrow(te_tempStats)+1,] <- format(nrow(subset(te_annoUniq, status == 1)), scientific = FALSE)
te_tempStats[nrow(te_tempStats)+1,] <- format(te_fraction$mCpgCount/te_fraction$cpgCount, scientific = FALSE)
te_tempStats[nrow(te_tempStats)+1,] <- format(te_fraction$sumAdjCReads/te_fraction$sumTotal, scientific = FALSE)

te_statsFile <- data.frame(description=character(), stringsAsFactors = FALSE)
te_statsFile[nrow(te_statsFile)+1,] <- c("No of loci analysed")
te_statsFile[nrow(te_statsFile)+1,] <- c("No of methylated loci")
te_statsFile[nrow(te_statsFile)+1,] <- c("Proportion of methylated CpGs")
te_statsFile[nrow(te_statsFile)+1,] <- c("Proportion of methylated reads")
te_statsFile <- cbind(te_statsFile, te_tempStats)
colnames(te_statsFile)[colnames(statsFile) == "count"] <- "ancestral"

te_featureLevel <- te_annoUniq %>% group_by(feature) %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = sum(adjCReads)/sum(cov))
te_featureLevel$proportionMethLoci <- te_featureLevel$methCLoci / sum(te_featureLevel$methCLoci)
te_featureLevel$fractionMethLoci <- te_featureLevel$methCLoci / te_featureLevel$totalCLoci

te_spfeatureLevel <- te_annoUniq %>% group_by(specific_feature, feature) %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = sum(adjCReads)/sum(cov))
te_spfeatureLevel$proportionMethLoci <- te_spfeatureLevel$methCLoci / sum(te_spfeatureLevel$methCLoci)
te_spfeatureLevel$fractionMethLoci <- te_spfeatureLevel$methCLoci / te_spfeatureLevel$totalCLoci

#####typed the necessary values generated above into a text file for the bar chart ####
barChart <- read.delim(file=paste0(wd,"/full_bar_chart.txt"), header = TRUE, sep = "\t")

ggplot(barChart, aes(x=Feature, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(labels = c("01_genome" = "Genome", "02_exons" = "Exon", "03_introns" = "Intron", "04_promoters" = "Promoters", "05_intergenic" = "Intergenic", "06_DNA"="DNA", "07_LINE"="LINE","08_LTR"="LTR","09_MITE"="MITE", "10_RC"="RC","11_Retroposon"="Retroposon", "12_SINE"="SINE")) +
  labs(x = "Feature", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust = 1, family=("serif")),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, family=("serif")),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family=("serif")),
        legend.position = "top")

#wmpg and expression analysis
annoUniq_wmpg <- subset(annoUniq, feature != "promoter" & feature != "intergenic")
wmpg <- select(annoUniq_wmpg %>% count(gene), gene=gene)
geneCat <- wmpg
catSum <- data.frame(category=character(), stringsAsFactors = FALSE)
catSum[nrow(catSum)+1,] <- c("High")
catSum[nrow(catSum)+1,] <- c("Low")
catSum[nrow(catSum)+1,] <- c("Medium")
catSum[nrow(catSum)+1,] <- c("Unmethylated")

geneLevel <- annoUniq_wmpg %>% group_by(gene) %>% summarise(sumCReads = sum(totalC), sumAdjCReads = sum(adjCReads), sumTReads = sum(totalT), sumTotal = sum(cov), methCLoci = sum(status), totalCLoci = sum(cCount), weightedMeth = format(sum(adjCReads)/sum(cov), scientific = FALSE))
geneLevel$fractionMethLoci <- geneLevel$methCLoci / geneLevel$totalCLoci
wmpg <- left_join(wmpg, select(geneLevel, gene, weightedMeth), by = "gene")
colnames(wmpg)[colnames(wmpg) == "weightedMeth"] <- "ancestral"
geneLevel$cat <- "U"
geneLevel$cat[geneLevel$weightedMeth >= 0.05 & geneLevel$weightedMeth < 0.3] <- "L"
geneLevel$cat[geneLevel$weightedMeth >= 0.3 & geneLevel$weightedMeth < 0.7] <- "M"
geneLevel$cat[geneLevel$weightedMeth >= 0.7] <- "H"
sampleCat <- geneLevel %>% count(cat)
catSum <- cbind(catSum, select(sampleCat, n))
colnames(catSum)[colnames(catSum) == "n"] <- "sus"
geneCat <- left_join(geneCat, select(geneLevel, gene, cat), by = "gene")
colnames(geneCat)[colnames(geneCat) == "cat"] <- "ancestral"
write.table(geneLevel, file=paste0(wd,"/geneLevel.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(featureLevel, file=paste0(wd,"/featureLevel.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(catSum, file=paste0(wd,"/catSum.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#violin plots
#uses wmpg and FPKM 
wmpgNA <- na.omit(wmpg)
FPKM <- read.delim(file=paste0(wd,"/FPKM.txt"), sep="\t", header = TRUE)
FPKMNa <- na.omit(wmpg)
both <- merge(wmpgNA, FPKM, by = "gene")
colnames(both) <- c("gene","meth","exp")
#both <- subset(both, both$exp != 0)
both <- subset(both, both$exp >= 1)
both$exp <- log10(both$exp)
both$exp[both$exp == -Inf] <- 0
both$bin[both$meth <= 0.05] <- "0_5"
both$bin[both$meth > 0.05 & both$meth <= 0.3] <- "5_30"
both$bin[both$meth > 0.3 & both$meth <= 0.7] <- "30_70"
both$bin[both$meth > 0.7] <- "70_100"
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
level_order <- c("0_5","5_30","30_70","70_100")
ggplot(both, aes(x=factor(bin, levels = level_order), y=exp, fill=bin)) + 
  geom_violin(trim=FALSE) + stat_summary(fun.data=data_summary) + theme_minimal() +labs(x="Methylation bin", y = "log10(FPKM)", fill = "Mean weighted methylation per gene") +
  scale_x_discrete(labels = c("0_5" = "<= 5%", "5_30" = "6 - 30%", "30_70" = "31 - 70%", "70_100" = "71% +")) +
  labs(x = "Weighted Methylation") +
  theme(axis.text.x = element_text(size = 12, family=("serif")),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_text(size = 12, family=("serif")),
        axis.title.y = element_text(size = 12, family=("serif")),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family=("serif")),
        legend.position = "none")
kruskal.test(exp ~ bin, data = both)
library(FSA)
dunnTest(exp ~ bin, data=both, method = "bonferroni")

#metagene plot
#with thanks to https://github.com/groverj3/genomics_visualizations/blob/master/metaplotteR.r 
#bedgraphs for ancestral sample created using CHH_CHG.R, and analysed with deeptools in metagene_all_contexts.sh
#.tab output replotted here.
read_deeptools_table <- function(file) {
  
  n <- max(count.fields(file, sep = '\t'), na.rm = TRUE)
  x <- readLines(file)
  
  .splitvar <- function(x, sep, n) {
    var <- unlist(strsplit(x, split = sep))
    length(var) <- n
    return(var)
  }
  
  x <- do.call(cbind, lapply(x, .splitvar, sep = '\t', n = n))
  x <- apply(x, 1, paste, collapse = '\t')
  plot_table <- na.omit(read.csv(text = x, sep = '\t')[-1,])  # Remove first row with "gene" label
  
  return(plot_table)
}
tab <- read_deeptools_table(paste0(wd,"/all_contexts.tab"))
long_table <- gather(tab, 'sample', 'score', -bin.labels, -bins)
ggplot(long_table, aes(x = bins, y = as.numeric(score), color = sample)) +
  geom_line() +
  labs(y="Methylation") +
  scale_x_continuous(breaks = long_table$bins, labels = long_table$bin.labels) +
  theme(axis.text.x = element_text(size = 12, family=("serif")),
        axis.text.y = element_text(size = 12, family=("serif")),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, family=("serif")),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 12, family=("serif")),
        legend.position = "right")

tab <- read.delim(paste0(wd,"/all_contexts.tab"), header=F, skip=2)

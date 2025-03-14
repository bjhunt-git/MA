#script for creating separate bedGraphs for each context using the CX info 
#to distinguish

library(tidyverse) #general wrangling
library(data.table) #also wrangling and processing the input data

import_path <- "projects/MA/data/cpg"

#NEED TO SKIP THE HEADER
bedgraph <- read.table(file=paste0(import_path,"/D59_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE,skip = 1)
cx_report <- read.table(file=paste0(import_path,"/D59_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#make id - the site in CX is the end of the bedgraph
bedgraph$id <- paste0(bedgraph$V1,":",bedgraph$V3)
cx_report$id <- paste0(cx_report$V1,":",cx_report$V2)
#filter CX for coverage (10 reads) and calculate methylation
cx_report$cov <- cx_report$V4 + cx_report$V5
cx_report_filter <- subset(cx_report, cx_report$cov >= 10)
#cx_report_filter$meth <- cx_report_filter$V4 / cx_report_filter$cov # not needed
#split CX into the three contexts
cx_cpg <- subset(cx_report_filter, V6 == "CG")
cx_chg <- subset(cx_report_filter, V6 == "CHG")
cx_chh <- subset(cx_report_filter, V6 == "CHH")
#subset bedgraph where id = for each context
cg_bedgraph <- merge(bedgraph, cx_cpg, by = "id")
cg_bedgraph <- select(cg_bedgraph, V1.x, V2.x, V3.x, V4.x)
cg_bedgraph <- arrange(cg_bedgraph,V1.x,V2.x)
chg_bedgraph <- merge(bedgraph, cx_chg, by = "id")
chg_bedgraph <- select(chg_bedgraph, V1.x, V2.x, V3.x, V4.x)
chg_bedgraph <- arrange(chg_bedgraph,V1.x,V2.x)
chh_bedgraph <- merge(bedgraph, cx_chh, by = "id")
chh_bedgraph <- select(chh_bedgraph, V1.x, V2.x, V3.x, V4.x)
chh_bedgraph <- arrange(chh_bedgraph,V1.x,V2.x)
#save bedgraphs
write.table(cg_bedgraph, file=paste0(import_path,"/cg.bedGraph"), col.names = FALSE, sep="\t", quote = FALSE, row.names = FALSE)
write.table(chg_bedgraph, file=paste0(import_path,"/chg.bedGraph"), col.names = FALSE, sep="\t", quote = FALSE, row.names = FALSE)
write.table(chh_bedgraph, file=paste0(import_path,"/chh.bedGraph"), col.names = FALSE, sep="\t", quote = FALSE, row.names = FALSE)

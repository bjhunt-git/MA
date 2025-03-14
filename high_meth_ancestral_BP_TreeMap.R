# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes
# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006281","DNA repair",2.895534847137308,4.6439741428068775,0.6055926239456005,0,"DNA repair"),
                     c("GO:0006364","rRNA processing",1.5444135826112384,1.6210292990973005,0.4382295832574098,0.69501634,"DNA repair"),
                     c("GO:0006397","mRNA processing",1.242261085587905,4.041436116778033,0.519821773752522,0.43597524,"DNA repair"),
                     c("GO:0006468","protein phosphorylation",1.0150465485898383,1.9711627217697396,0.6404029328825717,0.66477721,"DNA repair"),
                     c("GO:0006486","protein glycosylation",0.7526798873357411,2.8135318847986124,0.6259644837423973,0.20449964,"DNA repair"),
                     c("GO:0006497","protein lipidation",0.23154015599044747,1.31738280157417,0.6572758344734632,0.57511239,"DNA repair"),
                     c("GO:0006511","ubiquitin-dependent protein catabolic process",1.238068562564521,2.57984481908204,0.681707889361453,0.35409187,"DNA repair"),
                     c("GO:0008380","RNA splicing",0.9742871404547958,1.785247893110021,0.5569234792694018,0.65951647,"DNA repair"),
                     c("GO:0043039","tRNA aminoacylation",0.8820269865386736,1.4767973397455239,0.6049062587798882,0.50524739,"DNA repair"),
                     c("GO:0006661","phosphatidylinositol biosynthetic process",0.44089272100117355,2.2307264590699574,0.8010163225414003,0.09583789,"phosphatidylinositol biosynthetic process"),
                     c("GO:0006888","endoplasmic reticulum to Golgi vesicle-mediated transport",0.41624285459502003,2.409271377303531,0.9442554925855289,0.01167906,"endoplasmic reticulum to Golgi vesicle-mediated transport"),
                     c("GO:0071705","nitrogen compound transport",6.023736236521199,1.3465756600906729,0.9614039794992553,0.33129429,"endoplasmic reticulum to Golgi vesicle-mediated transport"),
                     c("GO:0017148","negative regulation of translation",0.21779183534115995,1.514042231198597,0.9129020085244949,-0,"negative regulation of translation"),
                     c("GO:0022402","cell cycle process",2.0165345615232173,1.5018351592170025,0.9794307657015331,0.01415401,"cell cycle process"),
                     c("GO:0065003","protein-containing complex assembly",2.2933051643109175,1.6052123437617338,0.8485352925590329,0.0144028,"protein-containing complex assembly"),
                     c("GO:0006996","organelle organization",6.781836087483291,3.779891911959945,0.8616385308053325,0.56756569,"protein-containing complex assembly"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_070524.pdf", width=10, height=6 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
#output this as a pdf and the annoying traces of the categorical stuff vanishes.
treemap(
  stuff,
  index = c("representative","description"),
  #index = c("description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  #title = "Highly methylated genes (Biological Process)",
  fontsize.title=c(0),
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = 0,   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "right",
  fontsize.legend=c(10),
  fontsize.labels=c(0,12),
  title.legend = "Representative",
  palette="Set3"
)

dev.off()



library(AlphaBeta)
library(dplyr)

#### unexposed lines
setwd("projects/MA/analysis/wgbs/alphabeta/control")

sampleFile <- "control_endpoint_nodelist.fn"
edgesFile <- "control_endpoint_edgelist.fn"
out.dir <- getwd()

output <- buildPedigree(nodelist = sampleFile, edgelist = edgesFile, cytosine = "CG",posteriorMaxFilter = 0.99)
head(output$Pdata)

plotPedigree(nodelist = sampleFile, edgelist = edgesFile, sampling.design = "progenitor.intermediate",
             output.dir = out.dir, plot.width = 5, plot.height = 5, aspect.ratio = 1, vertex.size = 6,
             vertex.label = TRUE, out.pdf = "control_pedigree")

pedigree <- output$Pdata
dt <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
pdf(paste0(out.dir,"/control_divergence.pdf"))
plot(dt, pedigree[, "D.value"], ylab = "Divergence value", xlab = expression(paste(Delta," t")))
dev.off()

p0uu_in <- output$tmpp0
p0uu_in
pedigree <- output$Pdata
output.data.dir <- out.dir
outputABneutral <- ABneutral(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in,
  eqp.weight = 1,
  Nstarts = 50,
  out.dir = output.data.dir,
  out.name = "control_ABneutral_CG_global_estimates"
)

head(outputABneutral$estimates)
#alpha beta
#3.244057e-05 0.0002746236
ABfile <- paste0(out.dir,"/control_ABneutral_CG_global_estimates.Rdata")
ABplot(pedigree.names = ABfile, output.dir = out.dir, out.name = "control_ABneutral", plot.height = 8,
       plot.width = 11)

outputABnull <- ABnull(pedigree.data = pedigree, out.dir = output.data.dir, out.name = "control_ABnull_CG_global_estimates")
outputABnull
nFile <- paste0(out.dir,"/control_ABnull_CG_global_estimates.Rdata")

outputABneutralSOMA <- ABneutralSOMA(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in, 
  eqp.weight = 0.001,
  Nstarts = 50, 
  out.dir = getwd(), 
  out.name = "control_ABneutralSOMA_CG_global_estimates")

head(outputABneutralSOMA$estimates)
#alpha        beta
#2.858445e-05 0.0004928439
somaFile <- paste0(out.dir,"/control_ABneutralSOMA_CG_global_estimates.Rdata")
ABplot(pedigree.names = ABfile, output.dir = out.dir, out.name = "control_ABneutralSOMA", plot.height = 8,
       plot.width = 11)

outFtest <- FtestRSS(pedigree.select = ABfile, pedigree.null = nFile)
outFtest$Ftest
#$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#2.143442e-02 2.291327e-02 3.200000e+02 3.240000e+02 5.519528e+00 2.618484e-04

outFtest <- FtestRSS(pedigree.select = somaFile, pedigree.null = nFile)
outFtest$Ftest
#$Ftest
#       RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#2.135540e-02 2.291327e-02 3.200000e+02 3.240000e+02 5.835951e+00 1.520821e-04

outFtest <- FtestRSS(pedigree.select = somaFile, pedigree.null = ABfile)
outFtest$Ftest

#### check for selection 
outputABselectMMSOMA <- ABselectMMSOMA(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in,
  eqp.weight = 0.001,
  Nstarts = 50,
  out.dir = output.data.dir,
  out.name = "control_ABselectMMSOMA_CG_global_estimates"
)
head(outputABselectMMSOMA$estimates)
ABselectMMSOMAfile <- paste0(out.dir,"/control_ABselectMMSOMA_CG_global_estimates.Rdata")
ABplot(pedigree.names = ABselectMMSOMAfile, output.dir = out.dir, out.name = "control_ABselectMMSOMA", plot.height = 8,
       plot.width = 11)

outFtest <- FtestRSS(pedigree.select = ABselectMMSOMAfile, pedigree.null = nFile)
outFtest$Ftest
#       RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#2.135187e-02 2.291327e-02 3.200000e+02 3.240000e+02 5.850139e+00 1.484203e-04
#therefore MM Soma better than null

outFtest <- FtestRSS(pedigree.select = ABselectMMSOMAfile, pedigree.null = ABfile)
outFtest$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#0.02135187   0.02143442 320.00000000 321.00000000   1.23709315   0.26686689
#therefore ABNeutral better than MM SOMA

outFtest <- FtestRSS(pedigree.select = ABselectMMSOMAfile, pedigree.null = somaFile)
outFtest$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#0.02135187   0.02135540 320.00000000 321.00000000   0.05289431   0.81824833
#Therefore SOMA better than MM SOMA

outputABselectUUSOMA <- ABselectUUSOMA(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in,
  eqp.weight = 0.001,
  Nstarts = 50,
  out.dir = output.data.dir,
  out.name = "control_ABselectUUSOMA_CG_global_estimates"
)
head(outputABselectUUSOMA$estimates)
ABselectUUSOMAfile <- paste0(out.dir,"/control_ABselectUUSOMA_CG_global_estimates.Rdata")
ABplot(pedigree.names = ABselectUUSOMAfile, output.dir = out.dir, out.name = "control_ABselectUUSOMA", plot.height = 8,
       plot.width = 11)

outFtest <- FtestRSS(pedigree.select = ABselectUUSOMAfile, pedigree.null = nFile)
outFtest$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#2.132594e-02 2.291327e-02 3.200000e+02 3.240000e+02 5.954529e+00 1.240532e-04
#Therefore UU SOMA better than neutral

outFtest <- FtestRSS(pedigree.select = ABselectUUSOMAfile, pedigree.null = ABfile)
outFtest$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#0.02132594   0.02143442 320.00000000 321.00000000   1.62770561   0.20294599
#Therefore ABNeutral better than UU SOMA

outFtest <- FtestRSS(pedigree.select = ABselectUUSOMAfile, pedigree.null = somaFile)
outFtest$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#0.02132594   0.02135540 320.00000000 321.00000000   0.44206682   0.50660544
#Therefore SOMA better than UU SOMA
#bootstrap the preferred model
inputModel <- "control_ABneutralSOMA_CG_global_estimates.Rdata"
Boutput <- BOOTmodel(pedigree.data = inputModel, Nboot = 50, out.dir = getwd(), out.name = "control_ABneutralSOMABOOT")
Boutput
#$standard.errors
#SE          2.5%        97.5%
#alpha      6.380283e-06  1.470306e-05 3.750605e-05
#beta       1.110963e-04  2.557320e-04 6.534922e-04
#beta/alpha 1.399004e-02  1.737938e+01 1.743094e+01
#weight     2.031482e-01 -2.647788e-01 5.611673e-01
#intercept  1.225719e-03  5.706570e-02 6.140502e-02
#PrMMinf    4.487292e-06  2.943429e-03 2.960319e-03
#PrUMinf    7.361759e-05  1.026199e-01 1.028970e-01
#PrUUinf    7.810488e-05  8.941427e-01 8.944367e-01

#### exposed lines
setwd("projects/MA/analysis/wgbs/alphabeta/exposed")
sampleFile <- "exposed_endpoint_nodelist.fn"
edgesFile <- "exposed_endpoint_edgelist.fn"
out.dir <- getwd()

#build the pedigree and plot
output <- buildPedigree(nodelist = sampleFile, edgelist = edgesFile, cytosine = "CG",posteriorMaxFilter = 0.99)
head(output$Pdata)

plotPedigree(nodelist = sampleFile, edgelist = edgesFile, sampling.design = "progenitor.intermediate",
             output.dir = out.dir, plot.width = 5, plot.height = 5, aspect.ratio = 1, vertex.size = 6,
             vertex.label = TRUE, out.pdf = "full_pedigree")

#plot the divergence values
pedigree <- output$Pdata
dt <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
pdf(paste0(out.dir,"/exposed_divergence.pdf"))
plot(dt, pedigree[, "D.value"], ylab = "Divergence value", xlab = expression(paste(Delta," t")))
dev.off()

#calculate alpha/beta under a model of neutral epimutation accumulation (no selection)
p0uu_in <- output$tmpp0
p0uu_in
pedigree <- output$Pdata
output.data.dir <- out.dir
outputABNeutral <- ABneutral(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in,
  eqp.weight = 1,
  Nstarts = 50,
  out.dir = output.data.dir,
  out.name = "exposed_ABneutral_CG_global_estimates"
)

outputABNeutral$estimates

#load the data into a file to plot 
ABfile <- paste0(out.dir,"/exposed_ABneutral_CG_global_estimates.Rdata")
ABplot(pedigree.names = ABfile, output.dir = out.dir, out.name = "exposed_ABneutral", plot.height = 8,
       plot.width = 11)

#calculate alpha/beta under a model of no epimutation accumulation
outputABnull <- ABnull(pedigree.data = pedigree, out.dir = output.data.dir, out.name = "exposed_ABnull_CG_global_estimates")
outputABnull$estimates

nFile <- paste0(out.dir,"/exposed_ABnull_CG_global_estimates.Rdata")

#compare the neutral and null model - significance depicts the superior model for explaining the data.
outFtest <- FtestRSS(pedigree.select = ABfile, pedigree.null = nFile)
outFtest
#$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#7.157728e-03 7.384080e-03 3.200000e+02 3.240000e+02 2.529868e+00 4.053281e-02

#calculate alpha/beta using a model designed for asexual, clonal and somatic systems (no selection)
outputABneutralSOMA <- ABneutralSOMA(
  pedigree.data = pedigree,
  p0uu = p0uu_in,
  eqp = p0uu_in, 
  eqp.weight = 0.001,
  Nstarts = 50, 
  out.dir = getwd(), 
  out.name = "exposed_ABneutralSOMA_CG_global_estimates")
outputABneutralSOMA$estimates
#alpha        beta
#1.542383e-05 0.0003218671

#compare the neutral SOMA and null model - significance depicts the superior model for explaining the data.
somaFile <- "exposed_ABneutralSOMA_CG_global_estimates.Rdata"
outFtest <- FtestRSS(pedigree.select = somaFile, pedigree.null = nFile)
outFtest
#$Ftest
#RSS_F        RSS_R         df_F         df_R       Fvalue       pvalue 
#7.044176e-03 7.384080e-03 3.200000e+02 3.240000e+02 3.860252e+00 4.445119e-03

#so as per the reporting in the alphabeta paper:
#context  annotation  alpha beta  beta/alpha  FM  RM  F-value df FM df RM pvalue
#CG global  1.542383e-05  3.218671e-04  
#bootstrap the preferred model
inputModel <- "exposed_ABneutralSOMA_CG_global_estimates.Rdata"
Boutput <- BOOTmodel(pedigree.data = inputModel, Nboot = 50, out.dir = getwd(), out.name = "ABneutralSOMABOOT")
Boutput
#$standard.errors
#                     SE          2.5%        97.5%
#alpha      5.146483e-06  3.541325e-06 2.143561e-05
#beta       1.077074e-04  7.411900e-05 4.486631e-04
#beta/alpha 4.550535e-03  2.092487e+01 2.093876e+01
#weight     1.638980e+00 -3.328822e+00 1.201090e+00
#intercept  6.516439e-04  4.484448e-02 4.690821e-02
#PrMMinf    8.622535e-07  2.077381e-03 2.080301e-03
#PrUMinf    1.718967e-05  8.700182e-02 8.706002e-02
#PrUUinf    1.805192e-05  9.108597e-01 9.109208e-01

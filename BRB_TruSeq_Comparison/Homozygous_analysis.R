##################################################################
### Analysis of homozygous line in round robin BRB-seq samples ###
##################################################################

### Script informations ####
# cf title

### Miscellaneous ####
## Paths ####
data.path <- "./Data"
figures.path <- "./Figures"
processed.path <- "./Processed"

## data
rna.cpm <- "./Data/cpm.table.76.lines.txt"
rna.DE <- "./Data/tt.condition.txt"

# First time running script only: create directories in the current one to store modified data and figures
if(!dir.exists(figures.path)) { dir.create(figures.path) }
if(!dir.exists(processed.path)) { dir.create(processed.path) }

## libraries ####
library(reshape2)
library(ggplot2)
library(cowplot)
library(gplots)
# library(FactoMineR)
# library(tidyr) #fct separate
# library(factoextra) #fviz fct
library(edgeR)
library(statmod)
library(viridis)
library(MASS)

## colPanels ####

## functions ####
plot.clustering = function(data.input, filename, type) 
{
  png(filename, width=1000, height=1000)
  par(mfrow = c(2,2), mar = c(7,4,1,1))
  
  # Distance Matrix
  cor.dist.matrix = as.dist((1 - cor(data.input))/2)
  
  # HC
  data.hc = hclust(cor.dist.matrix, method="ward.D2")
  plot(data.hc, main=paste0("HC on ",type,", WARD"))
  data.hc = hclust(cor.dist.matrix)
  plot(data.hc, main=paste0("HC on ",type,", COMPLETE"))
  
  # MDS
  #fit <- cmdscale(cor.dist.matrix, eig=TRUE, k=2) # k is the number of dim
  #plot(fit$points[,1:2], xlab="Coordinate 1", ylab="Coordinate 2", main=paste0("MDS on ",type),	type="n")
  #text(fit$points[,1:2], labels = colnames(data.input))
  
  # PCA
  data.pca = princomp(data.input)
  plot(data.pca$loadings[,1:2], main=paste0("PCA on ",type), type="n") # 2D PCA
  text(data.pca$loadings[,1:2], labels = colnames(data.input))
  dev.off()
}

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 20)
}

## vector with informations ####
# number of lines in experimental design
line.nb <- c("2","7","10","14")
treat.vect <- c("UC", "Pe")
replica.vect <- c("1","2")

### END ####

### Get all counts in one file and create design matrix ####
samples <- apply(expand.grid(line.nb, treat.vect, replica.vect), 1, paste, collapse="_")

for (sam in samples) {
    f.path <- file.path(data.path, sam, "output.counts.txt")
    #print the cross name to follow loop execution
    print(sam)
    print(f.path)
    #open files
    tmp.cnts <- read.delim(f.path, header=F)
    # add columns name
    colnames(tmp.cnts) <- c("FB.gene", sam)
    if (!exists("counts")){
      counts <- tmp.cnts
    }
    # if the merged dataset does exist, append to it
    else {
      counts <- merge(counts, tmp.cnts, by="FB.gene")
   }
}
rm(tmp.cnts)

row.names(counts) <- counts$FB.gene
counts$FB.gene <- NULL
saveRDS(counts, file.path(processed.path, "count_matrix.rds"))
rm(counts)

model.mat <- expand.grid(line.nb, treat.vect, replica.vect)
colnames(model.mat) <- c("Lines", "Treatment", "Replica")
rownames(model.mat) <- paste(model.mat$Lines, model.mat$Treatment, model.mat$Replica, sep="_")
model.mat$Lines.Treatment <- paste(model.mat$Lines, model.mat$Treatment, sep=".")
saveRDS(model.mat, file.path(processed.path, "model_matrix.rds"))
### END ####

### Analysis ####

## DE genes ####
#load data
counts <- readRDS(file.path(processed.path, "count_matrix.rds"))
model.mat <- readRDS(file.path(processed.path, "model_matrix.rds"))
RNA.data <- read.table(rna.cpm, header=T, row.names = 1)

# Filtering low expressed genes
data.dge = DGEList(counts = counts)

data.dge = data.dge[rowSums(log2(cpm(data.dge)+1)) > 0.5, ]
data.filtered = DGEList(counts = counts[rownames(data.dge$counts), ])
data.filtered = calcNormFactors(data.filtered)

# design matrix
groups <- model.mat$Lines.Treatment
names(groups) <- row.names(model.mat)
design <- model.matrix(~model.mat$Treatment)
design

# Voom normalization of the count matrix
data.voom = voom(counts = data.filtered, design, normalize.method="quantile", plot=F) #round 1
dc = duplicateCorrelation(data.voom, design = design, block = model.mat$Lines)
dt.noBatch = removeBatchEffect(data.voom$E, batch = model.mat$Replica, design = design)
plot.clustering(data.voom$E, file.path(figures.path, "voom_clustering.png"), "VOOM")

saveRDS(data.voom$E, file.path(processed.path, "Normalized_Cnt_Matrix.rds"))

fit= lmFit(data.voom, design = design, correlation = dc$consensus.correlation, block = model.mat$Lines)
#fit= lmFit(data.voom, design = design)

fit2 = eBayes(fit)

## Trueseq-BRB correlation on DE genes ####
BRB.DE <- topTable(fit2, coef = 2, n = Inf)
RNA.DE <- read.table(rna.DE, header=T, row.names = 1)

DE.merge <- merge(BRB.DE, RNA.DE, by=0, all = T)
DE.cor <- cor.test(DE.merge$logFC.x, DE.merge$logFC.y, method = "pearson", use="pairwise.complete.obs")
pdf(file.path(figures.path, "logFC_correlation.pdf"))
plot(DE.merge$logFC.x, DE.merge$logFC.y, main="logFC correlation", xlab="BRB", ylab="Truseq", 
     sub = paste("pearson:", round(DE.cor$estimate,3), sep=" "))
dev.off()

ggplot(data = DE.merge, aes(x=logFC.x, y=logFC.y)) + geom_point() + geom_density_2d() + 
  labs(title="LogFC correlation", subtitle=paste("pearson:", round(DE.cor$estimate,3), sep=" ")) +
  xlab("BRB") + ylab("TruSeq")
## END ####

### Comparisons on counts ####

RNA.data <- read.table(rna.cpm, header=T, row.names = 1)
BRB.data <- readRDS(file.path(processed.path, "Normalized_Cnt_Matrix.rds"))
## Heatmaps ####
# samples correlation
cor.mat <- round(cor(BRB.data, method = "spearman"),2)
heatmap.2(cor.mat,
          Rowv = T, Colv = T, #no dendogram and rearrangement
          trace = "none",
          denscol = "black",
          cexRow = 0.8, cexCol = 0.8,
          offsetRow = 0, offsetCol = 0, srtCol = 45,
          key.xlab = "Pearson correlation")

## END Heatmaps ####

## Barplot of nb of reads ####
jpeg(file.path(figures.path, "cnt_per_samples.jpg"))
barplot(colSums(counts), las=2)
dev.off()
## END ####

## Trueseq-BRB correlation on cpm ####

brb_tru <- merge(BRB.data, RNA.data, by=0, all=F)

# subset dataset
DGRP.386.UC <- c("2_UC_1","2_UC_2", "DGRP.386_Naive")
DGRP.386.Pe <- c("2_Pe_1","2_Pe_2", "DGRP.386_Treated")
DGRP.313.UC <- c("7_UC_1","7_UC_2", "DGRP.313_Naive")
DGRP.313.Pe <- c("7_Pe_1","7_Pe_2", "DGRP.313_Treated")
DGRP.486.UC <- c("10_UC_1","10_UC_2", "DGRP.486_Naive")
DGRP.486.Pe <- c("10_Pe_1","10_Pe_2", "DGRP.486_Treated")
DGRP.890.UC <- c("14_UC_1","14_UC_2", "DGRP.890_Naive")
DGRP.890.Pe <- c("14_Pe_1","14_Pe_2", "DGRP.890_Treated")
general.colnames <- c("BRB 1", "BRB 2", "TruSeq")

## Create the plots ####
df.plot <- brb_tru[, DGRP.386.UC]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr1.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 386 UC")
dev.off()

df.plot <- brb_tru[, DGRP.386.Pe]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr2.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 386 Pe")
dev.off()

df.plot <- brb_tru[, DGRP.313.UC]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr3.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 313 UC")
dev.off()

df.plot <- brb_tru[, DGRP.313.Pe]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr4.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 313 Pe")
dev.off()

df.plot <- brb_tru[, DGRP.486.UC]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr5.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 486 UC")
dev.off()

df.plot <- brb_tru[, DGRP.486.Pe]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr6.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 486 Pe")
dev.off()

df.plot <- brb_tru[, DGRP.890.UC]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr7.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 890 UC")
dev.off()

df.plot <- brb_tru[, DGRP.890.Pe]
colnames(df.plot) <- general.colnames
pdf(file.path(figures.path, "Corr8.pdf"))
pairs(df.plot, 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="DGRP 890 Pe")
dev.off()

### END ####
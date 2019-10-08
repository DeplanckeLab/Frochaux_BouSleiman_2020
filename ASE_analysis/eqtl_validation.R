#########################
### eQTL validation 1 ###
#########################

# Description
# General stats and analysis on the data
### useful stuff at the beginnings ####

## paths ####
fig.path <- file.path("./Figures")
processed.path <- file.path("./Processed")

# First script only: create two directories in the current one to store modified data and figures
if(!dir.exists(fig.path)) { dir.create(fig.path) }
if(!dir.exists(processed.path)) { dir.create(processed.path) }

## libraries ####
library(tictoc)
library(reshape2) # melt function
library(ggplot2)
library(tidyr)
library(cowplot)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(gplots)
library(dichromat)

## functions ####
# sub triangle
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# sup triangle
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#reordering correlation matrix
reorder_cormat <- function(cormat){
  # Utiliser la corr?lation entre les variables
  # comme mesure de distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

## vectors with informations ####
cross.vect <- c("c1x2", "c2x3", "c3x4", "c4x5", "c5x6", "c6x7", "c7x8", "c8x9",
                "c9x10", "c10x11", "c11x12", "c12x13", "c13x14", "c14x15",
                "c17x18", "c18x19", "c19x20", "c20x1")
treat.vect.pe <- c("Pe_1", "Pe_2")
treat.vect.uc <- c("UC_1", "UC_2")
#type.vect <- c("alt", "ref")
type.vect <- c("alt", "ref", "total")
all.treat <- c("Pe_1", "Pe_2", "UC_1", "UC_2")


## Color panels ####
colPanel.infection <- c("#c0c0c0", "#cc6633")
### END Beginning ####

#############################
###--- Data aggregation --###
# Need to be done once only #
#############################

# aggregate all the cnt from  all the count files into one file
# create count matrix

# output: file with 4 columns, all data appended
for (cr in cross.vect) {
  cross.filename <- gsub("c", "",cr)
  for (t in all.treat) {
    f.pat <- paste(cross.filename, t, sep="_")
    #print the cross name to follow loop execution
    print(f.pat)
    #open files
    tmp.cnts <- read.delim(file.path(Gene_cnts.path, list.files(Gene_cnts.path, pattern = f.pat)), header=T)
    # add sample columns
    tmp.cnts$sample <- f.pat
    # sum cnts for both lines
    colnames(tmp.cnts)[4:5] <- c("l1", "l2")
    tmp.cnts$assigned <- tmp.cnts$l1 + tmp.cnts$l2
    if (!exists("counts")){
      counts <- tmp.cnts[,c(1,8,9,10)]
    }
    # if the merged dataset does exist, append to it
    else {
      counts<- rbind(counts, tmp.cnts[,c(1,8,9,10)])
    }
  }
}
rm(tmp.cnts)

# remove sample c19x20_Pe_1 and c20x1_Pe_2
counts <- counts[ ! counts$sample %in% c("19x20_Pe_1"), ]
#counts <- counts[ ! counts$sample %in% c("19x20_Pe_1"), ]
if (nrow(counts) / 71 != 14869) {print("error in removing samples")}

saveRDS(counts, file = "./Processed/all_counts.rds")

# make count matrix
count.matrix <- unstack(counts, form=covering_reads~sample)
row.names(count.matrix) <- unique(counts$gene)
colnames(count.matrix) <- gsub(pattern = "X", replacement = "c", x = colnames(count.matrix))

saveRDS(count.matrix, "./Processed/Count_matrix_all.rds")

total.cnt <- as.data.frame(colSums(count.matrix))
total.cnt$sample <- row.names(total.cnt)
total.cnt <- total.cnt[order(total.cnt$`colSums(count.matrix)`),]
colnames(total.cnt)[1] <- "cnt"

pdf("./Figures/sample_cnts.pdf", width=11, height=7)
barplot(total.cnt$cnt, names.arg = total.cnt$sample, las=2)
dev.off()

rm(counts)
rm(count.matrix)


### END Data aggregation ####


############################
### stats on read counts ###
############################

cnt.mat <- readRDS("./Processed/Count_matrix_all.rds")

## correlation matrix ####
cor.mat <- round(cor(cnt.mat),2)

### Plot 1: heatmap of correlation matrix between samples ####

heatmap(cor.mat)

# heatmap.2
column.treat <- sapply(strsplit(x = colnames(cor.mat), "_"), "[[", 2)
column.treat <- gsub("UC", "#c0c0c0", column.treat)
column.treat <- gsub("Pe", "#cc6633", column.treat)

row.treat <- sapply(strsplit(x = row.names(cor.mat), "_"), "[[", 2)
row.treat <- gsub("UC", "#c0c0c0", row.treat)
row.treat <- gsub("Pe", "#cc6633", row.treat)

color.vect <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
column.label <- sapply(strsplit(x = colnames(cor.mat), "_"), "[[", 1)
color.v2 <- color.vect[as.factor(column.label)]

row.label <- sapply(strsplit(x = row.names(cor.mat), "_"), "[[", 1)
color.v3 <- color.vect[as.factor(row.label)]

heatmap.2(cor.mat, trace="none", #remove trace in cell
          ColSideColors = column.treat, RowSideColors = row.treat, #color bar next to histogram
          density.info = "none", #remove histogram in color key
          key.title = NA, key.xlab = "Correlation",
          #colRow = color.v2, colCol = color.v3,
          col = colorschemes$DarkRedtoBlue.12
          )
### End Plot 1 ####

### Plot 2: pca plot on samples ####

# data preparation
pca.mat <- as.data.frame(t(cnt.mat))
pca.mat$samples <- colnames(cnt.mat)
pca.mat <- separate(pca.mat, samples, into=c("cross", "treatment", "replica"), sep="_", remove=T)
pca.mat$replica <- paste("rep", pca.mat$replica, sep="_")

# PCA
res.pca <- PCA(pca.mat, quali.sup=c(14870:14872), scale.unit = T, graph=F)
fviz_pca_var(res.pca, col.var = "black", geom = "point")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
pca1 <- fviz_pca_ind(res.pca, 
                     geom="point", 
                     col.ind=pca.mat$treatment, 
                     addEllipses = TRUE,
                     palette = c("#cc6633", "#c0c0c0"),
                     legend.title = "Treatments",
                     title="Treatments separation")

pca2 <- fviz_pca_ind(res.pca, 
                     geom="point", 
                     col.ind=pca.mat$replica, 
                     addEllipses = TRUE,
                     palette = c("blue", "green"),
                     legend.title = "Replica",
                     title="Replica separation")


pca1 + theme_cowplot()
pca2 + theme_cowplot()

save_plot("./Figures/pca_treatment.pdf", pca1 + theme_cowplot())
save_plot("./Figures/pca_replica.pdf", pca2 + theme_cowplot())
### END Plot 2 ####
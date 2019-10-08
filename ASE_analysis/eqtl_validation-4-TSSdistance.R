#######################
### eQTL validation ###
#######################

### useful stuff at the beginnings ####

## paths ####
data.path <- "./Data"
fig.path <- "./Figures"
processed.path <- "./Processed"

## libraries ####
library(spatstat.utils) # for fct inside.range
library(tictoc)
library(reshape2)
library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(IRanges)
library(GenomicRanges)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene) # dataset of gene postion
library(compEpiTools) # function TSS
library(flux) # function auc
library(sm) # function sm.density.compare
library(RColorBrewer)
library(cowplot)

## functions #####
get_gene_end <- function(in.mat) {
  out <- tryCatch(
    {
      #model
      m1 <- glmer(cbind(alt, ref) ~ (1|cross), data=in.mat, family=binomial)
      #pval
      coef(summary(m1))
    },
    error=function(cond) {
      return(NA)
    }
  )
}

## colPanels ####
colPanel.infection <- c("#c0c0c0", "#cc6633")

### END useful stuff at the beginnings ####


################################################
### Calculate eQTL position relative to gene ###
################################################

### data loading and pre-processing ####
data.t <- readRDS(file.path(processed.path, "Infected_glmm_c_max6_quantile.rds"))
data.n <- readRDS(file.path(processed.path, "Ctrl_glmm_c_max6_quantile.rds"))

table(data.t$eqtl_type)
table(data.n$eqtl_type)

# select one of the data set
data <- data.t; set = "Infected"
data <- data.n; set = "Ctrl"

# get unique genes
genes <- as.character(unique(data$gene))

### add distance ####
# txdb object
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

# transcript to gene ID dataframe, subset for genes in data frame, get the strand information as well
tx.data <- select(txdb, keys=genes, columns=c("TXNAME", "TXSTRAND"), keytype = "GENEID")

# genes with no transcript name
no.tx <- tx.data[is.na(tx.data$TXNAME),]
no.tx.genes <- unique(no.tx$GENEID)

# list of TSS
tss <- TSS(txdb)
tss.df <- as.data.frame(tss, row.names = seq(1,length(tss),1))

# for each tss, add the gene name
tss.df$geneID <- NA
for (i in 1:nrow(tx.data)) {
  transcript.name <- tx.data[i, "TXNAME"]
  gene.name <- tx.data[i, "GENEID"]
  if (!is.na(transcript.name)) {
    tss.df[tss.df$tx_name == transcript.name, "geneID"] <- gene.name
  }
}

# remove tss with no gene name
tss.df2 <- tss.df[complete.cases(tss.df),]
rm(tss.df)

## compute closest TSS to eQTL ####
tic("computing closest tss")

for (i in 1:nrow(data)) {
  # initation
  if (i == 1) {
    data$closest.tss.distance <- NA
    data$gene.strand <- NA
    noTSS.nb <- 0
  }
  # status update
  if (i %% 500 == 0){
    print(paste(i, "line processed", sep=" "))
  }
  
  # eqtl gene and position
  gene <- data[i, "gene"]
  # if there is no transcript for the gene, skip
  if(gene %in% no.tx.genes) {
    noTSS.nb = noTSS.nb + 1
    next
  }
  eqtl.pos <- data[i, "position"]
  
  # for each eqtl, recover all TSS, put them in a vector
  tmp.df <- tss.df2[tss.df2$geneID == gene,]
  gene.strand <- unique(tmp.df$strand)
  
  # error if different strand
  if(length(gene.strand) != 1) {
    print("--------------------------")
    print("gene strand error")
    print(paste("line nb", i, sep=" "))
    print(tmp.df)
    readline(prompt="press enter")
  }

  # get all tss position
  pos.vect <- as.vector(tmp.df$start)
  
  # distance in bp to closest tss
  distance.vect <- eqtl.pos - pos.vect
  closest.zero <- distance.vect[which.min(abs(distance.vect))]
  
  # write opposite value if on opposite strand
  if(gene.strand == "+") {
    data[i,"closest.tss.distance"] <- closest.zero
  }
  else if(gene.strand == "-") {
    data[i,"closest.tss.distance"] <- -closest.zero
  }
  data[i, "gene.strand"] <- gene.strand
}
print(paste("nb of genes without TSS", noTSS.nb, sep=" "))

toc()

## add sig column for plotting
data$sig <- NA
for (i in 1:nrow(data)) {
  if (is.na(data[i,"postBH.signif"])){
    next
  }
  else if (data[i,"postBH.signif"] == 0) {
    data[i, "sig"] <- "NS"
  }
  else if (data[i,"postBH.signif"] == 1) {
    data[i, "sig"] <- "S"
  }
}
data$sig <- as.factor(data$sig)
# table(data$sig, useNA="always")
# table(data$postBH.signif, useNA="always")

## save ####
saveRDS(data, file.path(processed.path, paste(set, "adjusted_TSSdist.rds", sep="_")))


## clean-up ####
rm(data, data.t, data.n, no.tx, tmp.df, tx.data, tss, closest.zero, distance.vect, eqtl.pos, gene, gene.name, genes, i,
   noTSS.nb, pos.vect, transcript.name, txdb, tss.df2)

### END data processing ####

#################################################
### Plot Distribution of eQTL relative to TSS ###
#################################################

### load data + preprocessing ####
## choose one of the dataset ####
data <- readRDS(file.path(processed.path, "Infected_adjusted_TSSdist.rds")); set <- "Infected"
data <- readRDS(file.path(processed.path, "Ctrl_adjusted_TSSdist.rds")); set <- "Control"

dist.max <- max(data$closest.tss.distance, na.rm=T)
dist.min <- min(data$closest.tss.distance, na.rm=T)
dist.vect <- c(dist.min, dist.max)

## Plot 1: histogram ####
# hist(data.sub$closest.tss.distance, breaks=1000, xlim=c(-20000, 20000), 
#      main="all", xlab="distance to TSS (bp)")

## Plot 2: ggplot density plot with fct density ####

# KS test
p2.stat <- ks.test(data[data$sig=="S", "closest.tss.distance"], data[data$sig=="NS", "closest.tss.distance"], 
                     alternative = "two.sided", exact = T)

# plot
p2  <- ggplot(data, aes(closest.tss.distance, colour=sig, fill=sig)) + 
  geom_density(alpha = 0.4, color = NA) + 
  coord_cartesian(xlim = c(-15000, 15000)) + 
  scale_fill_manual(values = c("orange","blue", "green")) +
  labs(title=paste("distance to TSS", set, sep =" - "), x="distance to TSS in bp", fill="significant", 
       subtitle = paste("KS:", round(p2.stat$p.value, 3), sep=" "))
p2 + theme_cowplot()

# save
save_plot(file.path(fig.path, paste(set, "DistToTSS_All.pdf", sep="_")), p2 + theme_cowplot(), base_aspect_ratio = 1.5)

### END distribution of eQTL ####
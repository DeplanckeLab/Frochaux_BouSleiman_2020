#######################
### eQTL validation ###
#######################

### useful stuff at the beginnings ####

## paths ####
data.path <- "./Data"
fig.path <- "./Figures"
processed.path <- "./Processed"

## libraries ####
library(lme4) # for glmer fct, linear mixed model
library(tictoc)
library(edgeR) # for DGEList
library(limma) # for voom
library(lmerTest)
library(plotrix)

## functions ####
plot.clustering = function(data.input, filename, type){
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

fit_model <- function(in.mat) {
  out <- tryCatch(
    {
      #model
      m1 <- lmer(cnt ~ type + (1|cross), data = in.mat, na.action = "na.omit")
      m1.anova <- anova(m1)
      #pval
      m1.anova$`Pr(>F)`
    },
    error=function(cond) {
      return(NA)
    }
  )
}
## vectors with informations ####
cross.vect <- c("c1x2", "c2x3", "c3x4", "c4x5", "c5x6", "c6x7", "c7x8", "c8x9",
                "c9x10", "c10x11", "c11x12", "c12x13", "c13x14", "c14x15",
                "c17x18", "c18x19", "c19x20", "c20x1")
treat.vect <- c("Pe", "UC")
#type.vect <- c("alt", "ref")
replica.vect <- c("1", "2")

### END useful stuff at the beginnings ####

################
### Analysis ###
################

### Prepare Cnt matrix ####
## read counts
counts <- readRDS(file.path(data.path, "Count_matrix.rds"))
counts.all <- readRDS(file.path(data.path, "Count_matrix_all.rds"))

cor.mat <- round(cor(counts),2)
heatmap(cor.mat)

cor.mat.all <- round(cor(counts.all),2)
heatmap(cor.mat.all, main = "counts only")

cor.dist.matrix = as.dist((1 - cor(counts.all))/2)

# HC
data.hc = hclust(cor.dist.matrix, method="ward.D2")
plot(data.hc, main=paste0("HC on ",type,", WARD"))
data.hc = hclust(cor.dist.matrix)
plot(data.hc, main=paste0("HC on ",type,", COMPLETE"))

## model matrix
model.mat <- expand.grid(cross.vect, treat.vect, replica.vect)
colnames(model.mat) <- c("Lines", "Treatment", "Replica")
rownames(model.mat) <- paste(model.mat$Lines, model.mat$Treatment, model.mat$Replica, sep="_")
model.mat$Lines.Treatment <- paste(model.mat$Lines, model.mat$Treatment, sep=".")
model.mat <- model.mat[!rownames(model.mat) %in% c("c20x1_Pe_2"),]
saveRDS(model.mat, file.path(processed.path, "model_matrix.rds"))

# make dge list
counts.all$c20x1_Pe_2 <- NULL
data.dge = DGEList(counts = counts.all)

# filtering
data.dge = data.dge[rowSums(log2(cpm(data.dge)+1)) > 0.5, ]
data.filtered = DGEList(counts = counts.all[rownames(data.dge$counts), ])
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
plot.clustering(data.voom$E, file.path(fig.path, "voom_clustering_all_BatchRemoved.png"), "VOOM")

saveRDS(data.voom$E, file.path(processed.path, "Normalized_Cnt_Matrix.rds"))
### load data ####

## Choose data set ####
eqtl.Inf <- readRDS(file = file.path(processed.path, "Infected_glmm_c_max6_quantile.rds"))
eqtl.ctrl <- readRDS(file = file.path(processed.path, "Ctrl_glmm_c_max6_quantile.rds"))

data <- eqtl.Inf; set <- "Infected"
data <- eqtl.ctrl; set <- "Ctrl"

# Normalized Count matrix
brb.cnt.voom <- readRDS(file.path(processed.path, "Normalized_Cnt_Matrix.rds"))

## Main loop ####

#print(paste("cutoff =", cutoff, sep=" "))
for (line.it in 1:nrow(data)) {
  # initiate on loop start
  if (line.it == 1) {
    tic("Measurment of true trans")
    data$trans_pval <- NA
    
    # create error output
    err.output <- data[,1:4]
    err.output$error <- NA
  }
  
  # print(paste("line:", line.it, sep=" "))
  
  # loop information
  if (line.it %% 100 == 0) {
    print(line.it)
  }
  
  # 1. get the number of homozygous cross
  hom_alt.nb <- data[line.it, "nb_Hom_alt"]
  hom_ref.nb <- data[line.it, "nb_Hom_ref"]
  nb.Hom.cross <- hom_ref.nb + hom_alt.nb
  
  # if there is not at least on of each type of homozygote, skip
  if (hom_alt.nb <= 2 & hom_ref.nb <= 2) {
    err.output[line.it, "error"] <- "No_geno"
    next
  }
  
  # print("stage 1 complete")
  
  # 2. create matrix for analysis
  # get the all the homozygous cross
  sub.line <- as.vector(data[line.it, grep(pattern = "_status", colnames(data))])
  # sub.line
  # sub.line.ref <- sub.line[,which(sub.line == "Hom_ref")]
  # sub.line.ref
  # sub.line.alt <- sub.line[,which(sub.line == "Hom_alt")]
  # sub.line.alt
  sub.line2 <- sub.line[,which(sub.line == "Hom_ref" | sub.line == "Hom_alt")]
  sub.line2
  #rm(sub.line, sub.line.alt, sub.line.ref)
  
  # make the matrix
  eqtl.df <- as.data.frame(t((sub.line2)))
  eqtl.df[,2] <- row.names(eqtl.df)
  row.names(eqtl.df) <- c()
  colnames(eqtl.df) <- c("type", "cross")
  eqtl.df <- rbind(eqtl.df, eqtl.df)
  eqtl.df$type <- gsub("Hom_", "", eqtl.df$type)
  eqtl.df$cross <- gsub("_status", "", eqtl.df$cross)
  eqtl.df$replica <- c(rep("rep_1", length(sub.line2)), rep("rep_2", length(sub.line2)))
  
  # get the gene
  eqtl.gene <- as.character(data[line.it, "gene"])
  # check if gene is in the cnt matrix
  if (!eqtl.gene %in% row.names(brb.cnt.voom)) {
    err.output[line.it, "error"] <- "no_gene_cntMat"
    next
  }
  
  # print("matrix ok")
  
  # fill the matrix
  for (i in 1:nrow(eqtl.df)) {
    if (i == 1) { eqtl.df$cnt <- NA}
    df.cross <- eqtl.df[i, "cross"]
    df.replica <- eqtl.df[i, "replica"]
    if (set == "Infected") {df.replica <- gsub("rep", "Pe", df.replica)}
    if (set == "Ctrl") {df.replica <- gsub("rep", "UC", df.replica)}
    sample <- paste(df.cross, df.replica, sep = "_")
    if (sample == "c20x1_Pe_2") { value <- NA }
    else if (sample == "c19x20_Pe_1") { value <- NA }
    else { value <- brb.cnt.voom[eqtl.gene, sample] }
    eqtl.df[i, "cnt"] <- value
  }
  
  # print("fill ok")
  
  # # fit two models, difference between them is the result of the influence of alt and ref vs luck
  # m1 <- lmer(cnt ~ type + (1|cross), data = eqtl.df, na.action = "na.omit")
  # #m2 <- lmer(cnt ~ (1|cross), data = eqtl.df, na.action = "na.omit")
  # #lm.res <- anova(m1,m2)
  # lm.res <- anova(m1)
  p.val <- fit_model(eqtl.df)
  
  # add value to the data
  data[line.it, "trans_pval"] <- p.val
  err.output[line.it, "error"] <- "OK"
}

toc()

pdf(file = file.path(fig.path, paste(set, "min2ofeach.pdf", sep="_")))
par(mfrow=c(1,3))
hist(data$trans_pval, breaks=100, xlab = "P-value", ylab = "Count", main="Raw p-val, min 2 of each")
nb <- sum(data$trans_pval < 0.05, na.rm = T)
hist(p.adjust(data$trans_pval, method = "BH"), breaks=100, xlab = "P-value", ylab = "Count", 
     main=paste("BH p-val - ", as.character(nb), "sig", sep=" "))
barplot(table(err.output$error, useNA = "always"), las=2)
dev.off()

# keep important columns
data <- data[,-grep(pattern = "x", colnames(data))]

# add BH corrected p-values
data$trans_BH <- p.adjust(data$trans_pval, method = "BH")

saveRDS(data, file = file.path(processed.path, paste(set, "cis_trans.rds", sep="_")))
saveRDS(err.output, file = file.path(processed.path, paste(set, "trans_measurement_1model_Errors.rds", sep="_")))
## END ####

### Plot nb hom ####
hist(data$nb_Hom_ref)
abline(v=2, col="red", xlab="nb Hom ref", main="Nb Homozygous reference")
text(5, 4000, as.character(sum(data$nb_Hom_ref > 2, na.rm=T)), col="red")

hist(data$nb_Hom_ref)
abline(v=2, col="red", xlab="nb Hom ref", main="Nb Homozygous reference")
text(10, 1500, as.character(sum(data$nb_Hom_ref > 2, na.rm=T)), col="red")

sum(data$nb_Hom_alt > 1 & data$nb_Hom_ref > 2, na.rm=T)

### END ####


## PLot: barplot of repartition of tested and validated eQTLs Trans####
# get stats
error.n <- readRDS(file.path(processed.path, "Ctrl_trans_measurement_1model_Errors.rds"))
error.t <- readRDS(file.path(processed.path, "Infected_trans_measurement_1model_Errors.rds"))
data.n <- readRDS(file = file.path(processed.path, "Ctrl_cis_trans.rds"))
data.t <- readRDS(file = file.path(processed.path, "Infected_cis_trans.rds"))

table(data.t$trans_BH < 0.05, useNA="ifany")
table(data.n$trans_BH < 0.05, useNA="ifany")

summary.err.n <- table(error.n$error)
summary.err.t <- table(error.t$error)

summary.signif.t <- table(data.t$trans_BH < 0.05, useNA="ifany")
summary.signif.n <- table(data.n$trans_BH < 0.05, useNA="ifany")


error.df <- as.data.frame(rbind(summary.err.n, summary.err.t))
row.names(error.df) <- c("control", "infected")

# add line with number of validated eqtl

table.combined2 <- as.data.frame(rbind(summary.signif.n, summary.signif.t))
row.names(table.combined2) <- c("control", "infected")



combined1.melt <- melt(as.matrix(error.df))
colnames(combined1.melt) <- c("treatement", "variable", "all")
significant.n <- 1441 #in summary.signif.n
significant.t <- 2112
ctrl.NT <- 6459
inf.NT <- 6606
combined1.melt$signif <- 0
combined1.melt[combined1.melt$treatement=="control" & combined1.melt$variable=="ok", "signif"] <- significant.n
combined1.melt[combined1.melt$treatement=="infected" & combined1.melt$variable=="ok", "signif"] <- significant.t


simple.df <- as.data.frame(c("Control", "Infected", "Control", "Infected"))
colnames(simple.df) <- "treatement"
simple.df$variable <- c("OK", "OK", "NT", "NT")
simple.df$variable <- as.factor(simple.df$variable)
simple.df$all <- c(10343+727, 9944+1019,4781,4960)
simple.df$signif <- c(727,1019,0,0)

saveRDS(simple.df, file=file.path(data.path, "simple_trans.rds"))
simple.df <- readRDS(file.path(data.path, "simple_trans.rds"))
simple.df

p3 <- ggplot(data=simple.df, aes(x=variable)) +
  geom_bar(aes(y=all, fill=variable), stat="identity", position="dodge", alpha=0.5) +
  geom_bar(aes(y=signif, fill=variable), stat="identity", position="dodge") +
  facet_wrap(~treatement) +
  scale_fill_manual(values=c("red4", "darkgreen")) +
  labs(title="Distribution of eQTL based on number of Heterozygote cross",
       x="tested", y="Number of eQTL") +
  theme(legend.position="none")

p3 + theme_cowplot()
p3 + theme_cowplot() + theme(legend.position="none")

save_plot(file.path(fig.path, "errors_trans.pdf"), p3 + theme_cowplot(), 1.2)


# get numbers
727/(10343+727)*100
1019/(9944+1019)*100

### END ANALYSIS ####
library(edgeR)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicFeatures)
library(Mfuzz)
library(GOstats)
library(ggplot2)
library(Cairo)
library(MatrixEQTL)
library(plyr)
library(dplyr)
library(lattice)
library(sm)
source("./R/commonFunctions.R") # adapt based on the actual file location
##load experiment design data
targets <- read.table("./Data/design.txt", header=TRUE)
targets$condres = targets$condition:targets$resistance
counts <- read.table("./Data/count.table.formatted.onlyGenes.txt", header=TRUE, row.names=1)
colnames(counts) <- paste(targets$line,
                          substr(targets$condition, 1, 1),
                          substr(targets$resistance, 1, 1),
                          sep="_")

condres.colors = c("#c6e6ac", "#ff9595", "#529021","#c00000")
condres.colors.mapped = targets$condres
levels(condres.colors.mapped) = condres.colors
condres.colors.mapped = as.character(condres.colors.mapped)
resistance.colors = c("#529021","#c00000")
##load gene IDs
gene.ids = read.table("./Data/genesBDGP5.25.2014-05-23.geneIDsandNames.txt", header=FALSE) # prepared separately
colnames(gene.ids)<- c("FBgnID", "Gene.Name")
row.names(gene.ids)=gene.ids$FBgnID

#load txdb
txdb = loadDb("./Data/genesBDGP5.25.2014-05-23.gtf.txdb")
transcripts = transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
transcripts$gene_id = as.character(transcripts$gene_id)
genes = genes(txdb, columns="gene_id")

##Voom normalization of counts, followed by differential expression analysis
isexpr <- rowSums(counts >= 5) >= round(0.5*ncol(counts))
targets$condition
mod.mat = model.matrix( ~ condition  + condition:resistance, data=targets)
table(isexpr)
c <- counts[isexpr,]
nf <- calcNormFactors(c)
y.voom =  voom(c, mod.mat, plot=FALSE, normalize = "quantile", lib.size=colSums(c)*nf) # first round of voom : https://stat.ethz.ch/pipermail/bioconductor/attachments/20130530/4dcc9475/attachment.pl
corfit.y = duplicateCorrelation(y.voom, mod.mat, block = targets$line)
y.voom = voom(c, mod.mat, block = targets$line, correlation =
                corfit.y$consensus, normalize="quantile",lib.size=colSums(c)*nf, plot = FALSE) #second round of voom
fit = lmFit(y.voom, mod.mat, correlation = corfit.y$consensus.correlation, block = targets$line)
fit <- eBayes(fit)
summary(decideTests(fit, p.value = 0.05, lfc = 1) )
tt.n=topTable(fit, coef = "conditionNaive:resistancesusceptible", n = Inf)
tt.n$gene.name = gene.ids$Gene.Name[match(rownames(tt.n), gene.ids$FBgnID)]
head(tt.n, n=20)
write.table(tt.n, file = "./Data/tt.resistance.naive.txt")
tt.t=topTable(fit, coef = "conditionTreated:resistancesusceptible", n = Inf)
tt.t$gene.name = gene.ids$Gene.Name[match(rownames(tt.t), gene.ids$FBgnID)]
head(tt.t, n=20)
write.table(tt.t, file = "./Data/tt.resistance.treated.txt")
ttStatus=topTable(fit, coef = "conditionTreated", n = Inf)
ttStatus$gene.name = gene.ids$Gene.Name[match(rownames(ttStatus), gene.ids$FBgnID)]
write.table(ttStatus, file = "./Data/tt.condition.txt")
head(ttStatus, n=20)
save(y.voom, file = './Data/y.voom')

#PCA 
library(FactoMineR)
pca <- PCA(t(y.voom$E[,]), graph = FALSE)
pdf(file = "./Plots/cpm.pca.76.samples.pdf",width = 5, height = 5, pointsize = 0.6)

plot(pca$ind$coord[,1:2], col = ifelse(targets$resistance == "resistant","green","red"), pch = ifelse(targets$condition == 'Naive', 21,22), cex=3, lwd=2)
dev.off()
pca.naive = PCA(t(y.voom$E[,1:38]), graph = FALSE)
pdf(file = "./Plots/cpm.pca.38.samples.naive.pdf",width = 5, height = 5, pointsize = 0.6)

plot(pca.naive$ind$coord[,1:2], col = ifelse(targets$resistance[1:38] == "resistant","green","red"), pch = ifelse(targets$condition[1:38] == 'Naive', 21,22), cex=3, lwd=2)
dev.off()
pca.treated = PCA(t(y.voom$E[,39:76]), graph = FALSE)
pdf(file = "./Plots/cpm.pca.38.samples.treated.pdf",width = 5, height = 5, pointsize = 0.6)

plot(pca.treated$ind$coord[,1:2], col = ifelse(targets$resistance[1:38] == "resistant","green","red"), pch = ifelse(targets$condition[39:76] == 'Naive', 21,22), cex=3, lwd=2)
dev.off()


#### eQTL data loading ####
##now get the eQTL data
##first load the SNP data
SNP_file_name = "./Data/genoForMatrixeQTL.noS14.noS54.txt";
snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 5,000 rows
snps$LoadFile(SNP_file_name);

eQTL.result.files = list.files("./Data", "withCovariates.cis", full.names = TRUE)
eQTL.result.files <- eQTL.result.files[grepl("Naive|Treated", eQTL.result.files)]
eQTL.analysis.title = gsub("./Data/MatrixeQTL.", "", eQTL.result.files)
eQTL.analysis.title = gsub(".output.withCovariates.cis", "", eQTL.analysis.title)
eQTL.result = read.table(eQTL.result.files[1], header=TRUE)
eQTL.result = eQTL.result[eQTL.result$FDR < 0.05,]
eQTL.result$analysis.title = eQTL.analysis.title[1]
for(i in c(2,3)){
  toADD = read.table(eQTL.result.files[i], header=TRUE)
  toADD = toADD[toADD$FDR < 0.05,]
  toADD$analysis.title = eQTL.analysis.title[i]
  eQTL.result = rbind(eQTL.result, toADD)
}
head(eQTL.result)
length(unique(eQTL.result$gene[eQTL.result$analysis.title == "Naive"]))
eQTL.result$SNP = as.character(eQTL.result$SNP)
eQTL.result$chrom <- gsub(":(.)+$", "",eQTL.result$SNP)
eQTL.result$pos <- as.numeric(gsub("^(.)+:", "", gsub("_(.)+$", "", eQTL.result$SNP)))
eQTL.result$pos <- as.numeric(gsub("^(.)+:", "", gsub("_(.)+$", "", eQTL.result$SNP)))
#add distance from TSS

library(BiocGenerics)
library(GenomicRanges)
tss = genes
start(tss[strand(tss) == "-"]) = end(tss[strand(tss) == "-"])
end(tss[strand(tss) == "+"]) = start(tss[strand(tss) == "+"])
tes = genes
start(tes[strand(tes) == "+"]) = end(tes[strand(tes) == "+"])
end(tes[strand(tes) == "-"]) = start(tes[strand(tes) == "-"])

eQTL.gr = GRanges(seqnames = eQTL.result$chrom, ranges = IRanges(start = as.numeric(eQTL.result$pos) , width = 1 ) , 
                  strand = "*",
                  analysis = eQTL.result$analysis.title,
                  ID = eQTL.result$SNP,
                  gene = eQTL.result$gene,
                  gene_name = gene.ids$Gene.Name[match(eQTL.result$gene,gene.ids$FBgnID)],
                  beta = eQTL.result$beta,
                  t.stat = eQTL.result$t.stat,
                  p.value = eQTL.result$p.value,
                  FDR = eQTL.result$FDR,
                  Treatment.FC = ttStatus[eQTL.result$gene,"logFC"],
                  Treatment.FDR = ttStatus[eQTL.result$gene, "adj.P.Val"],
                  Resistance.FC.n = tt.n[eQTL.result$gene,"logFC"],
                  Resistance.FDR.n = tt.n[eQTL.result$gene, "adj.P.Val"],
                  Resistance.FC.t = tt.t[eQTL.result$gene,"logFC"],
                  Resistance.FDR.t = tt.t[eQTL.result$gene, "adj.P.Val"],
                  gene.pos = start(tss[eQTL.result$gene]),
                  gene.end = end(tes[eQTL.result$gene]),
                  gene.strand = as.character(strand(tss[eQTL.result$gene])))
eQTL.gr$Group <- ifelse(eQTL.gr$ID %in% intersect(eQTL.gr$ID[eQTL.gr$analysis == "Naive"], eQTL.gr$ID[eQTL.gr$analysis == "Treated"]), "Shared", "Unique")
eQTL.gr$Group[eQTL.gr$Group == "Unique"] <- paste0(eQTL.gr$Group[eQTL.gr$Group == "Unique"], "-",eQTL.gr$analysis[eQTL.gr$Group == "Unique"])
eQTL.gr$distance.from.gene.start = start(eQTL.gr) - eQTL.gr$gene.pos
eQTL.gr$distance.from.gene.start[eQTL.gr$gene.strand == "-"] = - (start(eQTL.gr) - eQTL.gr$gene.pos)[eQTL.gr$gene.strand == "-"]
eQTL.gr$distance.from.gene.end = start(eQTL.gr) - eQTL.gr$gene.end
eQTL.gr$distance.from.gene.end[eQTL.gr$gene.strand == "-"] = - (start(eQTL.gr) - eQTL.gr$gene.end)[eQTL.gr$gene.strand == "-"]
eQTL.gr$locationWrtGeneBody <- ifelse(sign(eQTL.gr$distance.from.gene.start) == sign(eQTL.gr$distance.from.gene.end), "Outside", "Within")
eQTL.gr$locationWrtGeneBody[eQTL.gr$locationWrtGeneBody == "Outside"] <- ifelse(eQTL.gr$distance.from.gene.start[eQTL.gr$locationWrtGeneBody=="Outside"] > 0, "Downstream", "Upstream")
table(eQTL.gr$locationWrtGeneBody)
eQTL.gr$GroupCompound <- paste0(eQTL.gr$Group,eQTL.gr$locationWrtGeneBody)
eQTL.gr$withinGene <- grepl("Within", eQTL.gr$locationWrtGeneBody)
eQTL.gr$gene.length <- width(genes)[match(as.character(eQTL.result$gene), genes$gene_id)]
#now calculate a metaplot value
normlength <- mean(eQTL.gr$gene.length[match(eQTL.gr$gene, unique(eQTL.gr$gene))])
normlength <- 10000
eQTL.gr$metaPlotValue <- ifelse(eQTL.gr$locationWrtGeneBody == "Upstream", eQTL.gr$distance.from.gene.start, eQTL.gr$distance.from.gene.end+normlength)
eQTL.gr$metaPlotValue[eQTL.gr$locationWrtGeneBody == "Within"] <- normlength*eQTL.gr$distance.from.gene.start[eQTL.gr$locationWrtGeneBody == "Within"]/eQTL.gr$gene.length[eQTL.gr$locationWrtGeneBody == "Within"]


### sample random variants from SNPS in a 10kb window on both sides
sampleFrom <-genes[as.character(unique(eQTL.gr$gene))]
start(sampleFrom) <- start(sampleFrom) -10000
end(sampleFrom) <- end(sampleFrom) +10000
snpsloc <- read.table("./Data/snpsloc.txt", header = T)
snpsloc$chr=as.character(plyr::revalue(as.character(snpsloc$chr), c("1"="2L", 
                                                                    "2"="2R",
                                                                    "3"="3L",
                                                                    "4"="3R",
                                                                    "5"="4",
                                                                    "6"="X")))

snpsloc.gr <- GRanges(seqnames = snpsloc$chr, ranges = IRanges(start = as.numeric(snpsloc$pos), width = 1))
snpsloc.gr$ID <- as.character(snpsloc$snp)
#annotate the snpsloc.gr with gene IDs
overlaps <- findOverlaps(snpsloc.gr, sampleFrom)
snpsloc.gr <- snpsloc.gr[queryHits(overlaps)]
snpsloc.gr$gene_id <- sampleFrom$gene_id[subjectHits(overlaps)]
snpsloc.gr$gene.pos = start(tss[snpsloc.gr$gene_id])
snpsloc.gr$gene.end = end(tes[snpsloc.gr$gene_id])
snpsloc.gr$gene.strand = as.character(strand(tss[snpsloc.gr$gene_id]))
snpsloc.gr$Group <- "Random"
snpsloc.gr$distance.from.gene.start = start(snpsloc.gr) - snpsloc.gr$gene.pos
snpsloc.gr$distance.from.gene.start[snpsloc.gr$gene.strand == "-"] = - (start(snpsloc.gr) - snpsloc.gr$gene.pos)[snpsloc.gr$gene.strand == "-"]
snpsloc.gr$distance.from.gene.end = start(snpsloc.gr) - snpsloc.gr$gene.end
snpsloc.gr$distance.from.gene.end[snpsloc.gr$gene.strand == "-"] = - (start(snpsloc.gr) - snpsloc.gr$gene.end)[snpsloc.gr$gene.strand == "-"]
snpsloc.gr$locationWrtGeneBody <- ifelse(sign(snpsloc.gr$distance.from.gene.start) == sign(snpsloc.gr$distance.from.gene.end), "Outside", "Within")
snpsloc.gr$locationWrtGeneBody[snpsloc.gr$locationWrtGeneBody == "Outside"] <- ifelse(snpsloc.gr$distance.from.gene.start[snpsloc.gr$locationWrtGeneBody=="Outside"] > 0, "Downstream", "Upstream")
snpsloc.gr$withinGene <- grepl("Within", snpsloc.gr$locationWrtGeneBody)
table(snpsloc.gr$locationWrtGeneBody)
snpsloc.gr$gene.length <- width(genes)[match(as.character(snpsloc.gr$gene_id), genes$gene_id)]
#now calculate a metaplot value, ranging from -10000 to +10000 for upstream and downstream, respectively
snpsloc.gr$metaPlotValue <- ifelse(snpsloc.gr$locationWrtGeneBody == "Upstream", snpsloc.gr$distance.from.gene.start, snpsloc.gr$distance.from.gene.end+normlength)
snpsloc.gr$metaPlotValue[snpsloc.gr$locationWrtGeneBody == "Within"] <- normlength*snpsloc.gr$distance.from.gene.start[snpsloc.gr$locationWrtGeneBody == "Within"]/snpsloc.gr$gene.length[snpsloc.gr$locationWrtGeneBody == "Within"]


# pull random snps
l<-length(eQTL.gr[eQTL.gr$analysis == "Treated"])
s <- lapply(as.list(1:100), function(x){
  s.gr<-snpsloc.gr[sample(1:length(snpsloc.gr), l)]
  s.gr$Group <- paste0("Random-",x)
  s.gr
})
s <- do.call("c", s)
s$GroupCompound <- paste0(s$Group, s$locationWrtGeneBody)
infection.colors = c("#c0c0c0", "#cc6633")

pdf(file = "./Plots/eQTL.metaplot.pdf", width = 4, height = 4) #not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(aes(x=distance.from.gene.start, y=..density.., col = analysis)) +
  geom_vline(xintercept = 5) +
  scale_x_continuous(limits = c(-10000, 10000)) + #scale_y_continuous(limits = c(4, 5)) +
  scale_color_manual(values= infection.colors) +
  scale_size_identity() + theme_light() + xlab("Distance from TSS") + ylab("Density") + theme(legend.position="none")
dev.off()

pdf(file = "./Plots/eQTL.metaplot.TES.pdf", width = 4, height = 4) # not included in paper

par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(aes(x=distance.from.gene.end, y=..density.., col = analysis)) +
  geom_vline(xintercept = 5) +
  scale_x_continuous(limits = c(-10000, 10000)) + #scale_y_continuous(limits = c(4, 5)) +
  scale_color_manual(values= infection.colors) +
  scale_size_identity() + theme_light() + xlab("Distance from TES") + ylab("Density") + theme(legend.position="none")
dev.off()

pdf(file = "./Plots/eQTL.metaplot.infectionSpecific.pdf", width = 7, height = 4) # not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(aes(x=distance.from.gene.start, y=..density.., col = Group)) +
  geom_vline(xintercept = 5) +
  scale_x_continuous(limits = c(-10000, 10000)) +
  scale_size_identity() + theme_light() + xlab("Distance from TSS") + ylab("Density")
dev.off()

pdf(file = "./Plots/eQTL.metaplot.infectionSpecific.count.pdf", width = 7, height = 4) # not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(aes(x=distance.from.gene.start, y=..count.., col = Group)) +
  geom_vline(xintercept = 5) +
  scale_x_continuous(limits = c(-10000, 10000)) + 
  scale_size_identity() + theme_light() + xlab("Distance from TSS") + ylab("Count")
dev.off()

pdf(file = "./Plots/eQTL.metaplot.infectionSpecific.TES.pdf", width = 7, height = 4) # not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(aes(x=distance.from.gene.end, y=..density.., col = Group)) +
  geom_vline(xintercept = 5) +
  scale_x_continuous(limits = c(-10000, 10000)) + 
  scale_size_identity() + theme_light() + xlab("Distance from TES") + ylab("Density")
dev.off()

pdf(file = "./Plots/eQTL.metaplot.infectionSpecific.wrtGene.pdf", width = 4, height = 4) # included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(data = as.data.frame(values(s)), aes(x=metaPlotValue, y=..density.., group = Group), trim = FALSE, col = "grey80", alpha = 0, size = 0.2, weight = 0.2) + 
  geom_density(aes(x=metaPlotValue, y=..density.., col = Group, group = Group), trim = FALSE) +
  geom_vline(xintercept = c(0,normlength), lty = 2) +
  scale_x_continuous(limits = c(-10000, 10000 + normlength)) + 
  scale_size_identity() + theme_light() + xlab("Position") + ylab("Density") + theme(legend.position="bottom")
dev.off()


pdf(file = "./Plots/eQTL.metaplot.infectionSpecific.wrtGene.count.pdf", width = 4, height = 4) # not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_density(data = as.data.frame(values(s)), aes(x=metaPlotValue, y=..count.., group = Group), trim = FALSE, col = "grey80", alpha = 0, size = 0.2, weight = 0.2) + 
  geom_density(aes(x=metaPlotValue, y=..count.., col = Group, group = Group), trim = FALSE) +
  geom_vline(xintercept = c(0,normlength), lty = 2) +
  scale_x_continuous(limits = c(-10000, 10000 + normlength)) + #scale_y_continuous(limits = c(4, 5)) +
  scale_size_identity() + theme_light() + xlab("Position") + ylab("Count") + theme(legend.position="bottom")
dev.off()

##also a metaplot of the p-values
library(splines)
library(MASS)
pdf(file = "./Plots/eQTL.metaplot.association.strength.pdf", width = 4, height = 4) # not included in paper
par(mfrow=c(1,1))
ggplot(data = as.data.frame(values(eQTL.gr[eQTL.gr$analysis != "Interaction"]))) +
  geom_point(aes(x=distance.from.gene.start, y=-log10(p.value), col = analysis), size = 1, stroke = 0, alpha = 0.1) +
  geom_smooth(formula = y ~ x, se = FALSE, aes(x=distance.from.gene.start, y=-log10(p.value), col = analysis)) +
  scale_x_continuous(limits = c(-10000, 10000)) + scale_y_continuous(limits = c(4, 5)) +
  scale_color_manual(values= infection.colors) + 
  scale_size_identity() + theme_light()  + theme(legend.position="none")
dev.off()


snps.matrix = as.matrix(snps)[eQTL.gr$ID,]
rownames(snps.matrix) = NULL
values(eQTL.gr) = cbind(values(eQTL.gr), as.data.frame(snps.matrix))
table(eQTL.gr$analysis)
ddply(as.data.frame(values(eQTL.gr)), .(analysis), summarize, nQTLs = length(unique(ID)), 
      nGenes = length(unique(gene)))
length(intersect(unique(eQTL.gr$ID[eQTL.gr$analysis == "Treated"]), unique(eQTL.gr$ID[eQTL.gr$analysis == "Naive"]) )) #common eQTLs

common.genes = intersect(unique(eQTL.gr$gene[eQTL.gr$analysis == "Treated"]), unique(eQTL.gr$gene[eQTL.gr$analysis == "Naive"]))
treated.genes = setdiff(unique(eQTL.gr$gene[eQTL.gr$analysis == "Treated"]), common.genes)
naive.genes =  setdiff(unique(eQTL.gr$gene[eQTL.gr$analysis == "Naive"]), common.genes)
universe = rownames(y.voom$E)
# perform GO enrichment for the three groups here
goresults <- lapply(c("common", "naive", "treated"), function(i){
    if(i == "common"){
      s = common.genes
    }else if(i == "naive"){
      s= naive.genes
    }else if(i == "treated"){
      s= treated.genes
    }
    GO.analyze(universe=universe,
               selected= s,
               ontology="BP",
               report=F,
               pv = 0.01,
               short.name=paste0("GO.BP.eQTL.",i),
               outputfile=paste0("GO.BP.eQTL.",i),
               kable.output=FALSE)

})
names(goresults) <- c("common", "naive", "treated")
saveRDS(goresults, file = "./Data/eQTL.GO.results.RDS")


goresults_summary <- do.call(rbind, mapply(x=goresults, y =names(goresults), function(x,y){
  out <- summary(x)
  out$Analysis <- y
  out
}, SIMPLIFY = F))
goresults_summary_revigo <- goresults_summary[,c("GOBPID","Pvalue", "Count", "Analysis")]
write.table(goresults_summary_revigo, file = './Reports/GO.eQTLs.pooled.for.revigo.txt', row.names = F, sep = "\t", quote = F)
#' the GO results are then sent to REVIGO: http://revigo.irb.hr/
#' Settings:
#' Allowed similarity = 0.7
#' The numbers associated to GO categoris are: p-values
#' Select a database with GO term sizes: Drosophila melanogaster
#' Select a semantic similarity measure to use: SimRel
#' The R script is downloaded and customized and found at: "./R/revigo.plot.eqtls.R

source("./R/revigo.plot.eqtls.R") # adapt path to your needs


##plot some genes - examples
gene.ids[gene.ids$Gene.Name == "ntc",]
eQTL.gr[eQTL.gr$gene == "FBgn0035461"]
# general.plot.eQTL(which(eQTL.gr$ID == "3L.3802119.D.1.I.1"))
# general.plot.eQTL(8140)
general.plot.eQTL(which(eQTL.gr$ID == "3L:3814085_G/T")) ### plot for ntc br eqtl
general.plot.eQTL(which(eQTL.gr$ID == "3L:3814115_T/C")) ### 
general.plot.eQTL(which(eQTL.gr$ID == "3L:3802119_A/T")) ### 
eQTL.gr[eQTL.gr$ID == "2L:157831_C/A"]
general.plot.eQTL(which(eQTL.gr$ID == "2L:157831_C/A")[1])
eQTL.gr[eQTL.gr$gene == "FBgn0031228"]
general.plot.eQTL(which(eQTL.gr$ID == "3R:13220527_G/C")[1])
general.plot("FBgn0031228")

save(eQTL.gr, file = "./Data/eQTL.gr")



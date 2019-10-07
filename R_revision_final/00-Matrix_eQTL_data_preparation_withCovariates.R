#############################################
##### Data preparation for eQTL analysis ####
#############################################

#' Files required for the analysis:
#' GTF file (BDGP5.25). Can be obtained from: http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Drosophila_melanogaster/Ensembl/BDGP5.25/Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
#' DGRP freeze 2 VCF file. Can be obtained from: http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf
#' Wolbachia infection status. Can be obtained from: http://dgrp2.gnets.ncsu.edu/data/website/wolbachia.xlsx 
#' RNAseq count table. Can be obtained from GEO GSE118142: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118142
#' Design matrix of the experiment: a tab-delimited file with information about strain and infection status
#' 


#############################################
#####   Gene locations and expression   #####
#############################################
library(GenomicFeatures)
if(length(list.files(path = "./Data/", pattern = ".txdb")) ==0){
  txdb <- makeTranscriptDbFromGFF("./Data/genesBDGP5.25.2014-05-23.gtf", format="gtf")
  saveDb(txdb)
}else{
  txdb = loadDb("./Data/genesBDGP5.25.2014-05-23.gtf.txdb")
}

genes <- genes(txdb)
genes<-as.data.frame(genes)[,-4]
head(genes)
##reformat as such: geneid chr s1 s2
genelocs <- data.frame(geneid=genes$gene_id, chr=genes$seqnames, s1=genes$start, s2=genes$end)
head(genelocs)
levels(genelocs$chr)
library(plyr)
genelocs$chr=as.character(revalue(as.character(genelocs$chr), c("2L"="1", 
                                                                "2R"="2",
                                                                "3L"="3",
                                                                "3R"="4",
                                                                "4"="5",
                                                                "X"="6")))

write.table(genelocs, file="./Data/genlocs.txt", col.names=T, sep="\t", quote=FALSE, row.names=FALSE)


#### prepare cpm table
library(limma)
library(edgeR)
library(rgl)
counts=read.table("./Data/count.table.formatted.onlyGenes.txt", header=TRUE)
rownames(counts)=counts[,1]
counts=counts[,-1]
head(counts)
for (i in 1:ncol(counts)){
  counts[,i]=as.numeric(counts[,i])
}

##first both conditions
targets <- read.table("./Data/design.txt", header=TRUE)
counts = counts
colnames(counts) <- targets$sampleID[match(colnames(counts), targets$sample)]
isexpr <- rowSums(counts >= 5) >= round(0.9*ncol(counts))
mod.mat = model.matrix( ~ condition + resistance, data=targets)
table(isexpr)
c <- counts[isexpr,]
nf <- calcNormFactors(c)
y.voom =  voom(c, mod.mat, plot=FALSE, normalize = "quantile", lib.size=colSums(c)*nf) # first round of voom : https://stat.ethz.ch/pipermail/bioconductor/attachments/20130530/4dcc9475/attachment.pl
corfit.y = duplicateCorrelation(y.voom, mod.mat, block = targets$line)
y.voom = voom(c, mod.mat, block = targets$line, correlation = corfit.y$consensus, 
              normalize="quantile",lib.size=colSums(c)*nf, plot = FALSE) #second round of voom

mds<-plotMDS(y.voom, labels=rownames(y.voom$samples), col=as.numeric(targets$resistance), ndim=10, top=500, gene.selection = "pairwise")
write.table(y.voom$E, file="./Data/cpm.table.76.lines.txt")
write.table(y.voom$weights, file="./Data/weights.table.76.lines.txt")

##now by condition
#naive first
targets <- read.table("./Data/design.txt", header=TRUE)[-c(39:76),]
c <- counts[,-c(39:76)]
isexpr <- rowSums(c >= 5) >= round(0.9*ncol(c))
table(isexpr)
c <- c[isexpr,]
nf <- calcNormFactors(c)
mod.mat = model.matrix( ~ resistance, data=targets)
y.voom =  voom(c, mod.mat, plot=FALSE, normalize = "quantile", lib.size=colSums(c)*nf) # first round of voom : https://stat.ethz.ch/pipermail/bioconductor/attachments/20130530/4dcc9475/attachment.pl
write.table(y.voom$E, file="./Data/cpm.table.Naive.txt")
write.table(y.voom$weights, file="./Data/weights.table.Naive.lines.txt")

#Treated next
targets <- read.table("./Data/design.txt", header=TRUE)[-c(1:38),]
c <- counts[,-c(1:38)]
isexpr <- rowSums(c >= 5) >= round(0.9*ncol(c))
table(isexpr)
c <- c[isexpr,]
nf <- calcNormFactors(c)
mod.mat = model.matrix( ~ resistance, data=targets)
y.voom =  voom(c, mod.mat, plot=FALSE, normalize = "quantile", lib.size=colSums(c)*nf) # first round of voom : https://stat.ethz.ch/pipermail/bioconductor/attachments/20130530/4dcc9475/attachment.pl
write.table(y.voom$E, file="./Data/cpm.table.Treated.txt")
write.table(y.voom$weights, file="./Data/weights.table.Treated.lines.txt")

#############################################
#####   Genotypes and variant location  #####
#############################################
library(VariantAnnotation)
library(BSgenome.Dmelanogaster.UCSC.dm3)
targets <- read.table("./Data/design.txt", header=TRUE)[,]
genome <- BSgenome.Dmelanogaster.UCSC.dm3
seqlevelsStyle(genome) <- "Ensembl"
genomeSeqInfo <- SeqinfoForBSGenome(genome)
param<-ScanVcfParam(samples = unique(as.character(targets$line)), trimEmpty=TRUE, info = NA, geno = "GT", fixed=c("ALT"))
vcf <- readVcf("./Data/freeze2.vcf.gz", param = param, genome = genomeSeqInfo)
geno <- geno(vcf)$GT
rr <- rowRanges(vcf)

rr$ALT <- unlist(rr$ALT)
filtwidth <- (width(rr) == 1) & (width(rr$ALT) == 1)
rr <- rr[filtwidth]
geno <- geno[filtwidth,]
head(geno)
# filter by allele counts and hets
nhet <- apply(geno, MAR = 1, function(x) {sum( (x == "1/0") | x == "0/1")})
nref <- apply(geno, MAR = 1, function(x) {sum( (x == "0/0") )})
nalt <- apply(geno, MAR = 1, function(x) {sum( (x == "1/1") )})
nmis <- apply(geno, MAR = 1, function(x) {sum( (x == ".") )})
filterallele <- (nref > 5) & (nalt > 5) & ((nhet + nmis) < 20)
table(filterallele)
geno <- geno[filterallele,]
rr <- rr[filterallele]
rm(nhet, nref, nmis, nalt, filtwidth, filterallele)

genoM <- t(apply(geno, MAR = 1, function(x){
  x[x == "."] <- "NA"
  x[x == "0/0"] <- "0"
  x[x == "1/1"] <- "2"
  x[x == "1/0"] <- "NA"
  x[x == "0/1"] <- "NA"
  as.numeric(x)
}))
colnames(genoM) <- colnames(geno)
genoM <- genoM[,match(as.character(targets$line), colnames(genoM))]
colnames(genoM) <- targets$sampleID

subsets <- list("./Data/genoForMatrixeQTL.noS14.noS54.txt" = "All",
                "./Data/genoForMatrixeQTL.noS14.noS54.Naive.txt" = "Naive",
                "./Data/genoForMatrixeQTL.noS14.noS54.Treated.txt" = "Treated")

mapply(x=subsets, y = names(subsets), function(x,y){
  if(x == 'All'){
    ind = 1:ncol(genoM)
  }else{
    ind = which(grepl(x,colnames(genoM)))
  }
  cursub <- cbind(data.frame(id = rownames(genoM)), genoM[,c(ind)])
  write.table(cursub, y, quote=FALSE, row.names = FALSE)
}, SIMPLIFY = F)

snploc <- data.frame(snp = rownames(genoM),
                     chr = seqnames(rr),
                     pos = start(rr))
snploc$chr <- as.character(plyr::revalue(as.character(snploc$chr), c("2L"="1", 
                                                   "2R"="2",
                                                   "3L"="3",
                                                   "3R"="4",
                                                   "4"="5",
                                                   "X"="6")))
write.table(snploc, file = "./Data/snpsloc.txt", quote=FALSE, row.names = FALSE)



#############################################
#####     Covariate file preparation    #####
#############################################
# first write out the vcf file for this study, then convert to gds format
VariantAnnotation::writeVcf(vcf, "./Data/freeze2.thisStudy.vcf.gz")
library(SNPRelate)
genofilepath <- "./Data/freeze2.filtered.eQTLs.gds"
snpgdsVCF2GDS("./Data/freeze2.thisStudy.vcf.gz", genofilepath, method="biallelic.only")
genofile <- snpgdsOpen(genofilepath)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, num.thread = 12)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=12)
save(pca, file = "./Data/pca.RData")
# plot(pca)
firstThreePC <- pca$eigenvect[,1:3]
rownames(firstThreePC) <- pca$sample.id
colnames(firstThreePC) <- paste0("PC", 1:3) 
write.csv(firstThreePC, "./Data/firstThreePC.csv")

targets <- read.table("./Data/design.txt", header=TRUE)
targets <- cbind(targets, firstThreePC[as.character(targets$line), ])

## add wolbachia status
wolb <- openxlsx::read.xlsx("./Data/DGRP2_website/wolbachia.xlsx")
wolb$DGRP.Line[nchar(wolb$DGRP.Line) == 8] <- gsub("__", "-0", wolb$DGRP.Line[nchar(wolb$DGRP.Line) == 8])
wolb$DGRP.Line[nchar(wolb$DGRP.Line) == 9] <- gsub("__", "-", wolb$DGRP.Line[nchar(wolb$DGRP.Line) == 9])
targets$wolbachia <- plyr::mapvalues(wolb$Infection.Status[match(targets$line, wolb$DGRP.Line)], c("n","y"), c(0,1))
targets$condition <- plyr::mapvalues(targets$condition, c("Naive", "Treated"), c(0,1))

covariateFile <- targets[,!colnames(targets) %in% c("line", "sample", "stock", "resistance", "D3.uncorrected.death")]
covariateFile <- covariateFile[,c("sampleID", setdiff(colnames(covariateFile), c("sampleID", "condition")), "condition")]
colnames(covariateFile) <- gsub("sampleID", "id", colnames(covariateFile))

write.table(t(covariateFile), "./Data/covariates.pc.wolb.inv.txt", quote = F, col.names = F, row.names = T, sep = ",")
write.table(t(covariateFile[covariateFile$condition == 0, colnames(covariateFile) != "condition"]), "./Data/covariates.pc.wolb.inv.Naive.txt", quote = F, col.names = F, row.names = T, sep = ",")
write.table(t(covariateFile[covariateFile$condition == 1, colnames(covariateFile) != "condition"]), "./Data/covariates.pc.wolb.inv.Treated.txt", quote = F, col.names = F, row.names = T, sep = ",")

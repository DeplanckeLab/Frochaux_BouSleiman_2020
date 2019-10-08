#######################
### eQTL validation ###
#######################

### useful stuff at the beginnings ####
# Use this script to xtract DNA stretch around each variant to check for TFBS using FIMO

## paths ####
data.path <- "./Data"
fig.path <- "./Figures"
processed.path <- "./Processed"

## libraries ####
library(tictoc)
library(lme4) # for glmer fct, linear mixed model
library(ggplot2)
library(reshape2)
library(cowplot)

## functions #####
get_lmm_pval2 <- function(in.mat) {
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
colPanel.snp <- c("lightblue", "orange")

### END useful stuff at the beginnings ####


############
### Main ###
############

# use to make plot of specific gene expression based on a specific snp

### make plot ####
## Analysis ####
data.t <- readRDS(file.path(processed.path, "Infected_glmm_c_max6_quantile.rds"))
data.n <- readRDS(file.path(processed.path, "Ctrl_glmm_c_max6_quantile.rds"))

col.vect <- c("chrom", "position", "variant", "unique", "gene", "beta", "t.stat", "p.value", "FDR", "nb_NU", "nb_Hom_alt", 
              "nb_Hom_ref", "nb_Het", "lmm.pval", "lmm.eSize", "usable.cross", "adj.lmm.pval", "preBH.signif", "postBH.signif")

# select one of the data set
# data <- data.t[,col.vect]
# data <- data.n[,col.vect]

nrow(data.t[data.t$gene=="FBgn0035461",])
nrow(data.n[data.n$gene=="FBgn0035461",])

# data.t <- data.t[data.t$gene=="FBgn0035461",col.vect]
# data.n <- data.n[data.n$gene=="FBgn0035461",col.vect]

SNP <- "3L_3814085_G/T"
current.treatment <- "Pe"

het.nb <- data.t[data.t$unique==SNP, "nb_Het"]

# 2. create empty matrix for each replica
rep_1 <- as.data.frame(matrix(NA, nrow = het.nb, ncol=4))
rep_2 <- as.data.frame(matrix(NA, nrow = het.nb, ncol=4))
colnames(rep_1) <- c("alt", "ref", "cross", "replica")
colnames(rep_2) <- c("alt", "ref", "cross", "replica")

# 3. get the name of het crosses
sub.line <- data.t[data.t$unique==SNP, as.vector(grep(pattern = "_status", colnames(data.t)))]
# if (het.nb == 1) {
#   sub.line[19] <- "Het"
#   sub.line <- sub.line[,which(sub.line == "Het")]
#   tmp.cr_vect <- gsub("_status", "", colnames(sub.line))
#   tmp.cr_vect <- tmp.cr_vect[-2]
# }
# else {
#   sub.line <- sub.line[,which(sub.line == "Het")]
#   tmp.cr_vect <- gsub("_status", "", colnames(sub.line))
# }

sub.line <- sub.line[,which(sub.line == "Het")]
tmp.cr_vect <- gsub("_status", "", colnames(sub.line))

# 4. fill the matrix with cross and replica
rep_1[,3] <- tmp.cr_vect
rep_2[,3] <- tmp.cr_vect
rep_1[,4] <- rep(paste(current.treatment, "1", sep="_"), het.nb)
rep_2[,4] <- rep(paste(current.treatment, "2", sep="_"), het.nb)


# 5. loop through all the het crosses to fill the matrix,
for (cr in tmp.cr_vect) {
  if (cr == "c20x1" & current.treatment == "Pe") {
    # IMPORTANT, remove c20x1 Pe 2 because sample failed, leave as NA and remove NA later
    rep_1[rep_1$cross == cr, "alt"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "1", "alt", sep="_"), colnames(data.t))]
    rep_1[rep_1$cross == cr, "ref"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "1", "ref", sep="_"), colnames(data.t))]
  }
  else {
    rep_1[rep_1$cross == cr, "alt"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "1", "alt", sep="_"), colnames(data.t))]
    rep_2[rep_2$cross == cr, "alt"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "2", "alt", sep="_"), colnames(data.t))]
    rep_1[rep_1$cross == cr, "ref"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "1", "ref", sep="_"), colnames(data.t))]
    rep_2[rep_2$cross == cr, "ref"] <- data.t[data.t$unique==SNP, grep(pattern = paste(cr, current.treatment, "2", "ref", sep="_"), colnames(data.t))]
  }
}


# 6. aggregate matrix 
tmp.mat <- rbind(rep_1, rep_2)
tmp.mat$total <- tmp.mat$alt + tmp.mat$ref

# remove c20x1 Pe_2 from tmp.mat
if (current.treatment == "Pe") {
  tmp.mat <- tmp.mat[complete.cases(tmp.mat),]
}

# cutoff
cutoff <- max(6, quantile(tmp.mat$total, na.rm = T)[2])
test.mat <- tmp.mat[tmp.mat$total >= cutoff,]

glmm.cutoff <- get_lmm_pval2(test.mat)
glmm.noCutoff <- get_lmm_pval2(tmp.mat)

glmm.cutoff[4]
glmm.noCutoff[4]

tmp.mat$ratio <- tmp.mat$alt / tmp.mat$ref
ggplot(tmp.mat, aes(x=ratio, y=log2(ratio+1))) + geom_violin()

plot.mat <- melt(tmp.mat, measure.vars = c("alt", "ref"))
p1 <- ggplot(plot.mat, aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  scale_fill_manual(values = colPanel.snp) +
  geom_jitter(aes(color=variable), shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values = c("blue", "darkorange3")) +
  xlab("Variant") + ylab("Counts")

save_plot(file.path(fig.path, "ntc_SNP_counts.pdf"), p1)
### END Analysis ####

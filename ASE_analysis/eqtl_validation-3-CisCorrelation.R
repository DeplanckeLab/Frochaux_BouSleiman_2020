#########################
### eQTL validation 3 ###
#########################

### useful stuff at the beginnings ####

## paths ####
# does not change
data.path <- "./Data"
fig.path <- "./Figures"
processed.path <- "./Processed"

## libraries ####
library(lme4) # for glmer fct, linear mixed model
library(tictoc)
library(tidyr)
library(ggplot2)
library(gtools) # for smartbind function
library(reshape2)
library(cowplot)
library(hexbin)
library(lattice)
library(RColorBrewer)

## functions #####

## vectors ####

## colPanels ####

# color dataframe for GenePlot fct
col.df <- as.data.frame(matrix(NA, nrow=4, ncol=2))
colnames(col.df) <- c("status", "col")
col.df$status <- c("Het", "Hom_alt", "Hom_ref", "NU")
col.df$col <- c("green", "black", "grey", "red")
colPanel.infection <- c("#c0c0c0", "#cc6633")

### END useful stuff at the beginnings ####


################
### Analysis ###
################

### Load data ####
data.t <- readRDS(file.path(processed.path, "Infected_glmm_c_max6_quantile.rds"))

data.n <- readRDS(file.path(processed.path, "Ctrl_glmm_c_max6_quantile.rds"))

### Histogram of number of Heterozygote cross per eQTL and nb of validated eqtl ####

#jpeg(file.path(figures.path, "data_Het_cnt.jpg"), 800, 600)
Het.table <- table(data.t[,c("nb_Het")])
Het.table <- table(data.n[,c("nb_Het")])
ylim <- c(0, 1.1*max(Het.table))
xx <- barplot(Het.table, main = "frequency of nb of Het cross", ylim = ylim,
              xlab="Number of Heterozygote cross",
              ylab="Number of eQTLs")
text(x = xx, y = Het.table, label = as.character(Het.table), 
     pos=3, cex=0.7, col="red")
#dev.off()

## combined plot with validated eqtl highlighted ####

Het.table.t <- table(data.t[,c("nb_Het")])
Het.table.n <- table(data.n[,c("nb_Het")])

table.combined1 <- as.data.frame(rbind(Het.table.n, Het.table.t))
row.names(table.combined1) <- c("Ctrl", "Infected")

# add line with number of validated eqtl
Het.table.val.t <- table(data.t[data.t$postBH.signif==1,"nb_Het"])
Het.table.val.n <- table(data.n[data.n$postBH.signif==1,"nb_Het"])

table.combined2 <- as.data.frame(rbind(Het.table.val.n, Het.table.val.t))
row.names(table.combined2) <- c("Ctrl", "Infected")



combined1.melt <- melt(as.matrix(table.combined1))
colnames(combined1.melt) <- c("treatement", "variable", "all")
combined2.melt <- melt(as.matrix(table.combined2))
colnames(combined2.melt) <- c("treatement", "variable", "validated")

combined.melt <- merge(combined1.melt, combined2.melt, all = T)
combined.melt[is.na(combined.melt)] <- 0

P1 <- ggplot(data=combined.melt) +
  geom_bar(aes(x=variable, y=all, fill=treatement), stat="identity", position="dodge", alpha=0.5) +
  geom_bar(aes(x=variable, y=validated, fill=treatement), stat="identity", position="dodge") +
  scale_fill_manual(values=colPanel.infection) +
  labs(title="Distribution of eQTL based on number of Heterozygote cross",
       x="Number of Heterozygote cross", y="Number of eQTL", fill="Treatment") +
  theme(legend.position="none")

P1
P1 + theme_cowplot()

save_plot(file.path(fig.path, "eqtl_distribution_bar_with_valid.pdf"), P1 + theme_cowplot(), base_aspect_ratio = 1.5)

# ggsave("eqtl_distribution_bar_with_valid.pdf", P1, device="pdf", 
#        path = "C:\\Users/Michael/Documents/Deplancke/Figures/Fig 3/")

## clean up
rm(Het.table, Het.table.n, Het.table.t, Het.table.val.n, Het.table.val.t, combined.melt, combined1.melt, combined2.melt, P1)

### END ####

### nb eqtl per gene ####
eqtl.table.t <- table(data.t$gene)
eqtl.table.n <- table(data.n$gene)

eqtl.table.combined <- as.data.frame(rbind(eqtl.table.n, eqtl.table.t))
row.names(eqtl.table.combined) <- c("naive", "treated")

eqtl.combined.melt <- melt(as.matrix(eqtl.table.combined))

P2 <- ggplot(data=eqtl.combined.melt, aes(x=value, fill=Var1)) +
  geom_histogram(binwidth=1, position="dodge") +
  scale_fill_manual(values=colPanel.infection) +
  labs(title="Distribution of number of eQTL per gene",
       x="Number of eQTL", y="count",
       tag="c", fill="Treatment")

P2 + theme_cowplot()

ggsave("nb_eqtl_perGene_bar.pdf", P2, device="pdf", path = fig.path)
### END ####


### P-value comparison ####
plot.mat <- data.t[,c("FDR", "adj.lmm.pval")]; set <- "Infected"
plot.mat <- data.n[,c("FDR", "adj.lmm.pval")]; set <- "Ctrl"

print(set)
round(cor(plot.mat$FDR, plot.mat$adj.lmm.pval, use = "complete.obs"), 5)

plot.mat2 <- -log10(plot.mat)
plot.mat2 <- na.omit(plot.mat2)
plot.mat2 <- plot.mat2[is.finite(plot.mat2[,2]),]
rf <- colorRampPalette(rev(brewer.pal(9,"BuPu")))

pdf(file.path(fig.path, paste(set, "pvalue_comparison.pdf", sep="_")))
hexbinplot(adj.lmm.pval ~ FDR, data=plot.mat2, 
           aspect="1", xbins=70, 
           cex.labels=1.0, 
           cex.title=1.0, 
           colramp=rf,
           panel=function(x, y, ...)
           {
             panel.hexbinplot(x, y, ...)
             panel.abline(v=-log10(0.05), h=-log10(0.05), col="red", lwd=1, lty=5)
           }
)
dev.off()
### END ####

### effect size correlation ####

## 1. subset data for relevant columns ####
plot.mat <- data.t[,c("beta", "lmm.eSize")]; set <- "Infected"
plot.mat <- data.n[,c("beta", "lmm.eSize")]; set <- "Ctrl"

plot.mat <- na.omit(plot.mat)
rf <- colorRampPalette(rev(brewer.pal(9,"BuPu")))

## Plot 1: hexbin plot with all eQTLs ####
hbin <- hexbin(plot.mat, xbins=70)
pdf(file.path(fig.path, paste(set, "eSize_comp.pdf", sep="_")))
plot(hbin,colramp=rf,border=gray(.75), main=paste("Comparison between predicted effect \nand measured effect - all samples", set, sep=" "),
     xlab="Predicted effect", ylab="Measured effect")
dev.off()

## Plot 2. only cis-eQTLs ####
plot.mat.s <- data.t[data.t$postBH.signif == 1,c("beta", "lmm.eSize")]; set <- "Infected"
plot.mat.s <- data.n[data.n$postBH.signif == 1,c("beta", "lmm.eSize")]; set <- "Ctrl"

plot.mat.s <- na.omit(plot.mat.s)

hbin2 <- hexbin(plot.mat.s, xbins=70)
pdf(file.path(fig.path, paste(set, "cis_esize_comp.pdf", sep="_")))
plot(hbin2,colramp=rf,border=gray(.75), main=paste("Comparison between predicted effect \nand measured effect - significant samples", set, sep=" "),
     xlab="Predicted effect", ylab="Measured effect")
dev.off()

## Plot 3: stripchart + boxplot based on beta value ####
plot.mat.s$b.stat <- NA
plot.mat$b.stat <- NA

plot.mat.s[plot.mat.s$beta < 0, "b.stat"] <- "neg"
plot.mat.s[plot.mat.s$beta > 0, "b.stat"] <- "pos"

p5 <- ggplot(plot.mat.s, aes(x=b.stat, y=lmm.eSize, color=b.stat)) +
  geom_boxplot(color="black", notch = TRUE)+
  geom_jitter(position=position_jitter(0.3)) +
  scale_color_manual(values=c("gold", "black")) +
  labs(main=paste("Comparison between predicted effect \nand measured effect - significant samples", set, sep=" "), 
       x="beta", y="measured effect size")
p5 + theme_cowplot()

save_plot(file.path(fig.path, paste(set, "uglyBoxplot.pdf", sep="_")), p5 + theme_cowplot())
### END ####
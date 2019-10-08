#######################
### eQTL validation ###
#######################

### useful stuff at the beginnings ####

## paths ####
data.infected <- "./Data/Treated_parsed_CrossInfo_CntData.txt"
data.ctrl <- "./Data/Naive_parsed_CrossInfo_CntData.txt"

fig.path <- "./Figures"
processed.path <- "./Processed"

## libraries ####
library(lme4) # for glmer fct, linear mixed model
library(tictoc)
library(reshape2)
library(ggplot2)
library(cowplot)

## functions ####
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
## vectors with informations ####
cross.vect <- c("c1x2", "c2x3", "c3x4", "c4x5", "c5x6", "c6x7", "c7x8", "c8x9",
                "c9x10", "c10x11", "c11x12", "c12x13", "c13x14", "c14x15",
                "c17x18", "c18x19", "c19x20", "c20x1")
treat.vect <- c("Pe_1", "Pe_2")
treat.vect <- c("UC_1", "UC_2")
#type.vect <- c("alt", "ref")
type.vect <- c("alt", "ref", "total")
all.treat <- c("Pe_1", "Pe_2", "UC_1", "UC_2")

### END useful stuff at the beginnings ####

################
### Analysis ###
################

### GLMM ####

## 1. load and pre-process data ####
# Infected
data <- read.delim(data.infected, header=T); set <- "Infected"

# Ctrl
data <- read.delim(data.ctrl, header=T); set <- "Ctrl"

# Interaction
#data <- read.delim(file.path(data.path, "Naive_parsed_CrossInfo_CntData.txt"), header=T); set <- "Interaction"

# remove columns with c15x16 and c16x17 because the cross did not work
# vector of colnames
c16.vect <- grep(pattern = "16", colnames(data))
data[,c16.vect] <- NULL

data <- data[complete.cases(data),]

# Add stats on number of cross type per eQTL
data$nb_NU <- rowSums(data == "NU", na.rm=T)
data$nb_Hom_alt <- rowSums(data == "Hom_alt", na.rm=T)
data$nb_Hom_ref <- rowSums(data == "Hom_ref", na.rm=T)
data$nb_Het <- rowSums(data == "Het", na.rm=T)

# add unique identifier to each eqtl
data$unique <- paste(data$chrom, data$position, data$variant, sep="_")

## 2. fitting model ####

#print(paste("cutoff =", cutoff, sep=" "))
for (line.it in 1:nrow(data)) {
  
  if (line.it == 1) {
    if (set == "Infected") {current.treatment <- "Pe"}
    if (set == "Ctrl") {current.treatment <- "UC"}
    
    print(paste("in use set:", current.treatment, sep=" "))
    tic("fitting model")
    
    data$lmm.pval <- NA
    data$lmm.eSize <- NA
    data$usable.cross <- NA
    data$zstat <- NA
    
    err.output <- data[,1:4]
    err.output$nb_het <- data$nb_Het
    err.output$tot.cross_t <- NA
    err.output$pass_cutoff <- NA
    err.output$cr_pass_cutoff <- NA
    err.output$error <- NA
    err.output$cutoff <- NA
  }
  
  #for (line.it in 1:10) {
  # loop information
  if (line.it %% 100 == 0) {
    print(line.it)
  }
  
  # 1. get the number of heterozygotes cross
  het.nb <- data[line.it, "nb_Het"]
  
  # if there is no heterozygote cross, stop the loop
  if (het.nb == 0) {
    err.output[line.it, "error"] <- "nb_Het"
    next
  }
  
  # 2. create empty matrix for each replica
  rep_1 <- as.data.frame(matrix(NA, nrow = het.nb, ncol=4))
  rep_2 <- as.data.frame(matrix(NA, nrow = het.nb, ncol=4))
  colnames(rep_1) <- c("alt", "ref", "cross", "replica")
  colnames(rep_2) <- c("alt", "ref", "cross", "replica")
  
  # 3. get the name of het crosses
  sub.line <- data[line.it, grep(pattern = "_status", colnames(data))]
  if (het.nb == 1) {
    sub.line[19] <- "Het"
    sub.line <- sub.line[,which(sub.line == "Het")]
    tmp.cr_vect <- gsub("_status", "", colnames(sub.line))
    tmp.cr_vect <- tmp.cr_vect[-2]
  }
  else {
    sub.line <- sub.line[,which(sub.line == "Het")]
    tmp.cr_vect <- gsub("_status", "", colnames(sub.line))
  }
  
  # 4. fill the matrix with cross and replica
  rep_1[,3] <- tmp.cr_vect
  rep_2[,3] <- tmp.cr_vect
  rep_1[,4] <- rep(paste(current.treatment, "1", sep="_"), het.nb)
  rep_2[,4] <- rep(paste(current.treatment, "2", sep="_"), het.nb)
  
  
  # 5. loop through all the het crosses to fill the matrix,
  for (cr in tmp.cr_vect) {
    if (cr == "c20x1" & current.treatment == "Pe") {
      # IMPORTANT, remove c20x1 Pe 2 because sample failed, leave as NA and remove NA later
      rep_1[rep_1$cross == cr, "alt"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "1", "alt", sep="_"), colnames(data))]
      rep_1[rep_1$cross == cr, "ref"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "1", "ref", sep="_"), colnames(data))]
    }
    else {
      rep_1[rep_1$cross == cr, "alt"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "1", "alt", sep="_"), colnames(data))]
      rep_2[rep_2$cross == cr, "alt"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "2", "alt", sep="_"), colnames(data))]
      rep_1[rep_1$cross == cr, "ref"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "1", "ref", sep="_"), colnames(data))]
      rep_2[rep_2$cross == cr, "ref"] <- data[line.it, grep(pattern = paste(cr, current.treatment, "2", "ref", sep="_"), colnames(data))]
    }
  }
  
  
  # 6. aggregate matrix 
  tmp.mat <- rbind(rep_1, rep_2)
  tmp.mat$total <- tmp.mat$alt + tmp.mat$ref
  
  # remove c20x1 Pe_2 from tmp.mat
  if (current.treatment == "Pe") {
    tmp.mat <- tmp.mat[complete.cases(tmp.mat),]
  }
  
  err.output[line.it, "tot.cross_t"] <- nrow(tmp.mat)
  
  # 7. apply simple linear mixed model, get p-value, store it in df
  # cutoff of the number of reads
  # there is different way to calculate cutoff, arbitrary cutoff or data-based cutoff
  # arbitrary
  #cutoff <- c_off
  
  cutoff <- max(6, quantile(tmp.mat$total, na.rm = T)[2])
  err.output[line.it, "cutoff"] <- cutoff
  
  test.mat <- tmp.mat[tmp.mat$total >= cutoff,]
  
  err.output[line.it, "pass_cutoff"] <- nrow(test.mat)
  err.output[line.it, "cr_pass_cutoff"] <- length(unique(test.mat$cross))
  
  data[line.it,"usable.cross"] <- nrow(test.mat)
  
  # 8. test if matrix usable
  if (nrow(test.mat) <= 1 ) { 
    err.output[line.it, "error"] <- "sample_cutoff"
    # print(tmp.mat)
    # print(test.mat)
    # readline(prompt="Press [enter] to continue")
    next 
  }
  if (length(unique(test.mat$cross)) <= 1 ) { 
    err.output[line.it, "error"] <- "cross_cutoff"
    next 
  }
  
  # 9. test
  res.out <- get_lmm_pval2(test.mat)
  if (is.na(res.out)) {
    err.output[line.it, "error"] <- "response_constant"
  }
  else {
    res.out.pval <- res.out[4]
    res.out.zstat <- res.out[3]
    res.out.eSize <- res.out[1]
    # p value
    data[line.it, "lmm.pval"] <- res.out.pval
    # z stat
    data[line.it, "zstat"] <- res.out.zstat
    # effect size
    eSize <- exp(res.out.eSize)/(1+exp(res.out.eSize))
    data[line.it, "lmm.eSize"] <-eSize
    # output ok in error
    err.output[line.it, "error"] <- "ok"
  }
}

toc()

## 3. adjusting p values ####
pdf(file.path(fig.path, paste(set, "Pvalue_Hist_PreBH.pdf", sep="_")))
hist(data$lmm.pval, breaks = 100, xlab="p-value", main=paste("Histogram of P-values", set, sep=" - "))
abline(v=0.05, col="red", ylim=c(0,600), lty=2)
text(x = 0.05, y=300, labels = "0.05", col="red", pos=4)
dev.off()

data$adj.lmm.pval <- p.adjust(data$lmm.pval, method="BH")

pdf(file.path(fig.path, paste(set, "Pvalue_Hist_PostBH.pdf", sep="_")))
hist(data$adj.lmm.pval, breaks = 100, xlab="p-value", main=paste("Histogram of P-values after correction", set, sep=" - "))
abline(v=0.05, col="red", ylim=c(0,600), lty=2)
text(x = 0.05, y=300, labels = "0.05", col="red", pos=4)
dev.off()

data$preBH.signif <- NA
data$postBH.signif <- NA

preBH.sig.nb <- 0
postBH.sig.nb <- 0

for (lit in 1:nrow(data)) {
  preBH <- data[lit, "lmm.pval"]
  postBH <- data[lit, "adj.lmm.pval"]
  if (!is.na(preBH)) {
    if (preBH < 0.05) {
      preBH.sig.nb <- preBH.sig.nb + 1
      data[lit, "preBH.signif"] <- 1
    }
    else {
      data[lit, "preBH.signif"] <- 0
    }
  }
  if (!is.na(postBH)) {
    if (postBH < 0.05) {
      postBH.sig.nb <- postBH.sig.nb + 1
      data[lit, "postBH.signif"] <- 1
    }
    else {
      data[lit, "postBH.signif"] <- 0
    }
  }
}


## 4. save output ####

# new
saveRDS(data, file = file.path(processed.path, paste(set, "glmm_c_max6_quantile.rds", sep="_")))
saveRDS(err.output, file = file.path(processed.path, paste(set, "glmm_c_max6_quantile_Error.rds", sep="_")))
### END GLMM ####

### error output analysis ####
## data loading and processing ####
error.n <- readRDS(file.path(processed.path, "Ctrl_glmm_c_max6_quantile_Error.rds"))
error.t <- readRDS(file.path(processed.path, "Infected_glmm_c_max6_quantile_Error.rds"))

# assign one error data set to error for further processing
err.output <- error.n; set="Ctrl"
err.output <- error.t; set="Infected"

significant.n <- 1250
significant.t <- 1301

postBH.sig.nb <- significant.n
postBH.sig.nb <- significant.t

# assess why eqtl cannot be measured
b <- table(err.output$error)
slices <- as.vector(b)
slices <- c(slices, postBH.sig.nb)
slices[3] <- slices[3] - postBH.sig.nb
pct <- round(slices/sum(slices)*100, 1)
lbls <- names(b)
lbls <- c(lbls, "significant")
lbls2 <- paste(lbls, pct)
lbls3 <- paste(lbls2,"%",sep="")

## Plot 1: standard pie chart ####
pie.palette <- c("darkorchid", "gold2", "greenyellow", "gray40", "firebrick3", "darkgreen")
pdf(file.path(fig.path, paste(set, "error_pie_chart.pdf", sep="_")))
pie(slices, labels = lbls3, col=pie.palette, 
    main=paste("Errors - ", set, sep=""))
title(sub = "Cutoff - max(6, 25th quantile)")
dev.off()

## Plot 2: ggplot bar chart ####
# combine slices and labels vectors
err.df <- melt(data.frame(lbls3,slices))
err.df$errors <- lbls
err.df <- err.df[order(err.df$value, na.last = T, decreasing = F),]
err.df$order <- c(1:nrow(err.df))
err.barplot <- ggplot(err.df, aes(x=order, y=value, fill=lbls3))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = paste(set, "errors", sep=" "), fill="errors", y = "Count", x="")

err.barplot + theme_cowplot()
### END Error ####

### Repartition of eQTLs depending on nb of het ####
## PLot 3: barplot of repartition of tested and validated eQTLs ####
# get stats
summary.err.n <- table(error.n$error)
summary.err.t <- table(error.t$error)

error.df <- as.data.frame(rbind(summary.err.n, summary.err.t))
row.names(error.df) <- c("control", "infected")

# add line with number of validated eqtl
combined1.melt <- melt(as.matrix(error.df))
colnames(combined1.melt) <- c("treatement", "variable", "all")
significant.n <- 1250
significant.t <- 1301
combined1.melt$signif <- 0
combined1.melt[combined1.melt$treatement=="control" & combined1.melt$variable=="ok", "signif"] <- significant.n
combined1.melt[combined1.melt$treatement=="infected" & combined1.melt$variable=="ok", "signif"] <- significant.t

simple.df <- combined1.melt
simple.df[11,] <- c("control", "not tested", 423+70+1+2518, 0)
simple.df[12,] <- c("infected", "not tested", 428+72+5+2366, 0)
simple.df <- simple.df[c(5,6,11,12),]
row.names(simple.df) <- c()
simple.df[,2] <- as.factor(c("ok", "ok", "NT", "NT"))
simple.df[,3] <- as.numeric(simple.df[,3])
simple.df[,4] <- as.numeric(simple.df[,4])

saveRDS(simple.df, file=file.path(processed.path, "simple.rds"))
simple.df <- readRDS(file.path(processed.path, "simple.rds"))
simple.df

p3 <- ggplot(data=simple.df, aes(x=variable)) +
  geom_bar(aes(y=all, fill=variable), stat="identity", position="dodge", alpha=0.5) +
  geom_bar(aes(y=signif, fill=variable), stat="identity", position="dodge") +
  facet_wrap(~treatement) +
  scale_fill_manual(values=c("red4", "darkgreen")) +
  labs(title="Results",
       x="tested", y="Number of eQTL") +
  theme(legend.position="none")

p3 + theme_cowplot()
p3 + theme_cowplot() + theme(legend.position="none")

save_plot(file.path(fig.path, "Results_Cis.pdf"), p3 + theme_cowplot(), 1.2)
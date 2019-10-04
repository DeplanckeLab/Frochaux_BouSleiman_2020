##############
### Mitomi ###
##############

### libraries ####
library(ggplot2)
library(cowplot)

### data ####
broad.data <- "./Data/broadMitomi_rep_param.txt"
DaRelSage.data <- "./Data/Mitomi_Rel_Da_Sage.txt"

### Broad ####
#brd.p <- read.delim(file.path(path, "broadMitomi_rep_param.txt"), header=T)
brd.p <- read.delim(broad.data, header=T)

## Variables ####
attach(brd.p)
p1.alt <- brd.p[replica==1 & allele=="alt", "p1"]
kd1.alt <- brd.p[replica==1 & allele=="alt", "kd"]
p2.alt <- brd.p[replica==2 & allele=="alt", "p1"]
kd2.alt <- brd.p[replica==2 & allele=="alt", "kd"]
p3.alt <- brd.p[replica==3 & allele=="alt", "p1"]
kd3.alt <- brd.p[replica==3 & allele=="alt", "kd"]
p1.ref <- brd.p[replica==1 & allele=="ref", "p1"]
kd1.ref <- brd.p[replica==1 & allele=="ref", "kd"]
p2.ref <- brd.p[replica==2 & allele=="ref", "p1"]
kd2.ref <- brd.p[replica==2 & allele=="ref", "kd"]
p3.ref <- brd.p[replica==3 & allele=="ref", "p1"]
kd3.ref <- brd.p[replica==3 & allele=="ref", "kd"]
p.mean.alt <- ((p1.alt + p2.alt + p3.alt)/3)
kd.mean.alt <- ((kd1.alt + kd2.alt + kd3.alt)/3)
p.mean.ref <- ((p1.ref + p2.ref + p3.ref)/3)
kd.mean.ref <- ((kd1.ref + kd2.ref + kd3.ref)/3)
detach(brd.p)

## Line plot ####
size <- 0.5
br.plot <- ggplot(data.frame(x=c(0, 41000)), aes(x)) + 
  stat_function(fun=function(x) (x*p1.ref)/(x+kd1.ref), colour="orange", size=size, linetype = "dotdash") +
  stat_function(fun=function(x) (x*p1.alt)/(x+kd1.alt), colour="lightblue", size=size, linetype = "twodash") +
  stat_function(fun=function(x) (x*p2.ref)/(x+kd2.ref), colour="orange", size=size, linetype = "dotdash") +
  stat_function(fun=function(x) (x*p2.alt)/(x+kd2.alt), colour="lightblue", size=size, linetype = "twodash") +
  stat_function(fun=function(x) (x*p3.ref)/(x+kd3.ref), colour="orange", size=size, linetype = "dotdash") +
  stat_function(fun=function(x) (x*p3.alt)/(x+kd3.alt), colour="lightblue", size=size, linetype = "twodash") +
  stat_function(fun=function(x) (x*p.mean.ref)/(x+kd.mean.ref), colour="orange", size=1) +
  stat_function(fun=function(x) (x*p.mean.alt)/(x+kd.mean.alt), colour="lightblue", size=1) +
  ggtitle(label = "Broad") + xlab(label = "Concentration") + ylab(label = "RFU") +
  theme_cowplot()
br.plot

save_plot("./Figures/Broad_differential_binding.pdf", br.plot)

## bosplot ####
kd.plot <- ggplot(brd.p, aes(x=allele, y=kd, fill=allele)) + geom_boxplot() + scale_fill_manual(values=c("lightblue", "orange")) +
  ggtitle("Dissociation constant - Broad") + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour="black", fill=NA)) +
  theme_cowplot() +
  theme(legend.position = "none")
kd.plot

save_plot("./Figures/Broad_kd_boxplot.pdf", kd.plot)

## stat ####
alt.val <- as.vector(brd.p[brd.p$allele=="alt", "p1"])
ref.val <- as.vector(brd.p[brd.p$allele=="ref", "p1"])
t.test(alt.val, ref.val, alternative = "two.sided", var.equal = FALSE)

### Other proteins ####
prots.dat <- read.table(DaRelSage.data, header=T)

## Variables ####
# Sage 
# alt
p1.sage.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="sage","p1"]
kd.sage.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="sage","kd"]

# ref
p1.sage.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="sage","p1"]
kd.sage.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="sage","kd"]

# Relish 
# alt
p1.rel.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="relish","p1"]
kd.rel.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="relish","kd"]

# ref
p1.rel.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="relish","p1"]
kd.rel.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="relish","kd"]

# Daughterless
# alt
p1.da.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="daughterless","p1"]
kd.da.alt <- prots.dat[prots.dat$allele=="alt" & prots.dat$protein=="daughterless","kd"]

# ref
p1.da.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="daughterless","p1"]
kd.da.ref <- prots.dat[prots.dat$allele=="ref" & prots.dat$protein=="daughterless","kd"]

## Plots ####

da.plot <- ggplot(data.frame(x=c(0, 40000)), aes(x)) + 
  stat_function(fun=function(x) (x*p1.da.ref)/(x+kd.da.ref), colour="orange", size=1) +
  stat_function(fun=function(x) (x*p1.da.alt)/(x+kd.da.alt), colour="lightblue", size=1) +
  ggtitle(label = "Daughterless") + xlab(label = "Concentration") + ylab(label = "RFU") +
  theme_cowplot()

rel.plot <- ggplot(data.frame(x=c(0, 40000)), aes(x)) + 
  stat_function(fun=function(x) (x*p1.rel.ref)/(x+kd.rel.ref), colour="orange", size=1) +
  stat_function(fun=function(x) (x*p1.rel.alt)/(x+kd.rel.alt), colour="lightblue", size=1) +
  ggtitle(label = "Relish") + xlab(label = "Concentration") + ylab(label = "RFU") +
  theme_cowplot()

sage.plot <- ggplot(data.frame(x=c(0, 40000)), aes(x)) + 
  stat_function(fun=function(x) (x*p1.sage.ref)/(x+kd.sage.ref), colour="orange", size=1) +
  stat_function(fun=function(x) (x*p1.sage.alt)/(x+kd.sage.alt), colour="lightblue", size=1) +
  ggtitle(label = "Sage") + xlab(label = "Concentration") + ylab(label = "RFU") +
  theme_cowplot()

save_plot("./Figures/da_plot.pdf", da.plot)
save_plot("./Figures/rel_plot.pdf", rel.plot)
save_plot("./Figures/sage_plot.pdf", sage.plot)

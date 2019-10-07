### association sharing analysis >> power implications

mateQTLfiles <- list.files("./Data/MatrixEQTL_noSigThresh/", pattern = "^MatrixeQTL(.)+\\.cis$", full.names = T) ## need to re-run matrix eqtl and remove the cis.pvalue threshold

keep <- grepl("withCovariates", mateQTLfiles) & !grepl("Interaction|FULL", mateQTLfiles)
mateQTLfiles <- mateQTLfiles[keep]
names(mateQTLfiles) <- ifelse(grepl("Naive", mateQTLfiles), "Control", "Infected")
library(data.table)

res <- lapply(mateQTLfiles, fread)

res <- mapply(x=res, y=names(res), function(x,y){
  assoc <- paste0(x$SNP, "&", x$gene)
  # colnames(x) <- paste0(y, colnames(x))
  x$assoc <- assoc
  x
}, SIMPLIFY = F)
\

allAssoc <- unique(res$Control$assoc, res$Infected$assoc)

allres <- do.call(cbind, lapply(res, function(x){
  out <- x[match(allAssoc, x$assoc),]
  out$sig <- out$FDR < 0.05
  out
}))


head(allres)

allres$sigBoth <- allres$Control.sig & allres$Infected.sig
allres$sigOne <- (allres$Control.sig | allres$Infected.sig) & !allres$sigBoth
allres$sigAny <- (allres$Control.sig | allres$Infected.sig)
allres$sigGroup <- "Not Sig"
allres$sigGroup[allres$sigBoth] <- "Shared"
allres$sigGroup[allres$sigOne & allres$Control.sig] <- "Sig Control"
allres$sigGroup[allres$sigOne & allres$Infected.sig] <- "Sig Infected"
allres$sigGroup <- factor(allres$sigGroup, levels = c("Shared", "Sig Control", "Sig Infected","Not Sig"))

table(allres$sigBoth)
table(allres$sigOne)
table(allres$Infected.sig)
table(allres$Control.sig)
table(allres$sigBoth, allres$sigOne)
### need to add maf

geno <- fread("./Data/genoForMatrixeQTL.noS14.noS54.Naive.txt")
nhet <- rowSums(geno == 2, na.rm = T)
nhomo <-  rowSums(geno == 0, na.rm = T)
nna <- rowSums(is.na(geno), na.rm = T)
nt <- 38-nna
mac <- mapply(x=nhet, y=nhomo, function(x,y){min(x,y)})
min(mac)
maf <- mac/nt
min(maf)
names(maf) <- geno$id

## append maf to the allres object
allres$maf <- maf[as.character(allres$Control.SNP)]
allres$mafGroup <- cut_width(allres$maf,width = 0.05, closed = 'left', boundary = 0.05)
library(ggplot2)
library(cowplot)
library(ggrepel)

tb <- table(allres$sigOne, allres$mafGroup)
cstb <- colSums(tb)
tbp <- t(apply(tb, MAR = 1, function(x){100*x/cstb}))


tbb <- table(allres$sigBoth, allres$mafGroup)
cstbb <- colSums(tbb)
tbbp <- t(apply(tbb, MAR = 1, function(x){100*x/cstbb}))


tbba <- table(allres$sigAny, allres$mafGroup)
cstbba <- colSums(tbba)
tbbpa <- t(apply(tbba, MAR = 1, function(x){100*x/cstbba}))


tbdt <- data.frame(mafGroup = factor(colnames(tb), levels = levels(allres$mafGroup )),
                   numberSigOne = tb['TRUE',],
                   percentageSigOne = tbp['TRUE',],
                   numberSigBoth = tbb['TRUE',],
                   percentageSigBoth = tbbp['TRUE',],
                   numberSigAny = tbba['TRUE',],
                   percentageSigAny = tbbpa['TRUE',])

ggplot(tbdt, aes(x = percentageSigOne, y = percentageSigBoth, color = mafGroup)) + 
  geom_text_repel(aes(label = mafGroup)) +
  geom_point() +
  geom_line(group = 1) +
  theme_cowplot() +
  scale_color_viridis_d() +
  scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA))


rp <- ggplot(tbdt, aes(x = mafGroup,   y = numberSigOne/numberSigBoth)) + 
  geom_point() +
  geom_line(group = 1) +
  theme_cowplot() +
  scale_color_viridis_d() +
  # scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("MAF group") +
  ylab("Sig. in one\n/ Sig. in both (number)")

rp_percentage <- ggplot(tbdt, aes(x = mafGroup,   y = percentageSigOne/percentageSigBoth)) + 
  geom_point() +
  geom_line(group = 1) +
  theme_cowplot() +
  scale_color_viridis_d() +
  # scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("MAF group") +
  ylab("Sig. in one\n/ Sig. in both (%)")

# arrange in another format
tbdtmelt <- melt(tbdt)

lp <- ggplot(tbdtmelt[tbdtmelt$variable %in% c("percentageSigOne", "percentageSigBoth", "percentageSigAny"),], aes(x = mafGroup,   y = value, color = variable, group = variable)) + 
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_y_continuous(limits = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("MAF group") +
  ylab("Percentage of all\ntested associations")

l <- ggplot(tbdtmelt[tbdtmelt$variable %in% c("numberSigOne", "numberSigBoth", "numberSigAny"),], aes(x = mafGroup,   y = value, color = variable, group = variable)) + 
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_y_continuous(limits = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("MAF group") +
  ylab("Number of associations")


plot_grid(l, rp, lp, rp_percentage, nrow = 2, rel_widths = c(2,1,2,1),
          labels = c("a", "c", "b", "d")) + 
  ggsave("./Plots/eQTL_sharing_power.pdf", width = 8, height = 8) +
  ggsave("./Plots/eQTL_sharing_power.png", width = 8, height = 8)

# ggplot(allres, aes(x = -log10(`Control.p-value`), y = -log10(`Infected.p-value`), color = sigGroup)) +
#   geom_hex(bins = 50) +
#   facet_wrap(~mafGroup) +
#   theme_cowplot()
# 
# 
# ggplot(allres[allres$sigAny,], aes(x = -log10(`Control.p-value`), y = -log10(`Infected.p-value`), color = sigGroup)) +
#   geom_point() +
#   facet_wrap(~mafGroup) +
#   theme_cowplot()
# 
# ggplot(allres[allres$sigAny,], aes(x = -log10(`Control.p-value`), y = -log10(`Infected.p-value`), color = sigGroup)) +
#   geom_point() +
#   facet_wrap(~mafGroup) +
#   theme_cowplot()
# 
# ggplot(allres[allres$sigAny,], aes(x = `Control.beta`, y = `Infected.beta`, color = sigGroup)) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_point() +
#   # geom_hex(bins = 50) +
#   # facet_wrap(~sigGroup) +
#   theme_cowplot()
# 
# 
# ggplot(allres[allres$sigAny,], aes(x = `Control.beta`, y = -log10(`Control.p-value`), color = sigGroup)) +
#   geom_point() +
#   theme_cowplot()
# ggplot(allres[allres$sigAny,], aes(x = `Infected.beta`, y = -log10(`Infected.p-value`), color = sigGroup)) +
#   geom_point() +
#   theme_cowplot()


# icistarget plots
comparison <- read.table("./Data/icistarget_analyses/i-cisTarget-1Naive2TreatedComparison.txt", sep = "\t", header = T, stringsAsFactors = F)
comparison$Control.NES <- comparison$NES.for.results1
comparison$Infected.NES <- comparison$NES.for.results2

comparison$NES.deviation <- apply(comparison[,c("Control.NES", "Infected.NES")], MAR = 1, FUN = function(x){
  (max(x) - min(x))
})

comparison$InfectedVsControl <- comparison$Infected.NES - comparison$Control.NES
comparison$Source <- gsub("_(.)+$","", comparison$Feature.ID)
comparison$Label <- paste0(comparison$Feature.ID,"\n", substr(comparison$Feature.description,1,15), "...")
comparison$LabelSimple <- comparison$Feature.annotation.for.results1


## since difference sources have different annotations/formats, create a new label for plot depending on source
comparison <- split(comparison, comparison$Source)

comparison <- lapply(comparison, function(x){
  if(x$Source[1] == "bergman"){
    x$LabelBySource <- paste("Bergman", x$Feature.description, sep = "\n")
  }else if(x$Source[1] == "cisbp"){
    f <- x$Feature.annotation.for.results1
    f[f==""] <- gsub("\\[(.)+", "", x$Feature.description)[f==""]
    x$LabelBySource <- paste("CIS-BP", f, sep = "\n")
  }else if(x$Source[1] == "elemento"){
    x$LabelBySource <- gsub("elemento__", "Elemento\n", x$Feature.ID)
  }else if(x$Source[1] == "factorbook"){
    x$LabelBySource <- gsub("factorbook__", "Factorbook\n", x$Feature.ID)
  }else if(x$Source[1] == "flyfactorsurvey"){
    x$LabelBySource <- gsub("_(.)+$", "", gsub("flyfactorsurvey__", "FlyFactorSurvey\n", x$Feature.ID))
  }else if(x$Source[1] == "hdpi"){
    x$LabelBySource <- gsub("hdpi__", "HDPI\n", x$Feature.ID)
  }else if(x$Source[1] == "hocomoco"){
    x$LabelBySource <- gsub("hocomoco__", "HOCOMOCO\n", gsub("\\.(.)+$", "", x$Feature.ID))
    x$LabelBySource <- gsub("_HUMAN", " (Human)", x$LabelBySource)
    x$LabelBySource <- gsub("_MOUSE", " (Mouse)", x$LabelBySource)
  }else if(x$Source[1] == "homer"){
    x$LabelBySource <- gsub("homer__", "Homer\n", gsub("__(.)+_", "__", x$Feature.ID))
  }else if(x$Source[1] == "idmmpmm"){
    x$LabelBySource <- gsub("idmmpmm__", "iDMMPMM\n", x$Feature.ID)
  }else if(x$Source[1] == "jaspar"){
    x$LabelBySource <- paste("Jaspar", x$Feature.description, sep = "\n")
  }else if(x$Source[1] == "neph"){
    x$LabelBySource <- paste("Neph", x$Feature.description, sep = "\n")
  }else if(x$Source[1] == "predrem"){
    x$LabelBySource <-  gsub("predrem__", "PreDREM\n", x$Feature.ID)
  }else if(x$Source[1] == "scertf"){
    x$LabelBySource <-  gsub("scertf__(.)+\\.", "ScerTF\n", x$Feature.ID)
  }else if(x$Source[1] == "stark"){
    f <- x$Feature.annotation.for.results1
    f[f == ""] <- x$Feature.description[f==""]
    x$LabelBySource <-  paste("Stark", f, sep = "\n")
  }else if(x$Source[1] == "swissregulon"){
    x$LabelBySource <-  gsub("swissregulon__(.)+__", "SwissRegulon\n", x$Feature.ID)
  }else if(x$Source[1] == "taipale"){
    x$LabelBySource <-  gsub("_(.)+$", "", paste("Taipale", gsub("taipale__", "", x$Feature.ID), sep = "\n"))
    tfpairs <- grepl("\ntaipale", x$LabelBySource) # TF pairs
    x$LabelBySource[tfpairs] <- paste0("Taipale\n",x$Feature.annotation.for.results1[tfpairs])
  }else if(x$Source[1] == "tiffin"){
    x$LabelBySource <-  gsub("tiffin__", "Tiffin\n", x$Feature.ID)
  }else if(x$Source[1] == "transfac"){
    x$LabelBySource <-  paste0("TRANSFAC\n", gsub("(.)+\\: ", "", x$Feature.description))
  }else if(x$Source[1] == "yetfasco"){
    x$LabelBySource <-  paste0("YeTFaSCo\n", gsub("(.)+\\: ", "", x$Feature.annotation.for.results1))
  }
  
  x
})
comparison <- do.call(rbind,comparison)
comparison$LabelBySource_noSource <- gsub("(.)+\n", "", comparison$LabelBySource)
library(ggplot2)
library(cowplot)
library(ggrepel)

gc <- ggplot(comparison, aes(x = Control.NES, y = Infected.NES, label = LabelBySource, Feature.ID = Feature.ID, Feature.description = Feature.description, color = Source)) +
  geom_point(aes()) +
  geom_abline() +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_label_repel(data = comparison[comparison$NES.deviation > quantile(comparison$NES.deviation, 0.90),], alpha = 0.5) +
  theme_cowplot()

# library(plotly)
# ggplotly(gc)


cp <- ggplot(comparison, aes(x = Control.NES, y = Infected.NES, label = LabelBySource, Feature.ID = Feature.ID, Feature.description = Feature.description, color = Source)) +
  geom_point(aes()) +
  geom_abline() +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_label_repel(data = comparison[abs(comparison$InfectedVsControl) > quantile(abs(comparison$InfectedVsControl), 0.90),], alpha = 0.5) +
  theme_cowplot()


comparison$InfectedVsControlRank <- rank(comparison$InfectedVsControl)
comparison$InfectedVsControlAbsRank <- rank(abs(comparison$InfectedVsControl))

nlabel <- 15

comparison$toLabel <- (comparison$InfectedVsControlRank<= nlabel) | 
  (comparison$InfectedVsControlRank >= (max(comparison$InfectedVsControlRank)- (nlabel - 1)) ) |
  (comparison$InfectedVsControlAbsRank <= nlabel )


comparison$Feature.ID.ordered <- factor(comparison$Feature.ID, levels = comparison$Feature.ID[order(comparison$InfectedVsControl)])

gcd <- ggplot(comparison, 
              aes(x = Feature.ID.ordered, 
                  y = InfectedVsControl, 
                  label = LabelBySource_noSource, 
                  Feature.ID = Feature.ID, 
                  Feature.description = Feature.description, 
                  color = Source)) +
  geom_point(aes()) +
  geom_hline(yintercept = 0, color = "grey")  +
  geom_label_repel(data = comparison[comparison$toLabel,], 
                   alpha = 0.8,
                   show.legend = FALSE,
                   size = 3) +
  # geom_label_repel(data = comparison[abs(comparison$InfectedVsControl) < 0.2,], 
  #                  alpha = 1,
  #                  show.legend = FALSE,
  #                  size = 3) +
  theme_cowplot()+
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
        legend.position = "right") +
  guides(color=guide_legend(ncol=1)) + 
  ylab("Infected NES - Control NES\n<-- enriched in Control / enriched in Infected -->")
gcd  + ggsave("./Plots/icistarget-comparison.pdf", width = 9, height = 5, useDingbats = F)


# library(plotly)
# ggplotly(gcd)

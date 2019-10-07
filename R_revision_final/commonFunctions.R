

### For GO analysis
library(GOstats)
library(org.Dm.eg.db)
library(GO.db)
library(ReportingTools)
GO.analyze = function(universe, selected, ontology, report, short.name,outputfile, kable.output, dir = "over", pv = 0.005, minSizeCutoff = 10){
  fb.entrez <- unlist(as.list(org.Dm.egFLYBASE2EG))
  
  
  universeIDs <- fb.entrez[match(universe, names(fb.entrez))]
  universeIDs <- universeIDs[!is.na(universeIDs)]
  selectedIDs <- fb.entrez[match(selected, names(fb.entrez))]
  selectedIDs <- selectedIDs[!is.na(selectedIDs)]
  goParams <- new("GOHyperGParams",
                  geneIds = selectedIDs,
                  universeGeneIds = universeIDs,
                  annotation ="org.Dm.eg" ,
                  ontology = ontology,
                  pvalueCutoff = pv,
                  conditional = TRUE,
                  testDirection = dir,
                  minSizeCutoff = minSizeCutoff)
  goResults <- hyperGTest(goParams)
  
  if(kable.output==TRUE){kable(summary(goResults)[1:10,])}
  if(report == TRUE){
    goReport <- HTMLReport(shortName = short.name,
                           title = outputfile,
                           reportDirectory = "Reports/");
    publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Dm.eg",
            pvalueCutoff= pv)
    finish(goReport)
    return(goReport)
  }
  return(goResults)
  
}


## legacy function for plotting multiple panels - can be replaced with cowplot::plot_grid
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#general plot function
general.plot = function(x, file = NULL){
  g<-gene.ids$Gene.Name[gene.ids$FBgnID == x]
  if(!is.null(file)){
    pdf(file = file, height = 5, width = 5.5)  
  }
  par(mai=c(1.5,0.75,0.1,0.1), mfrow = c(1,2))
  stripchart(y.voom$E[x, ] ~ targets$condres, 
             vertical = TRUE, pch = 19, col=condres.colors, 
             method = "jitter", las=2,
             ylab = paste(g," log2(cpm)"))
  medians <- tapply(y.voom$E[x, ],  targets$condres, median)
  loc <- 1:length(medians)
  segments(loc-0.2, medians, loc+0.2, medians, col="black", lwd=2)
  
  
  stripchart(y.voom$E[x, 39:76]-y.voom$E[x, 1:38] ~ targets$resistance[1:38], 
             vertical = TRUE, pch = 19, col=resistance.colors, 
             method = "jitter", las=2,
             ylab = paste(g," Fold Change"))
  medians <- tapply(y.voom$E[x, 39:76]-y.voom$E[x, 1:38],  targets$resistance[1:38], median)
  loc <- 1:length(medians)
  segments(loc-0.2, medians, loc+0.2, medians, col="black", lwd=2)
  
  
  if(!is.null(file)){dev.off()}
}


### to plot a local-eQTL association
general.plot.eQTL = function (index){ #index is the index in the GRanges
  ci = eQTL.gr[index]
  cols = grep("DGRP", colnames(values(eQTL.gr)))
  g<-ci$gene
  gn<- gene.ids$Gene.Name[gene.ids$FBgnID == as.character(g)]
  par(mai=c(1.5,0.75,0.1,0.1))
  f = factor(as.numeric(as.matrix(values(ci)[,cols])))
  s = factor(c(rep("N",38), rep("T",38))) #infection status
  
  #stripchart format
  # stripchart(y.voom$E[g, ] ~ f:s, 
  #            vertical = TRUE, pch = 19, 
  #            method = "jitter", las=2,
  #            ylab = paste(gn," log2(cpm)"))
  # means <- tapply(t(y.voom$E[g, ]), s:f, mean)
  # segments(1:nlevels(s:f)-0.25, means, 1:nlevels(s:f)+0.25, means, col = "red", lwd = 2)
  
  ##ggplot format
  d = data.frame(expression = y.voom$E[as.character(g),], eQTL.allele = f, treatment = s, class = targets$resistance)
  ggplot(d) + geom_boxplot(aes(x = treatment:eQTL.allele, y = expression, fill = class), alpha = 0.4, outlier.size = -1) +
    geom_point(aes(x = treatment:eQTL.allele, y = expression, fill = class), alpha = 0.4, position = position_jitterdodge(jitter.width = 0.5)) +
    scale_fill_manual(values = resistance.colors) + theme_light() + ylab("Expression level\nLog2(cpm)")

}
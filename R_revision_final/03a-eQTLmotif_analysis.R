library(GenomicRanges)
load(file = "./Data/eQTL.gr", verbose = T)
eQTL.gr$Group <- ifelse(eQTL.gr$ID %in% intersect(eQTL.gr$ID[eQTL.gr$analysis == "Naive"], eQTL.gr$ID[eQTL.gr$analysis == "Treated"]), "Shared", "Unique")
eQTL.gr$Group[eQTL.gr$Group == "Unique"] <- paste0(eQTL.gr$Group[eQTL.gr$Group == "Unique"], "-",eQTL.gr$analysis[eQTL.gr$Group == "Unique"])

icistarget.analyses <- unique(eQTL.gr$Group)


span <- 100 # number of bases to add from each side

for(i in 1:length(icistarget.analyses)){
  file.prefix = paste0("./Data/icistarget_analyses/", icistarget.analyses[i])
  gr.to.test = eQTL.gr[eQTL.gr$Group == icistarget.analyses[i]]
  start(gr.to.test) <- start(gr.to.test) - span
  end(gr.to.test) <- end(gr.to.test) + span
  export.bed(gr.to.test.bed, con = paste0(file.prefix, ".regions.bed"))
}


#' the bed files are then run in https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/
#' The following parameters are used
#' Species: Drosophila melanogaster (dm3)
#' Gene annotation version: Flybase r5.37
#' Input type = BED file
#' Database version = Version 5.0 of databases
#' Features = all features
#' 
#' After completing the Unique-Treated and the Unique-Naive analysis, a comparative analysis is performed here:
#' https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/results_comparison.php
#' The table is manually extracted to a text file and placed in ./Data/icistarget_analyses/i-cisTarget-1Naive2TreatedComparison.txt
#' 
#' The resulting output zip file is extracted and placed in ./Data/icistarget_analyses
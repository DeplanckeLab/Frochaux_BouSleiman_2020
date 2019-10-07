

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

##note to self, I tweeked the dispensability values in order to control the labels. Whenever I want something to show, I add 0.0000
# if I want to hide something, I add 1.0000
revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","userVal_2","userVal_3","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0009409","response to cold", 0.141,-2.702,-4.184, 1.230,-3.7101, 6.000,   NaN,0.932,0.000),
                     c("GO:0015800","acidic amino acid transport", 0.079, 0.407,-6.901, 1.000,-3.1643, 4.000,   NaN,0.829,0.000),
                     c("GO:0017001","antibiotic catabolic process", 0.527, 4.081,-6.352, 1.778,-3.9063, 6.000,   NaN,0.782,0.000),
                     c("GO:0031987","locomotion involved in locomotory behavior", 0.264,-0.903, 0.600, 1.491,-2.3156, 2.000,   NaN,0.956,0.000),
                     c("GO:0031023","microtubule organizing center organization", 1.239,-1.892, 6.595, 2.152,-3.0217,14.000,   NaN,0.829,0.037),
                     c("GO:0035050","embryonic heart tube development", 0.176,-4.912,-0.540, 1.322,-2.8177, 6.000,   NaN,0.905,0.040),
                     c("GO:0006044","N-acetylglucosamine metabolic process", 0.062, 6.821, 4.439, 0.903,-2.8878, 4.000,   NaN,0.836,0.050),
                     c("GO:0006353","DNA-templated transcription, termination", 0.053, 5.432, 5.243, 0.845,-2.8878, 4.000,   NaN,0.812,0.079),
                     c("GO:0005984","disaccharide metabolic process", 0.193, 7.218,-2.786, 1.362,-2.5264, 4.000,   NaN,0.721,0.091),
                     c("GO:0018054","peptidyl-lysine biotinylation", 0.018,-2.541, 4.371, 0.477,-2.3156, 2.000,   NaN,0.857,0.096),
                     c("GO:0071333","cellular response to glucose stimulus", 0.009,-1.027,-3.222, 0.301,-2.3156, 2.000,   NaN,0.849,0.117),
                     c("GO:0046209","nitric oxide metabolic process", 0.026, 1.762, 6.775, 0.602,-2.7919, 3.000,   NaN,0.864,0.121),
                     c("GO:2001057","reactive nitrogen species metabolic process", 0.035, 1.973, 5.557, 0.699,-2.4190, 3.000,   NaN,0.879,0.123),
                     c("GO:0001748","optic lobe placode development", 0.053,-3.909, 0.949, 0.845,-2.2451, 2.000,   NaN,0.919,0.124),
                     c("GO:0007516","hemocyte development", 0.070,-5.507, 0.954, 0.954,-2.1429, 3.000,   NaN,0.874,0.126),
                     c("GO:0000349","generation of catalytic spliceosome for first transesterification step", 0.009, 1.861, 7.343, 0.301,-2.2451, 2.000,   NaN,0.828,0.186),
                     c("GO:0030491","heteroduplex formation", 0.184, 4.157, 6.581, 1.322,-2.3156, 2.000,   NaN,0.828,0.207),
                     c("GO:0009305","protein biotinylation", 0.018,-1.710, 4.378, 0.477,-2.3156, 2.000,   NaN,0.879,0.211),
                     c("GO:0072359","circulatory system development", 0.958,-5.514,-0.445, 2.041,-2.0228,14.000,   NaN,0.902,0.226),
                     c("GO:0070124","mitochondrial translational initiation", 0.009, 3.335, 3.997, 0.301,-2.3156, 2.000,   NaN,0.772,0.234),
                     c("GO:0000052","citrulline metabolic process", 0.092, 6.786, 0.373, 1.041,-2.2451, 2.000,   NaN,0.684,0.242),
                     c("GO:0006393","termination of mitochondrial transcription", 0.009, 5.595, 5.999, 0.301,-2.3156, 2.000,   NaN,0.828,0.244),
                     c("GO:0051234","establishment of localization",17.913,-1.139,-6.655, 3.309,-2.5834,124.000,   NaN,0.918,0.273),
                     c("GO:0010312","detoxification of zinc ion", 0.097,-3.670,-3.667, 1.079,-2.2875, 4.000,   NaN,0.877,0.289),
                     c("GO:0072348","sulfur compound transport", 0.202,-0.560,-6.556, 1.380,-2.9656, 6.000,   NaN,0.874,0.307),
                     c("GO:1902017","regulation of cilium assembly", 0.009,-1.584, 6.987, 0.301,-2.2451, 2.000,   NaN,0.854,0.339),
                     c("GO:0006979","response to oxidative stress", 1.090,-3.315,-3.931, 2.097,-2.3325,16.000,   NaN,0.928,0.352),
                     c("GO:0071577","zinc II ion transmembrane transport", 0.114,-0.327,-7.001, 1.146,-2.0002, 4.000,   NaN,0.899,0.385),
                     c("GO:0043653","mitochondrial fragmentation involved in apoptotic process", 0.018,-3.583, 5.234, 0.477,-2.2451, 2.000,   NaN,0.848,0.399),
                     c("GO:0006105","succinate metabolic process", 0.009, 6.663,-1.645, 0.301,-2.2451, 2.000,   NaN,0.720,0.409),
                     c("GO:0046672","positive regulation of compound eye retinal cell programmed cell death", 0.079,-5.173, 2.002, 1.000,-2.1429, 3.000,   NaN,0.828,0.421),
                     c("GO:0006154","adenosine catabolic process", 0.053, 5.634, 0.452, 0.845,-3.3013, 3.000,   NaN,0.588,0.430),
                     c("GO:0055072","iron ion homeostasis", 0.149, 0.219,-3.204, 1.255,-2.1611, 4.000,   NaN,0.944,0.491),
                     c("GO:0019626","short-chain fatty acid catabolic process", 0.108, 6.062,-2.236, 1.114,-2.8957, 3.000,   NaN,0.566,0.499),
                     c("GO:0090325","regulation of locomotion involved in locomotory behavior", 0.009,-1.060, 1.272, 0.301,-2.3156, 2.000,   NaN,0.942,0.507),
                     c("GO:0042737","drug catabolic process", 0.527, 4.129,-6.207, 1.778,-2.6339,10.000,   NaN,0.782,0.525),
                     c("GO:0046102","inosine metabolic process", 0.244, 6.247, 1.691, 1.447,-3.3013, 3.000,   NaN,0.616,0.579),
                     c("GO:0006643","membrane lipid metabolic process", 0.668, 7.528,-1.769, 1.886,-2.6394,13.000,   NaN,0.665,0.584),
                     c("GO:0016051","carbohydrate biosynthetic process", 0.439, 7.649,-0.250, 1.708,-2.0631, 8.000,   NaN,0.682,0.586),
                     c("GO:0046327","glycerol biosynthetic process from pyruvate", 0.053, 7.065,-0.720, 0.845,-2.3156, 2.000,   NaN,0.574,0.588),
                     c("GO:0046103","inosine biosynthetic process", 0.067, 6.387, 1.857, 0.903,-3.3013, 3.000,   NaN,0.611,0.590),
                     c("GO:0046085","adenosine metabolic process", 0.070, 6.587, 1.893, 0.954,-3.3013, 3.000,   NaN,0.643,0.592),
                     c("GO:0030148","sphingolipid biosynthetic process", 0.149, 6.839, 0.855, 1.255,-2.3376, 7.000,   NaN,0.627,0.600),
                     c("GO:0032787","monocarboxylic acid metabolic process", 1.380, 7.165,-1.221, 2.199,-2.3685,19.000,   NaN,0.622,0.601),
                     c("GO:0009311","oligosaccharide metabolic process", 0.457, 7.150,-2.495, 1.724,-2.0349, 6.000,   NaN,0.722,0.633),
                     c("GO:0008610","lipid biosynthetic process", 1.679, 7.288,-0.093, 2.283,-2.4131,22.000,   NaN,0.630,0.671));



one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);

#### Add the number of genes per group

n.genes = read.table("Reports/GO.eQTLs.pooled.for.revigo.txt", header = T)

colnames(n.genes) <- c("GO", "p-value", "n.genes", "analysis")
one.data$n.genes = n.genes$n.genes[match(one.data$term_ID, n.genes$GO)]
one.data$analysis = n.genes$analysis[match(one.data$term_ID, n.genes$GO)]
one.data$analysis <- plyr::mapvalues(one.data$analysis, c("naive", "treated", "common"), c("Control", "Infected", "Shared"))

one.data$analysis[is.na(one.data$analysis)] = "Shared"
# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, fill = log10_p_value, color = analysis, size = n.genes, ID = description),  shape = 21, stroke = 2, alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_fill_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + scale_color_manual(values = c("#727272", "#026f91", "#ef926a"))
# p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = n.genes, shape = analysis), fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1 <- p1 + guides(shape = guide_legend(override.aes = list(size=5)))




# --------------------------------------------------------------------------
# Output the plot to screen

p1;

pdf(file = "./Plots/Revigo.eQTL.pdf", width = 9, height = 6)
p1;
dev.off()

library(plotly)
ggplotly(p1)

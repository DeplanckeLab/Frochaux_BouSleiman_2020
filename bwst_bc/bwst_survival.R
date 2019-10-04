###########################
### ms771 x w- Survival ###
###########################

### START ####

## data ####
# path <- "C:/Users/Michael/Documents/Deplancke/data/ntc/bwst Survival/"
# fig.path <- file.path(path, "Figures")
replica1 <- "./Data/160803_KM.txt"
replica2 <- "./Data/160817_KM.txt"
replica3 <- "./Data/160824_KM.txt"

## colors ####
col.palette <- c("black", "goldenrod") #control and mutant, black 0,0,0, goldenrod, r 218, g 165, b 32

## libraries ####
# library("Hmisc") # needed for fct errbar
library("survival") # KM and cox model
library("ggplot2") # nicer graphs
# library("broom") # tidy
# library("GGally")
library(survminer) # for ggsurvplot
library(reshape2)
library(cowplot)

## functions ####
KMrename <- function(data.frame) {
  data.frame$short <- NA
  data.frame$short[data.frame$Genotype=="w-"] <- "w"
  data.frame$short[data.frame$Genotype=="rel"] <- "rel"
  data.frame$short[data.frame$Genotype=="ms771/tm6b"] <- "ms.tm"
  data.frame$short[data.frame$Genotype=="bw;st;+/ms771"] <- "ms"
  data.frame$short[data.frame$Genotype=="bw;st;+/tm6b"] <- "tm"
  data.frame$short[data.frame$Genotype=="bw;st"] <- "bwst"
  return(data.frame)
}
### END START ####

##################################
### data analysis and plotting ###
##################################



### 1. data preprocessing ####
## 1.1 load data ####
KM.bio1 <- read.table(replica1, header=T, sep="\t")
KM.bio2 <- read.table(replica2, header=T, sep="\t")
KM.bio3 <- read.table(replica3, header=T, sep="\t")

## 1.2 add column ####
KM.bio1 <- KMrename(KM.bio1)
KM.bio2 <- KMrename(KM.bio2)
KM.bio3 <- KMrename(KM.bio3)

## 1.3 combine replica ####
KM.all <- rbind(KM.bio1, KM.bio2, KM.bio3)
KM.all.pe <- subset(KM.all,  Treatment=="pe")
KM.pe.TM <- subset(KM.all.pe, KM.all.pe$short=="ms.tm" | KM.all.pe$short=="tm")
KM.pe.wt <- subset(KM.all.pe, KM.all.pe$short=="ms" | KM.all.pe$short=="bwst")
KM.ctrl <- subset(KM.all.pe, KM.all.pe$short=="bwst" | KM.all.pe$short=="tm")
### END 1. data preprocessing ------------------------------------------------------------------------------------

### 2. fit survival model using Genotype as factor, extract values for plotting and save as RDS ####

## 2.1.1 strains with balancer ####
attach(KM.pe.TM)
TM.fit <- survfit(Surv(Time, Status) ~ Genotype, data=KM.pe.TM)
detach(KM.pe.TM)

## 2.1.2 transform survfit object into data frame ####
TM.summary <- summary(TM.fit)
TM.df <- as.data.frame(rbind(TM.summary$time, TM.summary$surv, 
                               TM.summary$upper, TM.summary$lower, gsub("Genotype=", "", as.character(TM.summary$strata))))
row.names(TM.df) <- c("time", "surv", "upper", "lower", "Genotype")

## 2.1.3 add two columns for starting time ####
geno <- gsub("Genotype=", "", levels(TM.summary$strata))
TM.df$t1 <- c(0,1,1,1, geno[1])
TM.df$t2 <- c(0,1,1,1, geno[2])

## 2.1.4 take transform of data frame and separate by genotype ####
surv.tm <- as.data.frame(t(TM.df))
tm.ctrl <- surv.tm[surv.tm$Genotype == "bw;st;+/tm6b",]
tm.ctrl <- tm.ctrl[order(as.numeric(as.character(tm.ctrl$time))),]
tm.test <- surv.tm[surv.tm$Genotype != "bw;st;+/tm6b",]
tm.test <- tm.test[order(as.numeric(as.character(tm.test$time))),]

## 2.1.5 saving data set -----------------------------------------------------------------------------------------
saveRDS(tm.ctrl, "./Processed/TM_ctrl.RDS")
saveRDS(tm.test, "./Processed/TM_test.RDS")

## 2.2.1 strains without balancer ####
attach(KM.pe.wt)
wt.fit <- survfit(Surv(Time, Status) ~ Genotype, data=KM.pe.wt)
detach(KM.pe.wt)

## 2.2.2 transform survfit object into data frame ####
wt.summary <- summary(wt.fit)
wt.df <- as.data.frame(rbind(wt.summary$time, wt.summary$surv, 
                             wt.summary$upper, wt.summary$lower, gsub("Genotype=", "", as.character(wt.summary$strata))))
row.names(wt.df) <- c("time", "surv", "upper", "lower", "Genotype")

## 2.2.3 add two columns for starting time ####
geno <- gsub("Genotype=", "", levels(wt.summary$strata))
wt.df$t1 <- c(0,1,1,1, geno[1])
wt.df$t2 <- c(0,1,1,1, geno[2])

## 2.2.4 take transform of data frame and separate by genotype ####
surv.wt <- as.data.frame(t(wt.df))
wt.ctrl <- surv.wt[surv.wt$Genotype == "bw;st",]
wt.ctrl <- wt.ctrl[order(as.numeric(as.character(wt.ctrl$time))),]
wt.test <- surv.wt[surv.wt$Genotype != "bw;st",]
wt.test <- wt.test[order(as.numeric(as.character(wt.test$time))),]

## 2.2.5 saving data set -----------------------------------------------------------------------------------------
saveRDS(wt.ctrl, "./Processed/wt_ctrl.RDS")
saveRDS(wt.test, "./Processed/wt_test.RDS")

## 2.3.1 control strain ####
attach(KM.ctrl)
ctrl.fit <- survfit(Surv(Time, Status) ~ Genotype, data=KM.ctrl)
detach(KM.ctrl)

## 2.3.2 transform survfit object into data frame ####
ctrl.summary <- summary(ctrl.fit)
ctrl.df <- as.data.frame(rbind(ctrl.summary$time, ctrl.summary$surv, 
                               ctrl.summary$upper, ctrl.summary$lower, gsub("Genotype=", "", as.character(ctrl.summary$strata))))
row.names(ctrl.df) <- c("time", "surv", "upper", "lower", "Genotype")

## 2.3.3 add two columns for starting time ####
geno <- gsub("Genotype=", "", levels(ctrl.summary$strata))
ctrl.df$t1 <- c(0,1,1,1, geno[1])
ctrl.df$t2 <- c(0,1,1,1, geno[2])

## 2.3.4 take transform of data frame and separate by genotype ####
surv.ctrl <- as.data.frame(t(ctrl.df))
ctrl1 <- surv.ctrl[surv.ctrl$Genotype == "bw;st",]
ctrl1 <- ctrl1[order(as.numeric(as.character(ctrl1$time))),]
ctrl2 <- surv.ctrl[surv.ctrl$Genotype != "bw;st",]
ctrl2 <- ctrl2[order(as.numeric(as.character(ctrl2$time))),]

## 2.3.5 saving data set -----------------------------------------------------------------------------------------
saveRDS(ctrl1, "./Processed/ctrl_test_1.RDS")
saveRDS(ctrl2, "./Processed/ctrl_test_2.RDS")

## 2.3 clean up ####
rm("geno", "KM.all", "KM.all.pe", "KM.bio1", "KM.bio2", "KM.bio3", "KM.pe.TM", "KM.pe.wt", 
   "surv.tm", "surv.wt", "TM.df", "TM.fit", "TM.summary", "wt.df", "wt.fit", "wt.summary")
### END 2. statistical model -------------------------------------------------------------------------------------

### 3. Plotting ####
## 3.1 load data from 2.x.5 if needed ####
wt.ctrl <- readRDS("./Processed/wt_ctrl.RDS")
wt.test <- readRDS("./Processed/wt_test.RDS")
tm.ctrl <- readRDS("./Processed/TM_ctrl.RDS")
tm.test <- readRDS("./Processed/TM_test.RDS")
ctrl1 <- readRDS("./Processed/ctrl_test_1.RDS")
ctrl2 <- readRDS("./Processed/ctrl_test_2.RDS")

## 3.2 colors and title ####
col.ctrl <- rgb(130,130,169, 80, maxColorValue = 255)
col.test <- rgb(218,165,32, 80, maxColorValue = 255)
bg.col <- rgb(230,230,230, 40, maxColorValue = 255)
main.t <- bquote(.("Survival of") ~ italic("ntc")^italic("ms771") ~ "mutants")

## 3.3 strains with balancer ####
# 3.3.1 line vectors
x.v <- as.numeric(as.character(tm.test$time))
y.v <- as.numeric(as.character(tm.test$surv))
x.v2 <- as.numeric(as.character(tm.ctrl$time))
y.v2 <- as.numeric(as.character(tm.ctrl$surv))


# 3.3.2 error vectors
err.x <- c(x.v, rev(x.v))
err.y <- c(as.numeric(as.character(tm.test$lower)), rev(as.numeric(as.character(tm.test$upper))))
err.x2 <- c(x.v2, rev(x.v2))
err.y2 <- c(as.numeric(as.character(tm.ctrl$lower)), rev(as.numeric(as.character(tm.ctrl$upper))))


# 3.3.3 Plot
# df for line
line.df1 <- rbind(x.v, y.v)
row.names(line.df1) <- c("time", "test")
line.df1 <- t(line.df1)
line.df1 <- as.data.frame(line.df1)

line.df2 <- rbind(x.v2, y.v2)
row.names(line.df2) <- c("time", "test")
line.df2 <- t(line.df2)
line.df2 <- as.data.frame(line.df2)

test.df <- rbind(x.v, as.numeric(as.character(tm.test$upper)), as.numeric(as.character(tm.test$lower)))
ctrl.df <- rbind(x.v2, as.numeric(as.character(tm.ctrl$upper)), as.numeric(as.character(tm.ctrl$lower)))
row.names(test.df) <- c("time", "upper", "lower")
test.df <- t(test.df)
test.df <- as.data.frame(test.df)
row.names(ctrl.df) <- c("time", "upper", "lower")
ctrl.df <- t(ctrl.df)
ctrl.df <- as.data.frame(ctrl.df)

p1 <- ggplot() +
  geom_line(data=line.df1, aes(x=time, y=test), colour="goldenrod") +
  geom_line(data=line.df2, aes(x=time, y=test), colour="black") +  
  geom_ribbon(data=test.df, aes(x=time, ymin=lower, ymax=upper), fill="goldenrod", alpha=0.4) +
  geom_ribbon(data=ctrl.df, aes(x=time, ymin=lower, ymax=upper), fill="grey", alpha=0.4) +
  theme(legend.position = "none", 
        panel.border=element_rect(colour="black", fill=NA, size=1),
        panel.background=element_rect(fill=NA),
        axis.ticks= element_line(size=1),
        axis.ticks.length = unit(6, "points"),
        axis.text = element_text(size=20),
        axis.title = element_blank())

print(p1)

save_plot(filename = "./Figures/tm_simple.pdf", p1)
### END 3. Plotting ####

## 3.4 strains without balancer ####
# 3.4.1 line vectors
x.v <- as.numeric(as.character(wt.test$time))
y.v <- as.numeric(as.character(wt.test$surv))
x.v2 <- as.numeric(as.character(wt.ctrl$time))
y.v2 <- as.numeric(as.character(wt.ctrl$surv))

# 3.4.2 error vectors
err.x <- c(x.v, rev(x.v))
err.y <- c(as.numeric(as.character(wt.test$lower)), rev(as.numeric(as.character(wt.test$upper))))
err.x2 <- c(x.v2, rev(x.v2))
err.y2 <- c(as.numeric(as.character(wt.ctrl$lower)), rev(as.numeric(as.character(wt.ctrl$upper))))

# df for line
line.df <- rbind(x.v, y.v, y.v2)
row.names(line.df) <- c("time", "test", "ctrl")
line.df <- t(line.df)
line.df <- as.data.frame(line.df)
line.df <- melt(data=line.df, id.vars="time", measure.vars=c("test", "ctrl"))

test.df <- rbind(x.v, as.numeric(as.character(wt.test$upper)), as.numeric(as.character(wt.test$lower)))
ctrl.df <- rbind(x.v, as.numeric(as.character(wt.ctrl$upper)), as.numeric(as.character(wt.ctrl$lower)))
row.names(test.df) <- c("time", "upper", "lower")
test.df <- t(test.df)
test.df <- as.data.frame(test.df)
row.names(ctrl.df) <- c("time", "upper", "lower")
ctrl.df <- t(ctrl.df)
ctrl.df <- as.data.frame(ctrl.df)

p2 <- ggplot(data=line.df) +
  geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(values = c("goldenrod", "black")) +  
  geom_ribbon(data=test.df, aes(x=time, ymin=lower, ymax=upper), fill="goldenrod", alpha=0.4) +
  geom_ribbon(data=ctrl.df, aes(x=time, ymin=lower, ymax=upper), fill="grey", alpha=0.4) +
  theme(legend.position = "none", 
    panel.border=element_rect(colour="black", fill=NA, size=1),
    panel.background=element_rect(fill=NA),
    axis.ticks= element_line(size=1),
    axis.ticks.length = unit(6, "points"),
    axis.text = element_text(size=20),
    axis.title = element_blank())

print(p2)

save_plot(filename = "./Figures/wt_simple.pdf", p2)
### END 3. Plotting ####

## 3.5 comparison between reference strains ####
# 3.4.1 line vectors
x.v <- as.numeric(as.character(ctrl1$time))
y.v <- as.numeric(as.character(ctrl1$surv))
x.v2 <- as.numeric(as.character(ctrl2$time))
y.v2 <- as.numeric(as.character(ctrl2$surv))

# 3.4.2 error vectors
err.x <- c(x.v, rev(x.v))
err.y <- c(as.numeric(as.character(ctrl1$lower)), rev(as.numeric(as.character(ctrl1$upper))))
err.x2 <- c(x.v2, rev(x.v2))
err.y2 <- c(as.numeric(as.character(ctrl2$lower)), rev(as.numeric(as.character(ctrl2$upper))))

# 3.4.3 Plot
# df for line
line.df <- rbind(x.v, y.v, y.v2)
row.names(line.df) <- c("time", "test", "ctrl")
line.df <- t(line.df)
line.df <- as.data.frame(line.df)
line.df <- melt(data=line.df, id.vars="time", measure.vars=c("test", "ctrl"))

test.df <- rbind(x.v, as.numeric(as.character(ctrl1$upper)), as.numeric(as.character(ctrl1$lower)))
ctrl.df <- rbind(x.v, as.numeric(as.character(ctrl2$upper)), as.numeric(as.character(ctrl2$lower)))
row.names(test.df) <- c("time", "upper", "lower")
test.df <- t(test.df)
test.df <- as.data.frame(test.df)
row.names(ctrl.df) <- c("time", "upper", "lower")
ctrl.df <- t(ctrl.df)
ctrl.df <- as.data.frame(ctrl.df)

p3 <- ggplot(data=line.df) +
  geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(values = c("goldenrod", "black")) +  
  geom_ribbon(data=test.df, aes(x=time, ymin=lower, ymax=upper), fill="goldenrod", alpha=0.4) +
  geom_ribbon(data=ctrl.df, aes(x=time, ymin=lower, ymax=upper), fill="grey", alpha=0.4) +
  theme(legend.position = "none", 
        panel.border=element_rect(colour="black", fill=NA, size=1),
        panel.background=element_rect(fill=NA),
        axis.ticks= element_line(size=1),
        axis.ticks.length = unit(6, "points"),
        axis.text = element_text(size=20),
        axis.title = element_blank())

print(p3)

save_plot(filename = "./Figures/comparison_ctrl.pdf", p3)
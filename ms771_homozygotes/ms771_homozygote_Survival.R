####################################################################
### ANALYSIS OF THE SURVIVAL OF NTC HOMOZYGOTE NULL MUTANT ms771 ###
####################################################################

### Start ####
## colors ####
col.palette <- c("black", "goldenrod") #control and mutant, black 0,0,0, goldenrod, r 218, g 165, b 32

## data ####
replica1 <- "./Data/151130_KM.txt"
replica2 <- "./Data/151202_KM.txt"
replica3 <- "./Data/151216_KM.txt"
replica4 <- "./Data/160211_KM.txt"

## libraries ####
library("Hmisc") # needed for fct errbar
library("survival") # KM and cox model
library(ggplot2)
library(cowplot)


### END Start ####

##################################
### data analysis and plotting ###
##################################

### 1. data preprocessing ####
## 1.1 load data ####
bio1 <- read.delim(replica1, header=T)
bio2 <- read.delim(replica2, header=T)
bio3 <- read.delim(replica3, header=T)
bio4 <- read.delim(replica4, header=T)

## 1.2 combine replica ####
ms.all <- rbind(bio1, bio2, bio3, bio4)
ms.pe <- ms.all[ms.all$Treatment=="pe",]
ms.uc <- ms.all[ms.all$Treatment=="uc",]

### END 1. data preprocessing ------------------------------------------------------------------------------------

### 2. fit survival model using Genotype as factor, extract values for plotting and save as RDS ####

## 2.1.1 Unchallenged ####
attach(ms.uc)
uc.fit <- survfit(Surv(Time, Status) ~ Genotype, data=ms.uc)
detach(ms.uc)

## 2.1.2 transform survfit object into data frame ####
uc.summary <- summary(uc.fit)
uc.df <- as.data.frame(rbind(uc.summary$time, uc.summary$surv, 
                             uc.summary$upper, uc.summary$lower, gsub("Genotype=", "", as.character(uc.summary$strata))))
row.names(uc.df) <- c("time", "surv", "upper", "lower", "Genotype")

## 2.1.3 add two columns for starting time ####
geno <- gsub("Genotype=", "", levels(uc.summary$strata))
uc.df$t1 <- c(0,1,1,1, geno[1])
uc.df$t2 <- c(0,1,1,1, geno[2])

## 2.1.4 take transform of data frame and separate by genotype ####
surv.uc <- as.data.frame(t(uc.df), stringsAsFactors = FALSE)
uc.ctrl <- surv.uc[surv.uc$Genotype == "w-",]
uc.ctrl <- uc.ctrl[order(as.numeric(as.character(uc.ctrl$time))),]
uc.test <- surv.uc[surv.uc$Genotype != "w-",]
uc.test <- uc.test[order(as.numeric(as.character(uc.test$time))),]

## 2.1.5 saving data set -----------------------------------------------------------------------------------------
saveRDS(uc.ctrl, "./Processed/uc_ctrl.RDS")
saveRDS(uc.test, "./Processed/uc_test.RDS")

## 2.2.1 strains without balancer ####
attach(ms.pe)
pe.fit <- survfit(Surv(Time, Status) ~ Genotype, data=ms.pe)
detach(ms.pe)

## 2.2.2 transform survfit object into data frame ####
pe.summary <- summary(pe.fit)
pe.df <- as.data.frame(rbind(pe.summary$time, pe.summary$surv, 
                             pe.summary$upper, pe.summary$lower, gsub("Genotype=", "", as.character(pe.summary$strata))))
row.names(pe.df) <- c("time", "surv", "upper", "lower", "Genotype")

## 2.2.3 add two columns for starting time ####
geno <- gsub("Genotype=", "", levels(pe.summary$strata))
pe.df$t1 <- c(0,1,1,1, geno[1])
pe.df$t2 <- c(0,1,1,1, geno[2])

## 2.2.4 take transform of data frame and separate by genotype ####
surv.pe <- as.data.frame(t(pe.df))
pe.ctrl <- surv.pe[surv.pe$Genotype == "w-",]
pe.ctrl <- pe.ctrl[order(as.numeric(as.character(pe.ctrl$time))),]
pe.test <- surv.pe[surv.pe$Genotype != "w-",]
pe.test <- pe.test[order(as.numeric(as.character(pe.test$time))),]

## 2.2.5 saving data set -----------------------------------------------------------------------------------------
saveRDS(pe.ctrl, "./Processed/pe_ctrl.RDS")
saveRDS(pe.test, "./Processed/pe_test.RDS")

## 2.3 clean up ####
rm("geno", "ms.all", "bio1", "bio2", "bio3", "bio4", "ms.uc", "ms.pe", 
   "surv.uc", "surv.pe", "uc.df", "uc.fit", "uc.summary", "pe.df", "pe.fit", "pe.summary", pe.ctrl, pe.test, uc.ctrl, uc.test)
### END 2. statistical model -------------------------------------------------------------------------------------

### 3. Plotting ####
## 3.1 load data from 2.x.5 if needed ####
pe.ctrl <- readRDS("./Processed/pe_ctrl.RDS")
pe.test <- readRDS("./Processed/pe_test.RDS")
uc.ctrl <- readRDS("./Processed/uc_ctrl.RDS")
uc.test <- readRDS("./Processed/uc_test.RDS")

# rm last timepoint from uc.test
uc.test <- uc.test[-5,]

## 3.2 colors and title ####
col.ctrl <- rgb(190,190,190, 80, maxColorValue = 255)
col.test <- rgb(218,165,32, 80, maxColorValue = 255)
bg.col <- rgb(230,230,230, 40, maxColorValue = 255)
main.t <- bquote(.("Survival of") ~ italic("ntc")^italic("ms771") ~ "homozygotes compared to" ~ italic("w")^italic("1118"))

## 3.3 unchallenged ####
# 3.3.1 line vectors
x.v <- as.numeric(as.character(uc.test$time))
y.v <- as.numeric(as.character(uc.test$surv))
x.v2 <- as.numeric(as.character(uc.ctrl$time))
y.v2 <- as.numeric(as.character(uc.ctrl$surv))

# 3.3.2 error vectors
err.x <- c(x.v, rev(x.v))
err.y <- c(as.numeric(as.character(uc.test$lower)), rev(as.numeric(as.character(uc.test$upper))))
err.x2 <- c(x.v2, rev(x.v2))
err.y2 <- c(as.numeric(as.character(uc.ctrl$lower)), rev(as.numeric(as.character(uc.ctrl$upper))))

# df for line
line.df1 <- rbind(x.v, y.v)
row.names(line.df1) <- c("time", "test")
line.df1 <- t(line.df1)
line.df1 <- as.data.frame(line.df1)

line.df2 <- rbind(x.v2, y.v2)
row.names(line.df2) <- c("time", "test")
line.df2 <- t(line.df2)
line.df2 <- as.data.frame(line.df2)

test.df <- rbind(x.v, as.numeric(as.character(uc.test$upper)), as.numeric(as.character(uc.test$lower)))
ctrl.df <- rbind(x.v2, as.numeric(as.character(uc.ctrl$upper)), as.numeric(as.character(uc.ctrl$lower)))
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

save_plot("./Figures/homozygous_uc.pdf", p1)

## 3.4 pe ####
# 3.4.1 line vectors
x.v <- as.numeric(as.character(pe.test$time))
y.v <- as.numeric(as.character(pe.test$surv))
x.v2 <- as.numeric(as.character(pe.ctrl$time))
y.v2 <- as.numeric(as.character(pe.ctrl$surv))

# 3.4.2 error vectors
err.x <- c(x.v, rev(x.v))
err.y <- c(as.numeric(as.character(pe.test$lower)), rev(as.numeric(as.character(pe.test$upper))))
err.x2 <- c(x.v2, rev(x.v2))
err.y2 <- c(as.numeric(as.character(pe.ctrl$lower)), rev(as.numeric(as.character(pe.ctrl$upper))))

# 3.4.3 Plot
# df for line
line.df1 <- rbind(x.v, y.v)
row.names(line.df1) <- c("time", "test")
line.df1 <- t(line.df1)
line.df1 <- as.data.frame(line.df1)

line.df2 <- rbind(x.v2, y.v2)
row.names(line.df2) <- c("time", "test")
line.df2 <- t(line.df2)
line.df2 <- as.data.frame(line.df2)

test.df <- rbind(x.v, as.numeric(as.character(pe.test$upper)), as.numeric(as.character(pe.test$lower)))
ctrl.df <- rbind(x.v2, as.numeric(as.character(pe.ctrl$upper)), as.numeric(as.character(pe.ctrl$lower)))
row.names(test.df) <- c("time", "upper", "lower")
test.df <- t(test.df)
test.df <- as.data.frame(test.df)
row.names(ctrl.df) <- c("time", "upper", "lower")
ctrl.df <- t(ctrl.df)
ctrl.df <- as.data.frame(ctrl.df)

p2 <- ggplot() +
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

print(p2)

save_plot("./Figures/homozygous_pe.pdf", p2)
### END 3. Plotting ####
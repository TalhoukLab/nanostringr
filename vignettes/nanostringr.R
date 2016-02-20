## ----setupQualityMeasures, message = FALSE, echo = FALSE, warning = FALSE----
library(knitr)
library(dplyr)
opts_chunk$set(message = FALSE, echo = FALSE, warning = FALSE)


# Define colors constants----
COL.OVD <- "#66C2A5"
COL.OVO <- "#A6D854"
COL.OVCL <- "#FC8D62"
COL.HLD <-  "#8DA0CB"
COL.HLO <-"#E78AC3"

## ----message=TRUE,echo=TRUE----------------------------------------------
library(nanostringr)
expOVD <- NanoStringQC(ovd.r,subset(expQC,OVD=="Yes"))
expOVO <- NanoStringQC(ovo.r,subset(expQC,OVO=="Yes"))
expOVCL <- NanoStringQC(ovc.r,subset(expQC,OVCL=="Yes"))
expHLD <- NanoStringQC(hld.r,subset(expQC,HLD=="Yes"))
expHLO <- NanoStringQC(hlo.r,subset(expQC,HLO=="Yes"))
expQC <- rbind(expHLD,expOVD,expHLO,expOVO,expOVCL)
expQC$cohort <- factor(c(rep("HLD", nrow(expHLD)), 
                rep("OVD", nrow(expOVD)),rep("HLO",nrow(expHLO)),
                rep("OVO", nrow(expOVO)), rep("OVCL", nrow(expOVCL))))
expQC <- expQC %>% 
  mutate(cohort = factor(stringr::str_replace_all(cohort,
                                                   c("HLD" = "HL",
                                                     "OVD" = "OC")),
                          levels = c("HL", "OC", "OVCL", "HLO", "OVO")))


## ----perFOVPlot, fig.cap =  "Samples that failed imaging QC based on percent fields of view (FOV) counted across cohorts"----

boxplot(perFOV~cohort, ylab = "% FOV", main = "% FOV by Cohort", data = expQC, pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO));
abline(h = 75, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)


## ----linPCPlot, fig.cap ="Plot of $R^2$ of postive control probles from samples across all cohorts."----
boxplot(linPC~cohort, ylab = expression(R^2), main = "Linearity of Positive Controls by Cohort", data = expQC, pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO), ylim=c(0,1));
abline(h = 0.95, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)

## ----averageHKPlot, fig.cap ="Average log expression of Housekeeping genes by Cohort."----
boxplot(averageHK~cohort, ylab = "Average log HK expression", main = "Average log expression of Housekeeping genes by Cohort", data = expQC, pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO));
abline(h = 50, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)

## ----lodPlot, fig.cap ="Limit of detection by cohort."-------------------
boxplot(lod~cohort, ylab = "LOD", main = "Limit of detection (LOD) by Cohort", data = expQC, pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO));
abline(h = 50, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)

## ----pergdPlot,fig.cap ="Percent genes of total (excluding controls) detected above the limit of detection"----
boxplot(pergd~cohort, data = expQC, border = "white",
       ylab = "% Genes Detected", 
       main = "Percent of Genes Detected Above \n the Limit of Detection", 
       pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO));
abline(h = 50, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)
stripchart(pergd ~ cohort, data = expQC, 
           vertical = TRUE, method = "jitter",  
           pch = 20, cex = 0.4 , col = "#3A6EE3", 
           add = TRUE) 

## ----snPlot,fig.cap ="Signal to Noise versus % Gene Detected by cohort"----
sn <- 100
detect <- 60

plot(expOVD$sn, expOVD$pergd, pch = 20, col = COL.OVD, xaxt = "n", ylim = c(0,100), xlim = range(expOVD$sn),
     xlab = "Signal to Noise Ratio",ylab = "% Genes Detected")
points(expOVO$sn, expOVO$pergd, pch = 20, col = COL.OVO)
points(expOVCL$sn, expOVCL$pergd, pch = 20, col = COL.OVCL)
points(expHLD$sn, expHLD$pergd, pch = 20, col = COL.HLD)
points(expHLO$sn, expHLO$pergd, pch = 20, col = COL.HLO)
axis(1, at = seq(0, max(expQC$sn) + 1, 300))
abline(v = sn, col = "red", lwd = 2)
abline(h = detect, lty = 2)
title("Signal to Noise vs \n Ratio of Genes Detected")
legend("bottomright", c("HL","OC","OVCL","HLO","OVO"), pch = 20, bty = 'n', col = c(COL.HLD,COL.OVD, COL.OVCL,COL.HLO, COL.OVO))

## ----snZoom, fig.cap ="Signal to Noise versus % Gene Detected by cohort zoomed in to the area of possible failures"----

plot(expOVD$sn, expOVD$pergd, pch = 20, col = COL.OVD, xaxt = "n", ylim = c(0,100), xlim = c(0,6000),
     xlab = "Signal to Noise Ratio ",ylab = "Ratio of Genes Detected")
points(expOVO$sn, expOVO$pergd, pch = 20, col = COL.OVO)
points(expOVCL$sn, expOVCL$pergd, pch = 20, col = COL.OVCL)
points(expHLD$sn, expHLD$pergd, pch = 20, col = COL.HLD)
points(expHLO$sn, expHLO$pergd, pch = 20, col = COL.HLO)
axis(1, at = seq(0, max(expQC$sn) + 1, 300))

abline(v = sn, col = "red", lwd = 2)
abline(h = detect, lty = 2)
title("Signal to Noise vs \n Ratio of Genes Detected (Zooming-in)")
legend("bottomright", c("HL","OC","OVCL","HLO","OVO"), pch = 20, bty = 'n', col = c(COL.HLD,COL.OVD, COL.OVCL,COL.HLO, COL.OVO))


## ----bdPlot, fig.cap = "Binding density by cohort. Samples outside of the dashed lines are flagged for having failed the binding density QC measure"----
boxplot(binding.density~cohort, ylab = "Binding Density", main = "Binding Density (BD) by Cohort", data = expQC, pch = 20, col = c(COL.HLD, COL.OVD, COL.OVCL,COL.HLO,COL.OVO));
abline(h = 0.05, lty = 2, col = "red")
abline(h = 2.25, lty = 2, col = "red")
grid(NULL, NULL, lwd = 1)

## ----lodbd, echo = FALSE-------------------------------------------------
plot(expQC$lod, type = "l", lwd = 1, col = "coral",
     ylab = "LOD", xlab = "Samples", main = "LOD and Binding Density")
points(which(expQC$bdFlag == "Failed"), expQC$lod[expQC$bdFlag == "Failed"], col = "blue", pch = 19)
grid(NA, 5, lwd = 1)
legend("topleft", legend = "Failed Binding Density QC", col = "blue", pch = 19, bty = 'n')

plot(expQC$binding.density, expQC$lod, pch = 20, col = "deepskyblue", lwd = 1, xlim = c(0,5), ylim = c(-6,60), xaxt = "n", panel.first = grid(NULL,NULL, lwd = 1), ylab = "LOD", xlab = "Binding Density", main = "LOD vs Binding Density")
axis(1, at = seq(0, 5, by = .5), las = 2)
points(expQC$binding.density[expQC$bdFlag == "Failed"], expQC$lod[expQC$bdFlag == "Failed"], col = "blue", pch = 19)
legend("topright", legend = "Failed Binding Density QC", col = "blue", pch = 20, bty = 'n')

## ----setupNormalization, message = FALSE, echo = FALSE, warning = FALSE----
library(nanostringr)
library(NanoStringNorm)
library(dplyr)
library(knitr)
library(batchAdj)
opts_chunk$set(message = FALSE, echo = FALSE, warning = FALSE)

expOVD <- NanoStringQC(ovd.r,subset(expQC,OVD=="Yes"))


## ----fig.height=11, fig.width=7------------------------------------------
HK=t(ovd.r[ovd.r$Code.Class=="Housekeeping",-(1:3)])
par(mfrow=c(5,2))
genes=colnames(HK)[order(apply(HK,2,mean))]
for(i in 1: length(genes)){
  hist(HK[,genes[i]], xlab=genes[i], main="Raw", probability = T)
  hist(log(HK[,genes[i]]), xlab=genes[i], main="Log Base 2", probability = T)
}

## ----NanoStringNorm, results = 'hide'------------------------------------
#Normalized to Positive controls and housekeeping genes using geometric mean
MAPC.d <- NanoStringNorm(ovd.r, CodeCount = 'geo.mean', SampleContent = 'housekeeping.geo.mean', round.values = TRUE, take.log = T)
normPC <- t(MAPC.d$normalized.data[MAPC.d$normalized.data$Code.Class == "Endogenous", -c(1:3)])

#Normalized to housekeeping genes using geometric mean
MAHK.d <- NanoStringNorm(ovd.r, SampleContent = 'housekeeping.geo.mean', round.values = TRUE, take.log = T)
normHK <- t(MAHK.d$normalized.data[MAHK.d$normalized.data$Code.Class == "Endogenous", -c(1:3)])

#Normalized to housekeeping genes using geometric mean and background corrected
MANC.d = NanoStringNorm(ovd.r,Background = 'mean.2sd' ,SampleContent='housekeeping.geo.mean',round.values=TRUE, take.log=T)
normNC=t(MANC.d$normalized.data[MANC.d$normalized.data$Code.Class == "Endogenous", -c(1:3)])

set.seed(13)
gene="Gene 49 "

## ----PCNorm, fig.height = 8, fig.width = 6, caption = "Ovarian Cancer data comparisons of two normalization procedures of the data of the gene `r gene`. The method that performs only housekeeping normalization is denoted by HK1 and HK2 in CodeSet 1 and 2 respectively. On the other hand, the method that normalizes both to positive controls and housekeeping genes at the same time, is denoted by PC1 and PC2 when ran on CodeSet 1 and CodeSet 2 respectively."----

pc1 <- normPC[expOVD$geneRLF=="CS1", gene]
pc2 <- normPC[expOVD$geneRLF=="CS2", gene]
hk1 <- normHK[expOVD$geneRLF=="CS1", gene]
hk2 <- normHK[expOVD$geneRLF=="CS2", gene]
par(mfrow = c(2, 1), oma = c(.5, .5, 2, .5), mar = c(5.1,4.1,1.1,2.1))
CCplot(hk1, pc1, Ptype = "scatter", xlabel = "HK1", ylabel = "PC1", subtitle = paste( letters[1],":","CodeSet 1"))
CCplot(hk2, pc2, Ptype = "scatter", xlabel = "HK2", ylabel =  "PC2", subtitle =paste( letters[2], ":","CodeSet 2"))
title(gene, outer = T)

## ----MABackground, fig.height = 8, fig.width = 6, fig.cap="Ovarian Cancer data comparisons of two normalization procedures of the data of the gene `r gene`. The method that performs only housekeeping normalization is denoted by HK1 and HK2 in CodeSet 1 and 2 respectively. On the other hand, the method that normalizes housekeeping genes and performs background correction, is denoted by NC1 and NC2 when ran on CodeSet 1 and CodeSet 2 respectively."----
nc1 <- normNC[expOVD$geneRLF=="CS1", gene]
nc2 <- normNC[expOVD$geneRLF=="CS2", gene]
hk1 <- normHK[expOVD$geneRLF=="CS1", gene]
hk2 <- normHK[expOVD$geneRLF=="CS2", gene]
par(mfrow = c(2, 1), oma = c(.5, .5, 2, .5), mar = c(5.1,4.1,1.1,2.1))
CCplot(hk1, nc1, Ptype = "scatter", xlabel = "HK1", ylabel = "NC1", subtitle = paste( letters[1],":","CodeSet 1"))
CCplot(hk2, nc2, Ptype = "scatter", xlabel = "HK2", ylabel =  "NC2", subtitle = paste( letters[2],":","CodeSet 2"))
title(gene, outer = T)


## ----MAmeansdPlot,message=FALSE,echo=FALSE,fig.height=8,fig.width=6, fig.cap="In both CodeSets we  see that where the negative control genes are in the lower limit of detection. The housekeeping genes have constant varianceand it's low relative to other genes."----
Plot.NanoStringNorm(x = MAHK.d,label.best.guess = FALSE, plot.type = 'mean.sd',title=FALSE);
title("OC Clinical Samples", outer=T)

## ----fig.height=8,fig.width=6--------------------------------------------
par(mfrow = c(2, 1), oma = c(.5, .5, 2, .5), mar = c(5.1,4.1,1.1,2.1))
ovd.n <- t(HKnorm(ovd.r)[,-(1:3)])
ovd.1 <- ovd.n[expOVD$geneRLF=="CS1",]
ovd.2 <- ovd.n[expOVD$geneRLF=="CS2",]

CCplot(hk1, ovd.1[,gene], Ptype = "scatter", xlabel = "Normalized with scaling", ylabel = "Normalized without scaling", subtitle = paste( letters[1],":","CodeSet 1"),xrange=c(-10,15),yrange=c(-10,15))

CCplot(hk2, ovd.2[, gene], Ptype = "scatter", xlabel = "Normalized with scaling", ylabel =  "Normalized without scaling", subtitle = paste( letters[2],":","CodeSet 2"),xrange=c(-10,15),yrange=c(-10,15))
title(gene, outer = T)


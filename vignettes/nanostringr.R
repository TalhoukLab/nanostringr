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

## ----message = FALSE, echo = FALSE, warning = FALSE----------------------
library(CHL26predictor)
library(dplyr)
library(nanostringr)
getNum <- function(str.vect){
  sapply(strsplit(str.vect,"[_]"),"[[",2)
}


## ----message = FALSE, echo = FALSE, warning = FALSE----------------------
# Normalize to HK 
hld.n <- HKnorm(hld.r)
expHLD <- subset(expQC,HLD=="Yes")

hld1 <- hld.n[,grep("HL1",colnames(hld.n))]
exp.hld1 <- subset(expHLD,geneRLF=="HL1")

hld2 <- hld.n[,grep("HL2",colnames(hld.n))]
exp.hld2 <- subset(expHLD,geneRLF=="HL2")

CHL26.HL1.exprs=hld.n[rownames(hld.n)%in%CHL26.model.coef.df$geneName,grep("HL1",colnames(hld.n))]+log(1000,2)
CHL26.HL2.exprs=hld.n[rownames(hld.n)%in%CHL26.model.coef.df$geneName,grep("HL2",colnames(hld.n))]+log(1000,2)
risk.thres <- 0.6235

scores.df1 <- get_CHL26_scores(as.matrix(CHL26.HL1.exprs))
scores.df2 <- get_CHL26_scores(as.matrix(CHL26.HL2.exprs))


scores.risk.df1 <- scores.df1 %>%
  mutate(riskClass = ifelse(score >= risk.thres, "High", "Low"))
scores.risk.df2 <- scores.df2 %>%
  mutate(riskClass = ifelse(score >= risk.thres, "High", "Low"))
tabRisk=table(scores.risk.df1$riskClass,scores.risk.df2$riskClass)
ind.mis <- which(scores.risk.df1$riskClass!=scores.risk.df2$riskClass)
n=(dim(CHL26.HL1.exprs)[2])
mis=n-sum(diag(tabRisk))


## ----RawScores, message = FALSE, echo = FALSE, warning = FALSE-----------
CCplot(scores.risk.df1$score,scores.risk.df2$score,Ptype = "scatter",xrange = range(scores.df1$score), yrange = range(scores.df2$score), xlabel = "HL1", ylabel = "HL2", subtitle = "Scores without Batch Adjustment")
abline(h = risk.thres, col=2, lty = 1)
text (0.7,0.58,paste("Misclassified:",mis),col="red")
points(scores.risk.df1$score[ind.mis],scores.risk.df2$score[ind.mis],col="red")

## ----BootStrapped, message = FALSE, echo = FALSE, warning = FALSE--------
set.seed(40)
r=3 #Number of reference samples
nB=5000# of bootstrap samples

nth=1000
misCount=rep(0,nB)
ind.mis.count=rep(0,nB)
Cmetrix=matrix(0, nrow = nB, ncol = 3)
for(i in 1:nB){
choice.refs <- exp.hld1$sampleID[sample((1:dim(exp.hld1)[1]), r, replace = F)] 
DSR1 <- t(hld1[,choice.refs]+log2(1000))
DSR2 <- t(hld2[,paste("HL2",getNum(choice.refs),sep="_")]+log2(1000))
DSY <- t(hld2[,!colnames(hld2)%in%paste("HL2",getNum(choice.refs),sep="_")]+log2(1000))

DSS2.r <- t(refMethod(DSY,DSR1, DSR2))

CHL26.HL2.r.SS.exprs=DSS2.r[rownames(DSS2.r)%in%CHL26.model.coef.df$geneName,]
scores.ss.df2.r <- get_CHL26_scores(as.matrix(CHL26.HL2.r.SS.exprs))
scores.risk.ss.df2.r <- scores.ss.df2.r %>%
  mutate(riskClass = ifelse(score >= risk.thres, "High", "Low"))
ndx2.r <- substring(scores.df1$sampleID,5) %in% substring(scores.ss.df2.r$sampleID,5)
misCount[i] <- (n-r)-sum(diag(table(scores.risk.df1$riskClass[ndx2.r],scores.risk.ss.df2.r$riskClass)))
ind.mis.t <- scores.risk.ss.df2.r$sampleID[which(scores.risk.df1$riskClass[ndx2.r]!=scores.risk.ss.df2.r$riskClass)]

ind.mis.count[i] <- (sum(ind.mis.t%in%c("HL2_30","HL2_32"))==misCount[i])

Cmetrix[i,]=CCplot(scores.df1$score[ndx2.r],scores.risk.ss.df2.r$score, metrics = TRUE)

if(i %% nth == 0){
CCplot(scores.df1$score[ndx2.r],scores.risk.ss.df2.r$score,Ptype = "scatter", xrange = range(scores.df1$score), yrange = range(scores.df2$score), xlabel = "HL1", ylabel = "HL2")
abline(h = risk.thres, col=2, lty = 1)
text(0.7,0.6,paste("Misclassified",misCount[i]))
}
}
TotalMisCount <- table(misCount)/nB
Accuracy <- Cmetrix[,2]
table(Accuracy)
if (length(TotalMisCount) > 3) {print("more than 3 misclassifications observed")}


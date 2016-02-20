# Introduction to using the nanostringr package
Aline Talhouk  








```r
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
```


# Metrics for Quality Assurance


## Fields Of View (FOV) 
 

<div class="figure">
<img src="nanostringr_files/figure-html/perFOVPlot-1.png" alt="Samples that failed imaging QC based on percent fields of view (FOV) counted across cohorts"  />
<p class="caption">Samples that failed imaging QC based on percent fields of view (FOV) counted across cohorts</p>
</div>


## Postive Controls
    
### Linearity of Positive Controls

<div class="figure">
<img src="nanostringr_files/figure-html/linPCPlot-1.png" alt="Plot of $R^2$ of postive control probles from samples across all cohorts."  />
<p class="caption">Plot of $R^2$ of postive control probles from samples across all cohorts.</p>
</div>


## Signal to Noise Ratio (SNR)

### Level of Housekeeping Genes

<div class="figure">
<img src="nanostringr_files/figure-html/averageHKPlot-1.png" alt="Average log expression of Housekeeping genes by Cohort."  />
<p class="caption">Average log expression of Housekeeping genes by Cohort.</p>
</div>


### Limit of Detection (LOD)

<div class="figure">
<img src="nanostringr_files/figure-html/lodPlot-1.png" alt="Limit of detection by cohort."  />
<p class="caption">Limit of detection by cohort.</p>
</div>


<div class="figure">
<img src="nanostringr_files/figure-html/pergdPlot-1.png" alt="Percent genes of total (excluding controls) detected above the limit of detection"  />
<p class="caption">Percent genes of total (excluding controls) detected above the limit of detection</p>
</div>



<div class="figure">
<img src="nanostringr_files/figure-html/snPlot-1.png" alt="Signal to Noise versus % Gene Detected by cohort"  />
<p class="caption">Signal to Noise versus % Gene Detected by cohort</p>
</div>

<div class="figure">
<img src="nanostringr_files/figure-html/snZoom-1.png" alt="Signal to Noise versus % Gene Detected by cohort zoomed in to the area of possible failures"  />
<p class="caption">Signal to Noise versus % Gene Detected by cohort zoomed in to the area of possible failures</p>
</div>


## Binding Density (BD)

<div class="figure">
<img src="nanostringr_files/figure-html/bdPlot-1.png" alt="Binding density by cohort. Samples outside of the dashed lines are flagged for having failed the binding density QC measure"  />
<p class="caption">Binding density by cohort. Samples outside of the dashed lines are flagged for having failed the binding density QC measure</p>
</div>


![](nanostringr_files/figure-html/lodbd-1.png)![](nanostringr_files/figure-html/lodbd-2.png)





# Normalization for Single-patient Samples
##  Distribution of HK genes
![](nanostringr_files/figure-html/unnamed-chunk-4-1.png)


## Manufacturer Suggested Normalization  



![](nanostringr_files/figure-html/PCNorm-1.png)



#### Background correction is not necessary


<div class="figure">
<img src="nanostringr_files/figure-html/MABackground-1.png" alt="Ovarian Cancer data comparisons of two normalization procedures of the data of the gene `r gene`. The method that performs only housekeeping normalization is denoted by HK1 and HK2 in CodeSet 1 and 2 respectively. On the other hand, the method that normalizes housekeeping genes and performs background correction, is denoted by NC1 and NC2 when ran on CodeSet 1 and CodeSet 2 respectively."  />
<p class="caption">Ovarian Cancer data comparisons of two normalization procedures of the data of the gene `r gene`. The method that performs only housekeeping normalization is denoted by HK1 and HK2 in CodeSet 1 and 2 respectively. On the other hand, the method that normalizes housekeeping genes and performs background correction, is denoted by NC1 and NC2 when ran on CodeSet 1 and CodeSet 2 respectively.</p>
</div>



#### Choice of Housekeeping Genes

<div class="figure">
<img src="nanostringr_files/figure-html/MAmeansdPlot-1.png" alt="In both CodeSets we  see that where the negative control genes are in the lower limit of detection. The housekeeping genes have constant varianceand it's low relative to other genes."  />
<p class="caption">In both CodeSets we  see that where the negative control genes are in the lower limit of detection. The housekeeping genes have constant varianceand it's low relative to other genes.</p>
</div>


#### Correction using Housekeeping Genes without scaling by overall mean.

![](nanostringr_files/figure-html/unnamed-chunk-5-1.png)


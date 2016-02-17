# Introduction to using the NanoString package
Derek Chiu  
`r Sys.Date()`  

## Introduction


```
## Warning: package 'ggplot2' was built under R version 3.2.3
```

We demonstrate how to use some of the functions in this package below.

## NanoStringQC


```r
data(rawOVCA2, rawPROT, rawOTTA, annot)

# Codeset 1, 2, 3 and annotations
cs1 <- rawOVCA2
cs2 <- rawPROT
cs3 <- rawOTTA
exp0 <- annot
exp0$geneRLF <- as.character(factor(exp0$geneRLF,
                     labels = c("HL1", "HL2", "HL3", "HuRef", "CS3", "mini", "CS1", "CS2")))

# Compute NanoString QC
exp.CS1 <- NanoStringQC(cs1, exp0[exp0$geneRLF == "CS1", ],
                        plots = FALSE, detect = 50, ttl = "CodeSet 1")
exp.CS2 <- NanoStringQC(cs2, exp0[exp0$geneRLF == "CS2", ],
                        plots = FALSE, sn = 100, ttl = "CodeSet 2")
exp.CS3 <- NanoStringQC(cs3, exp0[exp0$geneRLF == "CS3", ],
                        plots = TRUE, detect = 50, sn = 100, ttl = "CodeSet 3")
```

![](NanoString_files/figure-html/codesets-1.png)\

We can also compute the quality control parameters *just* for the the 3 POOLs:


```r
pool1 <- exp0$File.Name[which(exp0$POOL1 == "Yes")]
pool2 <- exp0$File.Name[which(exp0$POOL2 == "Yes")]
pool3 <- exp0$File.Name[which(exp0$POOL3 == "Yes")]
cs3.names <- str_sub(names(cs3), start = 2)
cs2.names <- str_sub(names(cs2), start = 2)
cs3.pools <- cs3[, c(1:3, sort(match(c(pool1, pool2, pool3), cs3.names)))]
cs3.pool1 <- cs3[, sort(match(pool1, cs3.names))]
cs3.pool2 <- cs3[, sort(match(pool2, cs3.names))]
cs3.pool3 <- cs3[, sort(match(pool3, cs3.names))]
cs2.pools <- cs2[, c(1:3, sort(match(c(pool1, pool2, pool3), cs2.names)))]

# Compute NanoString QC
exp.CS3.pools <- NanoStringQC(cs3.pools, exp0[sort(match(c(pool1, pool2, pool3), cs3.names)) - 3, ],
                              plots = TRUE, detect = 50, sn = 100,
                              ttl = "CodeSet 3 in POOL 1-3")
```

![](NanoString_files/figure-html/pools-1.png)\

```r
exp.CS2.pools <- NanoStringQC(cs2.pools, exp0[sort(match(c(pool1, pool2, pool3), cs2.names)) - 3, ],
                              plots = TRUE, detect = 50, sn = 100,
                              ttl = "CodeSet 2 in POOL 1-3")
```

![](NanoString_files/figure-html/pools-2.png)\

More importantly, we want to verify the concordance between CodeSet2 and CodeSet3 in terms of the gene expression in common genes for each of POOL1, POOL2, and POOL3.


```r
cs3.pool1 <- cs3[, c(1:3, sort(match(pool1, cs3.names)))]
cs3.pool2 <- cs3[, c(1:3, sort(match(pool2, cs3.names)))]
cs3.pool3 <- cs3[, c(1:3, sort(match(pool3, cs3.names)))]

cs2.pool1 <- cs2[, c(1:3, sort(match(pool1, cs2.names)))]
cs2.pool2 <- cs2[, c(1:3, sort(match(pool2, cs2.names)))]
cs2.pool3 <- cs2[, c(1:3, sort(match(pool3, cs2.names)))]

cs3.pool1.c <- cs3.pool1[cs3.pool1$Name %in% intersect(cs3.pool1$Name, cs2.pool1$Name), ]
cs2.pool1.c <- cs2.pool1[cs2.pool1$Name %in% intersect(cs3.pool1$Name, cs2.pool1$Name), ]

cs3.pool2.c <- cs3.pool2[cs3.pool2$Name %in% intersect(cs3.pool2$Name, cs2.pool2$Name), ]
cs2.pool2.c <- cs2.pool2[cs2.pool2$Name %in% intersect(cs3.pool2$Name, cs2.pool2$Name), ]

cs3.pool3.c <- cs3.pool3[cs3.pool3$Name %in% intersect(cs3.pool3$Name, cs2.pool3$Name), ]
cs2.pool3.c <- cs2.pool3[cs2.pool3$Name %in% intersect(cs3.pool3$Name, cs2.pool3$Name), ]

pool1.dat <- merge(cs3.pool1.c, cs2.pool1.c, by = c("Name", "Code.Class")) %>% 
  tidyr::gather(CodeSets, Expr, which(startsWith(names(.), "X"))) %>% 
  mutate(CodeSets = ifelse(grepl("RNA", CodeSets), "CS2", "CS3")) %>% 
  arrange(Name) %>% 
  select(Name, CodeSets, Expr)

pool2.dat <- merge(cs3.pool2.c, cs2.pool2.c, by = c("Name", "Code.Class")) %>% 
  tidyr::gather(CodeSets, Expr, which(startsWith(names(.), "X"))) %>% 
  mutate(CodeSets = ifelse(grepl("RNA", CodeSets), "CS2", "CS3")) %>% 
  arrange(Name) %>% 
  select(Name, CodeSets, Expr)

pool3.dat <- merge(cs3.pool3.c, cs2.pool3.c, by = c("Name", "Code.Class")) %>% 
  tidyr::gather(CodeSets, Expr, which(startsWith(names(.), "X"))) %>% 
  mutate(CodeSets = ifelse(grepl("RNA", CodeSets), "CS2", "CS3")) %>% 
  arrange(Name) %>% 
  select(Name, CodeSets, Expr)

# Randomly select 60 genes from each pool
set.seed(2016)
ngenes <- 60
pool1.rand <- pool1.dat %>% 
  magrittr::extract(.$Name %in% sample(unique(.$Name), ngenes), )
pool2.rand <- pool2.dat %>% 
  magrittr::extract(.$Name %in% sample(unique(.$Name), ngenes), )
pool3.rand <- pool3.dat %>% 
  magrittr::extract(.$Name %in% sample(unique(.$Name), ngenes), )

# Create side-by-side boxplots
ggplot(pool1.rand, aes(x = CodeSets, y = Expr)) +
  geom_boxplot() +
  facet_wrap(~ Name, ncol = 6) +
  scale_y_continuous(trans = "log2")
#> Warning: Removed 1 rows containing non-finite values (stat_boxplot).
```

![](NanoString_files/figure-html/pool_concord-1.png)\

```r

ggplot(pool2.rand, aes(x = CodeSets, y = Expr)) +
  geom_boxplot() +
  facet_wrap(~ Name, ncol = 6) +
  scale_y_continuous(trans = "log2")
#> Warning: Removed 3 rows containing non-finite values (stat_boxplot).
```

![](NanoString_files/figure-html/pool_concord-2.png)\

```r

ggplot(pool3.rand, aes(x = CodeSets, y = Expr)) +
  geom_boxplot() +
  facet_wrap(~ Name, ncol = 6) +
  scale_y_continuous(trans = "log2")
#> Warning: Removed 1 rows containing non-finite values (stat_boxplot).
```

![](NanoString_files/figure-html/pool_concord-3.png)\

## ratioMethod


```r
set.seed(12)
A <- matrix(rnorm(120), ncol = 10)
B <- matrix(rnorm(80), ncol = 10)
C <- matrix(rnorm(50), ncol = 10)
pander::pandoc.table(ratioMethod(A, B, C))
```


------- -------- ------- ------- ------- -------- --------- -------- ------- -------
-0.4833  -1.029  -0.4359 -0.6705 0.5392  0.08516  -0.006878 -0.08182 0.5791  -1.612 

 2.574  -0.2378  0.3219  -1.425  -0.6848  0.673     0.24     0.8387  -0.5679  1.076 

0.04048 -0.4022  0.3902  -0.1471 -0.2349  0.5344   0.3192   -0.2352  -0.1122 -0.3695

0.07722 -0.9532  0.7205   -1.95  -0.3049 -0.1243   -0.2113   1.055   -0.8805 -0.9702

  -1     0.9391  0.7351  -0.8403 0.2646   0.5129   -0.8485  -0.6003   2.054  -0.967 

0.7249  0.09076  0.9514  -2.101   1.828   0.5254   0.03858   -1.022  0.1782  0.3044 

0.6819   0.2572   1.263  -0.3671 -1.243   1.633    0.9393    0.318    1.927  -1.625 

 0.369  -0.5431   2.661  -2.541  0.5424   -2.471   0.2446    -1.043  -0.4247 -1.797 

0.8908  -0.02611 0.0483  -1.254   0.347   0.6497   0.4268    0.8391  0.2339  -0.6548

 1.425   1.757   -0.4812 -0.4958 -1.507   0.8237   -0.8515  0.07082  0.9159  -0.532 

0.2195   0.7622  0.2169  -1.922  -0.4423 -0.8468   0.2204    0.9753  -0.1175 -0.1855

-0.2967 -0.5522  0.1042  -0.7552  0.122  -0.07106   1.587   -0.4007   0.511  0.7618 
------- -------- ------- ------- ------- -------- --------- -------- ------- -------

## HKnorm


```r
data(NanoString)
NanoString.mRNA[NanoString.mRNA$Name %in%
c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'), 'Code.Class'] <- 'Housekeeping'
out <- HKnorm(NanoString.mRNA)
pander::pandoc.table(head(out))
```


----------------------------------------------------------------------------
   &nbsp;      Code.Class   Name     Accession     HW1D.a   HW2B.b   HW3B.a 
------------- ------------ ------- -------------- -------- -------- --------
  **Ahrr**     Endogenous   Ahrr   NM_001024285.1  -9.062   -12.37   -12.21 

 **Aldh3a1**   Endogenous  Aldh3a1  NM_031972.1    -10.38   -9.783   -9.626 

  **Bbs2**     Endogenous   Bbs2    NM_053618.1    -7.683   -7.368   -7.567 

  **Ccbl1**    Endogenous   Ccbl1  NM_001013164.3  -4.876   -4.95    -5.123 

 **Cmkor1**    Endogenous  Cmkor1   NM_053352.1    -8.214   -8.12    -8.626 

 **Col18a1**   Endogenous  Col18a1  NM_053489.1    -2.04    -1.97    -2.735 
----------------------------------------------------------------------------

Table: Table continues below

 
-------------------------------------------------------------------------------
   &nbsp;      HW4C.b   HW5D.b   HW6A.b   WW41B.a   WW42.a   WW45A.b   WW46D.a 
------------- -------- -------- -------- --------- -------- --------- ---------
  **Ahrr**     -9.846   -10.1    -10.43    -7.73    -7.294   -8.525     -8.04  

 **Aldh3a1**   -4.898   -6.98    -6.52    -0.5884   0.503    -1.188    -2.409  

  **Bbs2**     -7.038   -7.288   -7.289   -6.122    -6.012   -6.747    -6.318  

  **Ccbl1**    -5.053   -4.697   -4.941   -2.905    -2.751   -3.561    -3.495  

 **Cmkor1**    -7.068   -7.925   -7.539    -7.23    -7.108   -7.474    -7.153  

 **Col18a1**   -1.953   -2.299   -1.646   -1.553    -1.61    -2.331    -2.359  
-------------------------------------------------------------------------------

Table: Table continues below

 
----------------------------------------------------------------------------
   &nbsp;      WW47A.b   WW48B.a   LE100C.a   LE101A.a   LE102.a   LE4D51.a 
------------- --------- --------- ---------- ---------- --------- ----------
  **Ahrr**     -8.019    -9.482     -9.051     -11.82    -9.461     -9.912  

 **Aldh3a1**   -1.714    -0.9331    -6.865     -5.284    -5.891     -6.461  

  **Bbs2**     -7.149    -7.023     -7.285     -7.823    -7.139     -7.295  

  **Ccbl1**    -3.517    -3.389     -5.073     -4.834    -4.189     -5.143  

 **Cmkor1**    -6.648    -6.867     -7.285     -7.238    -7.175     -6.609  

 **Col18a1**   -2.069    -2.467     -2.258     -1.933    -1.746     -1.758  
----------------------------------------------------------------------------

Table: Table continues below

 
--------------------------------------------------------------------------------
   &nbsp;      LE4D53.a   LE4D56.b   WW37.b   WW38.b   WW39.a   WW40.b   WW43.b 
------------- ---------- ---------- -------- -------- -------- -------- --------
  **Ahrr**      -9.409     -9.592    -7.742   -8.062   -8.224   -8.251   -8.259 

 **Aldh3a1**    -5.642     -4.11    0.01903  -0.6194   0.6325  0.09724  -0.6031 

  **Bbs2**      -6.824     -7.025    -5.967   -6.206   -5.818   -6.264   -6.585 

  **Ccbl1**     -4.484     -4.595    -2.813   -2.656   -2.787   -3.171   -2.557 

 **Cmkor1**     -6.824     -6.737    -8.235   -7.966   -7.26    -7.849   -7.225 

 **Col18a1**    -1.866     -1.477    -1.885   -2.095   -1.599   -1.82    -1.642 
--------------------------------------------------------------------------------

Table: Table continues below

 
---------------------------------
   &nbsp;      WW44.a   HW4D12.B 
------------- -------- ----------
  **Ahrr**       -8      -10.8   

 **Aldh3a1**  -0.9055    -10.8   

  **Bbs2**     -6.104    -7.418  

  **Ccbl1**    -2.651    -4.587  

 **Cmkor1**    -7.737    -7.949  

 **Col18a1**   -1.61     -1.508  
---------------------------------

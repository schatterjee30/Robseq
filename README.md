# `Robseq` <a style="position: relative; display: inline-block;"><img src='man/figures/RobseqLogo.png' align="right"  height="160"/></a>

[![Build
Status](https://img.shields.io/badge/build-ok-brightgreen)](%5Bhttps://example-link-to-build-status-page%5D(https://github.com/schatterjee30/Robseq/blob/main/README.md))
 
[![HitCount](http://hits.dwyl.com/schatterjee30/Robseq.svg)](http://hits.dwyl.com/schatterjee30/Robseq "Get hits on your repository!")
  [![GitHub
stars](https://img.shields.io/github/stars/schatterjee30/Robseq.svg?style=social&color=green&label=Stars&cacheBust=1)](https://github.com/schatterjee30/Robseq/stargazers)

# Overview

`Robseq` A Robust Statistical Model for Differential Gene Expression
Analysis in RNA-Seq Studies

<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/dearseq?color=blue)](https://www.r-pkg.org/pkg/dearseq) --
&#10;
# Robseq: A Robust Statistical Model for Differential Gene Expression Analysis in RNA-Seq Studies
&#10;<!-- badges: start -->
<!-- badges: end -->

## Robseq Model Pipeline

<figure>
<img src="Robseq%20Pipeline.png" alt="Pipeline" />
</figure>

## Installation

For the installation, the R package devtools is needed.

``` r
install.packages('devtools')
library(devtools)
```

Robseq has a couple of Bioconductor dependencies that need to be
installed before hand. We recommend to install first the dependencies
manually and then Robseq.

Bioconductor Libraries can be installed in the following manner:

``` r
install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
```

After installing the dependencies, Robseq can be installed by using
devtools.

``` r
devtools::install_github("schatterjee30/Robseq")
```

## Arguments

`Robseq` has six main arguments which the user needs to supply/define.

### The table below details the required arguments:

| Parameter   | Default  | Description                                                                                                                                                                                                                                                         |
|:------------|:--------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| features    |          | A dataframe with gene expression counts in the row and samples in the column.                                                                                                                                                                                       |
| metadata    |          | A dataframe with information on the subjects such as disease status, gender and etc.                                                                                                                                                                                |
| norm.method |   TMM    | The normalization method to be used. The user can choose from 5 different methods such as TMM, RLE, CPM, Upper quartile and Qauntile.                                                                                                                               |
| expVar      | Exposure | The name of the variable on which the differential expression will be evaluated. If the user provides no name then the metadata should have a column named as exposure which should have information on things such as disease status, treatment conditions or etc. |
| coVars      |   NULL   | The names of the covariates/confounders that needs to adjusted for in the differential expression analysis.                                                                                                                                                         |
| parallel    |  FALSE   | If true, the analysis will be performed on multiple cores with faster runtimes.                                                                                                                                                                                     |
| ncores      |    1     | The number of cores on which the analysis will be serially performed. The user needs to specify this only when parallel = TRUE.                                                                                                                                     |
| verbose     |  FALSE   | If true, it will print the progress report.                                                                                                                                                                                                                         |

## Values

`Robseq` outputs has 7 value arguments.

### The table below details the values returned by Robseq:

| Value   | Description                                                                                       |
|:--------|:--------------------------------------------------------------------------------------------------|
| Genes   | Gene that was analysed for differential expression between treatment conditions or disease status |
| Log2FC  | The estimated log2 fold change for the genes that were analysed                                   |
| SE      | The standard error for the genes that were analysed                                               |
| LCI     | The lower confidence interval for the genes that were analysed                                    |
| UCI     | The upper confidence interval for the genes that were analysed                                    |
| Pval    | The pvalue for the genes that were analysed                                                       |
| adjPval | The BH adjustd pvalue after correcting for multipe testing for the genes that were analysed       |

## Working Example

### Importing libraries

``` r
library(Robseq)
library(edgeR)
library(doParallel)
library(EnhancedVolcano)
```

### Loading Example Data

Loading Colon cancer data

``` r
load("~/Current Data Path/Colon Cancer.RData") 
```

Extracting Colon cancer gene expression data and metadata

``` r
features = data$counts
metadata = data$metadata
```

### Snapshot of data

A typical bulk RNA-seq gene expression count data frame looks something
like below. Note, the genes should be in the rows and the samples in
columns

``` r
features[1:5, 1:5]
```

    ##          GSM4731674 GSM4731675 GSM4731676 GSM4731677 GSM4731678
    ## TSPAN6          103         44         76        417        630
    ## TNMD              1          0          0          0         18
    ## DPM1             86         63         97        173        309
    ## SCYL3            32         40         77         56         97
    ## C1orf112         28         35         25         60         55

A typical metadata data frame looks something like below. Note, in our
pipeline the user must include a column labeled as “Exposure” which
should be the variable which will be used by Robseq in performing
differential expression analysis. If the user has a different label for
the treatment/condition or disease status variable then the user should
supply that name via “expVar” argument in Robseq

``` r
metadata[1:5, ]
```

    ##       Sample Exposure
    ## 1 GSM4731674    Tumor
    ## 2 GSM4731675    Tumor
    ## 3 GSM4731676    Tumor
    ## 4 GSM4731677    Tumor
    ## 5 GSM4731678    Tumor

### Preprocessing

A typical preprocessing step is to filter lowly abundant genes. We do so
using edgeR’s “filterByExpr” function

``` r
keep.exprs <- filterByExpr(features, group = as.factor(metadata$Exposure))

```

    ## [1] "18625 lowly expressed genes were filtered out"

``` r
features <- features[keep.exprs, ]
```

### Performing differential expression analysis using Robseq

To perform differential gene expression (DGE) analyses using Robseq
please use the code below. Note, if your metadata has only the
“Exposure” variable (such as treatment groups, conditions, disease
status, and etc.) set the “coVars” argument to “NULL”. However, in this
working example our metadata contained three covariates such as
“post-mortem interval (pmi)”, “rna integrity number (rin)” and “age at
death” which needed to be adjusted for in our DGE analysis. Therefore,
we have supplied this three variables to the “coVars” argument

``` r
fit <- Robseq::robust.dge(features = features,
                          metadata = metadata,
                          norm.method = "RLE",
                          expVar =  "Exposure",
                          coVars = NULL,
                          parallel = TRUE,
                          ncores = detectCores() - 2,
                          verbose = FALSE)
```

### Obtaining results from Robseq

After performing DGE analysis you can extract the results table in the
following manner

``` r
results <- fit$res
```

The results table from Robseq should look something like the following

``` r
results[1:5,]
```

    ##       Genes log2FC    SE      L.CI      U.CI         Pval      adjPval
    ## 6872  BEST4 -6.631 0.255 -6.131209 -7.130791 1.565312e-54 1.617267e-50
    ## 10982  ETV4  5.184 0.202  5.579913  4.788087 2.121424e-54 1.617267e-50
    ## 11710 OTOP3 -6.227 0.244 -5.748769 -6.705231 1.430501e-53 7.270284e-50
    ## 11725 OTOP2 -7.246 0.298 -6.661931 -7.830069 2.196497e-51 8.372498e-48
    ## 14866  SPIB -5.385 0.239 -4.916569 -5.853431 7.877754e-48 2.402242e-44

### Volcano plot to visualize DGE results

Volcano plot to visualize the DGE results obtained from Robseq

``` r
  EnhancedVolcano(results,
    lab = results$Genes,
    x = 'log2FC',
    y = 'adjPval')
```
![](/man/volcano.png)<!-- -->

### Citation

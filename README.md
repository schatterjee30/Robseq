# `Robseq` <a style="position: relative; display: inline-block;"><img src='man/figures/RobseqLogo.png' align="right"  height="160"/></a>

   


[![Build Status](https://img.shields.io/badge/build-ok-brightgreen)]([https://example-link-to-build-status-page](https://github.com/schatterjee30/Robseq/blob/main/README.md))  &nbsp;
[![HitCount](http://hits.dwyl.com/schatterjee30/Robseq.svg)](http://hits.dwyl.com/schatterjee30/Robseq "Get hits on your repository!") 
 &nbsp;
[![GitHub stars](https://img.shields.io/github/stars/schatterjee30/Robseq.svg?style=social&color=green&label=Stars&cacheBust=1)](https://github.com/schatterjee30/Robseq/stargazers)


## Overview

`Robseq` A Robust Statistical Model for Differential Gene Expression Analysis in RNA-Seq Studies


<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/dearseq?color=blue)](https://www.r-pkg.org/pkg/dearseq) --


# Robseq: A Robust Statistical Model for Differential Gene Expression Analysis in RNA-Seq Studies

<!-- badges: start -->

<!-- badges: end -->

## Robseq Model Pipeline
![Pipeline](Pipeline%20Image.png)

## Installation

For the installation, the R package devtools is needed.
``` r
install.packages('devtools')
library(devtools)
```
Robseq has a couple of Bioconductor dependencies that need to be installed before hand. We recommend to install first the dependencies manually and then Robseq.

Bioconductor Libraries can be installed in the following manner:
``` r
install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
```
After installing the dependencies, Robseq can be installed by using devtools.
``` r
devtools::install_github("schatterjee30/Robseq")
```

## Arguments

`Robseq` has six main arguments which the user needs to supply/define.

### The table below details the required arguments:

| Parameter     | Default  | Description                                                                                                          |
|:--------------|:--------:|:---------------------------------------------------------------------------------------------------------------------|
| features |  | A dataframe with gene expression counts in the row and samples in the column.
| metadata |  | A dataframe with information on the subjects such as disease status, gender and etc.       
| norm.method | TMM | The normalization method to be used. The user can choose from 5 different methods such as TMM, RLE, CPM, Upper quartile and Qauntile.  
| expVar | Exposure | The name of the variable on which the differential expression will be evaluated. If the user provides no name then the metadata should have a column named as exposure which should have information on things such as disease status, treatment conditions or etc.
| coVars | NULL | The names of the covariates/confounders that needs to adjusted for in the differential expression analysis.
| parallel | FALSE | If true, the analysis will be performed on multiple cores with faster runtimes.
| ncores | 1 | The number of cores on which the analysis will be serially performed. The user needs to specify this only when parallel = TRUE.

## Values

`Robseq` outputs has 7 value arguments.

### The table below details the values returned by Robseq:

| Value     | Description                                                                                                              |
|:--------------|:---------------------------------------------------------------------------------------------------------------------|
| Genes | Gene that was analysed for differential expression between treatment conditions or disease status
| Log2FC | The estimated log2 fold change for the genes that were analysed
| SE | The standard error for the genes that were analysed
| LCI | The lower confidence interval for the genes that were analysed
| UCI | The upper confidence interval for the genes that were analysed
| Pval | The pvalue for the genes that were analysed
| adjPval | The BH adjustd pvalue after correcting for multipe testing for the genes that were analysed

## Example

This is a basic example of the functions in the package.

``` r
library(Robseq)
## basic example code

```

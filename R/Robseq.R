#### Load libraries ####

library(edgeR)
library(doSNOW)
library(dfadjust)
library(lattice)
library(DESeq2)

#### Normalizing using Quantile ####

quan_norm <- function(features,metadata){
  norm.y <- normalize.quantiles(as.matrix(features))
  norm.y <- data.frame(norm.y)
  names(norm.y) <- names(features)
  rownames(norm.y) <- rownames(features)
  return(norm.y)
}

#### Normalizing using Upper Quartile ####

uqrt_norm <- function(features,metadata){
  quant.exp <- apply(as.matrix(features), 2, function(x){quantile(x[x > 0], 0.75)})
  norm.y <- data.frame(t(t(as.matrix(features))/quant.exp))
  return(norm.y)
}

#### Normalizing using TMM (edgeR) ####

tmm_norm <- function(features, metadata){
  norm.y <- DGEList(features)
  norm.y <- edgeR::calcNormFactors(norm.y, method = "TMM")
  norm.y <- as.data.frame(cpm(norm.y, log = FALSE))
  return(norm.y)
}

#### Normalizing using Geometric Means (DESeq2) ####

rle_norm <- function(features, metadata){
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep = ''))
  x <- suppressMessages(DESeqDataSetFromMatrix(countData = as.matrix(features),
                                               colData = metadata,
                                               design = formula))
  gm_mean <- function(x, na.rm = TRUE){
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
  }
  geoMeans <- apply(counts(x), 1, gm_mean)
  s <- estimateSizeFactors(x,geoMeans = geoMeans)
  s <- s$sizeFactor
  norm.y <- data.frame(t(apply(features, 1, function(x)x/s)))
  return(norm.y)
}

#### Per-Gene Modeling Function ####

perGene.mod = function(expr, formula, regData, expVar){

  ### Fitting RLM ###

  regData <- data.frame(expr, regData)
  fit.rlm <- tryCatch({
    fit.rlm <- suppressWarnings(rlm(formula, data = regData))
  }, error=function(err){
    fit.rlm <- suppressWarnings(lm(formula, data = regData))
    return(fit.rlm)
  })

  ### Estimating BM adjusted SE ###

  fit.se <- tryCatch({
    fit.se <- dfadjustSE(fit.rlm)
  }, error = function(err){
    fit.lm <- suppressWarnings(lm(formula, data = regData))
    fit.se <- dfadjustSE(fit.lm)
    return(fit.se)
  })

  ### Collecting Outputs ###

  if(is.list(fit.se)){
    rDF <- as.data.frame(cbind(fit.se$coefficients,
                               "p-value" = 2*stats::pt(-abs(fit.se$coefficients[, "Estimate"] /
                                                              fit.se$coefficients[, "HC2 se"]),
                                                       df = fit.se$coefficients[, "df"])))
    log2FC <- round(rDF[grep(expVar, rownames(rDF)), 'Estimate'], 3)
    se <- round(rDF[grep(expVar, rownames(rDF)), 'HC2 se'], 3)
    U.CI <- log2FC + qnorm(.025) * se
    L.CI <- log2FC - qnorm(.025) * se
    pval <- rDF[grep(expVar, rownames(rDF)), 'p-value']
  }else{
    log2FC <- NA
    se <- NA
    U.CI <- NA
    L.CI <- NA
    pval <- NA
  }

  ### Summarizing Outputs ###

  est.df <- data.frame(log2FC = log2FC,
                       SE = se,
                       L.CI = L.CI,
                       U.CI = U.CI,
                       Pval = pval)
  return(est.df)
}

#### Main Function ####

Robseq <- function(features,
                       metadata,
                       norm.method = 'tmm',
                       expVar = 'Exposure',
                       coVars = NULL,
                       parallel = FALSE,
                       ncores = 1){

  #### Creating Regression Pre-requisites ####

  start.time <- Sys.time()
  if(is.null(coVars)){
    regData <- metadata[, c(expVar), drop = FALSE]
  }else{
    regData <- metadata[, c(expVar, coVars)]
  }
  regData[sapply(regData, is.character)] <- lapply(regData[sapply(regData, is.character)], as.factor)
  formula <- as.formula(paste("expr ~ ", paste(colnames(regData), collapse = "+")))

  #### Normalizing Expression counts ####

  if(norm.method == 'tmm'){
    norm.y = suppressMessages(tmm_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else if(norm.method == 'rle'){
    norm.y = suppressMessages(rle_norm(features, metadata))
    norm.y = log2(norm.y + 0.5)
  }else{
    norm.y = data.frame(cpm(features, log = TRUE, prior.count = 1))
  }

  #### Apply model ####

  if(parallel){
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    packages <- c('MASS', 'dfadjust')
    exports <- c('perGene.mod', 'expVar', 'coVars')
    pb <- txtProgressBar(max = nrow(norm.y), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    res <- foreach(j = 1:nrow(norm.y), .combine = rbind,
                   .packages = packages, .options.snow = opts,
                   .export = exports) %dopar% {
                     expr <- as.numeric(norm.y[j, ])
                     tmpfit <- perGene.mod(expr = expr,
                                           formula = formula,
                                           regData = regData,
                                           expVar = expVar)
                     return(tmpfit)
                   }
    close(pb)
    stopCluster(cl)
  }else{
    res <- pbapply(norm.y, 1, function(x) perGene.mod(expr = x,
                                                      formula = formula,
                                                      regData = regData,
                                                      expVar = expVar))
    res <- do.call('rbind', res)
  }

  #### Summarizing Output ####

  Genes <- rownames(norm.y)
  output <- data.frame(Genes = Genes, res)
  output$adjPval <- p.adjust(output$Pval, method = 'BH')
  output <- output[order(output$adjPval), ]
  stop.time <- Sys.time()
  Time.min = round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'Robseq',
                res = output,
                features = features,
                Time.min = Time.min)
  return(finres)
}
